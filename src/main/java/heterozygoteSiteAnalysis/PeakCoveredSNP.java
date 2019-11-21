package heterozygoteSiteAnalysis;

import GTFComponent.GeneExonIntervalTree;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.Arrays;
import java.util.HashMap;

public class PeakCoveredSNP {
    private File gtfFile, vcfFile, peakCallingRes, outputFile;
    private int readsInfimum;
    HashMap<String, HashMap<String, String>> peakGene, geneNames;
    private Logger logger;

    /**
     * Constructor
     * @param gtfFile GTF format file
     * @param vcfRecordFile vcf record file
     * @param peakCallingRes m6a peak calling result bed file
     * @param outputFile output file path
     * @param readsInfimum minimum reads coverage for a SNV site under m6A signal
     * @param logger Logger instance
     */
    public PeakCoveredSNP(String gtfFile, String vcfRecordFile, String peakCallingRes,
                          String outputFile, int readsInfimum, Logger logger) {
        this.logger = logger;
        this.gtfFile = new File(gtfFile);
        if (!this.gtfFile.exists() || !this.gtfFile.isFile()) {
            this.logger.error("invalid GTF file");
            System.exit(2);
        }
        this.vcfFile = new File(vcfRecordFile);
        if (!vcfFile.exists() || !this.vcfFile.isFile()) {
            this.logger.error("vcf file not exists");
            System.exit(2);
        }
        this.peakCallingRes = new File(peakCallingRes);
        if (!this.peakCallingRes.exists() || !this.peakCallingRes.isFile()) {
            this.logger.error("peak calling result bed file not exists");
            System.exit(2);
        }
        this.readsInfimum = readsInfimum;
        this.outputFile = new File(outputFile);
    }

    /**
     * filter m6A peaks which covered SNV site and output into file
     */
    public void filterSNPAndPeak() {
        // chr -> +/- -> m6AIntervalTree
        HashMap<String, HashMap<String, IntervalTree>> m6aTreeMap = this.getM6aPeakTree();
        // chr -> geneId -> exonIntervalTree
        HashMap<String, HashMap<String, IntervalTree>> geneExonTree = this.generateExonTree();
        // chr -> peakRange -> geneId
        this.peakMatchGeneId();
        // chr -> geneId -> geneName
        this.parseGTFFile();

        BufferedReader bfr = null;
        BufferedWriter bfw = null;
        try {
            bfr = new BufferedReader(
                    new InputStreamReader(new FileInputStream(this.vcfFile))
            );
            bfw = new BufferedWriter(
                    new OutputStreamWriter(new FileOutputStream(this.outputFile))
            );

            String line = "", writeOut;
            String[] info;
            String chrNum, refNc, altNc, majorNc, minorNc;
            int allele1Count, allele2Count, majorAlleleCount, minorAlleleCount;
            boolean specialMutation;
            int[] readsCount;
            int position;
            bfw.write("#chr\tgeneId\tgeneName\tpeakStrand\tposition\tpeakStart\tpeakEnd\tmajorAllele\tminorAllele\tmajorAlleleCount\tminorAlleleCount\n");
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;

                    info = line.split("\t");
                    chrNum = info[0];
                    refNc = info[3];
                    altNc = info[4];

                    HashMap<String, IntervalTree> chrTree = m6aTreeMap.getOrDefault(chrNum, null);
                    HashMap<String, IntervalTree> chrExon = geneExonTree.getOrDefault(chrNum, null);
                    HashMap<String, String> chrPeaks = peakGene.getOrDefault(chrNum, null);
                    if (chrTree == null || chrExon == null || chrPeaks == null)
                        continue;
                    position = Integer.parseInt(info[1]);
                    // whether SNV site contains multiple possible mutations
                    specialMutation = altNc.length() > 1 && altNc.contains(",");

                    if (info.length >= 10) {    // may contains FORMAT and sample
                        readsCount = this.parseAd(info[8], info[9]);   // parse AD(allele depth) from VCF record
                        if (readsCount != null) {
                            if (!specialMutation) {
                                allele1Count = readsCount[0];
                                allele2Count = readsCount[1];
                                if (allele1Count >= allele2Count) {
                                    majorAlleleCount = allele1Count;
                                    minorAlleleCount = allele2Count;
                                    majorNc = refNc;
                                    minorNc = altNc;
                                } else {
                                    majorAlleleCount = allele2Count;
                                    minorAlleleCount = allele1Count;
                                    majorNc = altNc;
                                    minorNc = refNc;
                                }
                            } else {
                                String[] possibleMutations = altNc.split(",");
                                int[] originReads = Arrays.copyOf(readsCount, readsCount.length);
                                Arrays.sort(readsCount);
                                majorAlleleCount = readsCount[readsCount.length - 1];
                                minorAlleleCount = readsCount[readsCount.length - 2];
                                int[] indexRes = findIndex(originReads, majorAlleleCount, minorAlleleCount);
                                int majorAlleleIdx = indexRes[0];
                                int minorAlleleIdx = indexRes[1];
                                majorNc = (majorAlleleIdx == 0)? refNc : possibleMutations[majorAlleleIdx - 1];
                                minorNc = (minorAlleleIdx == 0)? refNc : possibleMutations[minorAlleleIdx - 1];
                            }
                        } else {    // if FORMAT contains no AD information, try to parse DP4 of INFO column
                            readsCount = this.parseDp4(info[7]);
                            if (readsCount == null) {
                                System.out.println("invalid VCF record: " + line );
                                System.exit(2);
                            }
                            allele1Count = readsCount[0] + readsCount[1];
                            allele2Count = readsCount[2] + readsCount[3];
                            majorAlleleCount = (allele1Count >= allele2Count)? allele1Count: allele2Count;
                            minorAlleleCount = (allele1Count < allele2Count)? allele1Count: allele2Count;
                            majorNc = (allele1Count >= allele2Count)? refNc : altNc.split(",")[0];
                            minorNc = (allele1Count < allele2Count)? refNc: altNc.split(",")[0];
                        }
                    } else {  // contains INFO but no INFO and Sample
                        readsCount = this.parseDp4(info[7]);
                        if (readsCount == null) {
                            System.out.println("invalid VCF record: " + line );
                            System.exit(2);
                        }
                        allele1Count = readsCount[0] + readsCount[1];
                        allele2Count = readsCount[2] + readsCount[3];
                        majorAlleleCount = (allele1Count >= allele2Count)? allele1Count: allele2Count;
                        minorAlleleCount = (allele1Count < allele2Count)? allele1Count: allele2Count;
                        majorNc = (allele1Count >= allele2Count)? refNc : altNc.split(",")[0];
                        minorNc = (allele1Count < allele2Count)? refNc: altNc.split(",")[0];
                    }

                    // filter low coverage SNP sites
                    if (majorAlleleCount < this.readsInfimum)
                        continue;

                    IntervalTree posStrandTree = chrTree.getOrDefault("+", null);
                    IntervalTree negStrandTree = chrTree.getOrDefault("-", null);


                    if (posStrandTree != null) {
                        IntervalTreeNode posStrandSearchResult = posStrandTree.search(posStrandTree.root, position);
                        if (posStrandSearchResult == null)
                            continue;
                        String peakStart = Integer.toString(posStrandSearchResult.peakStart);
                        String peakEnd = Integer.toString(posStrandSearchResult.peakEnd);
                        String geneId = chrPeaks.get(String.join("->", new String[] {peakStart, peakEnd}));
                        if (!this.ifInExon(chrExon, geneId, position))
                            continue;
                        String geneName = this.geneNames.get(chrNum).getOrDefault(geneId, "unknown");
                        String[] newLine = new String[]{chrNum, geneId, geneName, "+", info[1], peakStart,peakEnd,
                                                        majorNc, minorNc, String.valueOf(majorAlleleCount),
                                                        String.valueOf(minorAlleleCount)};
                        writeOut = String.join("\t", newLine);
                        bfw.write(writeOut);
                        bfw.newLine();
                    }
                    if (negStrandTree != null) {
                        IntervalTreeNode negStrandSearchResult = negStrandTree.search(negStrandTree.root, position);
                        if (negStrandSearchResult == null)
                            continue;
                        String peakStart = Integer.toString(negStrandSearchResult.peakStart);
                        String peakEnd = Integer.toString(negStrandSearchResult.peakEnd);
                        String geneId = chrPeaks.get(String.join("->", new String[] {peakStart, peakEnd}));
                        if (!this.ifInExon(chrExon, geneId, position))
                            continue;
                        String geneName = this.geneNames.get(chrNum).getOrDefault(geneId, "unknown");
                        String[] newLine = new String[]{chrNum, geneId, geneName, "-", info[1], peakStart, peakEnd,
                                                        majorNc, minorNc, String.valueOf(majorAlleleCount),
                                                        String.valueOf(minorAlleleCount)};
                        writeOut = String.join("\t", newLine);
                        bfw.write(writeOut);
                        bfw.newLine();
                    }
                }
            }
            m6aTreeMap = null;
            geneExonTree = null;
            this.geneNames = null;
            this.peakGene = null;
            bfw.flush();
        } catch (IOException ie) {
            this.logger.error("load peak calling result information failed.");
            this.logger.error(ie.getMessage());
        } finally {
            if (bfr != null) {
                try {
                    bfr.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
            if (bfw != null) {
                try {
                    bfw.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    private void peakMatchGeneId() {
        BufferedReader bfr = null;
        this.peakGene = new HashMap<>();
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(this.peakCallingRes)));
            String line = "", chrNum, peakStart, peakEnd, geneId, label;
            String[] info;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    chrNum = info[0];
                    peakStart = info[1];
                    peakEnd = info[2];
                    geneId = info[3];
                    HashMap<String, String> chrPeaks = peakGene.getOrDefault(chrNum, new HashMap<>());
                    label = String.join("->", new String[] {peakStart, peakEnd});
                    chrPeaks.put(label, geneId);
                    peakGene.put(chrNum, chrPeaks);
                }
                info = null;
            }
        } catch (IOException ie) {
            ie.printStackTrace();
            System.exit(2);
        } finally {
            if (bfr != null) {
                try {
                    bfr.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    private void parseGTFFile() {
        BufferedReader bfr = null;
        this.geneNames = new HashMap<>();
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(this.gtfFile)));
            String line = "", chrNum, geneId, geneName;
            String[] info, geneInfo;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    if (!info[2].equals("gene")) {
                        info = null;
                        continue;
                    }
                    chrNum = info[0];
                    geneInfo = this.getGeneInfo(info[8]);
                    geneId = geneInfo[0];
                    geneName = geneInfo[1];
                    HashMap<String, String> chrGenes = this.geneNames.getOrDefault(chrNum, new HashMap<>());
                    chrGenes.put(geneId, geneName);
                    this.geneNames.put(chrNum, chrGenes);
                }
                info = null;
            }
        } catch (IOException ie) {
            ie.printStackTrace();
            System.exit(2);
        } finally {
            if (bfr != null) {
                try {
                    bfr.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    private HashMap<String, HashMap<String, IntervalTree>> generateExonTree() {
        GeneExonIntervalTree geit = new GeneExonIntervalTree(this.gtfFile.getAbsolutePath());
        geit.generateExonTree();

        return geit.getGeneExonIntervalTree();
    }

    /**
     * build m6A signal interval tree using peak calling result
     * @return m6A interval tree for each strand on each chromosome [chrNum: [Strand: interval tree]]
     */
    private HashMap<String, HashMap<String, IntervalTree>> getM6aPeakTree() {
        PeakIntervalTree pit = new PeakIntervalTree(this.peakCallingRes.getAbsolutePath(), this.logger);
        HashMap<String, HashMap<String, IntervalTree>> treeMap = pit.getPeakTrees();
        pit = null;

        return treeMap;
    }

    /**
     * get reference and alternative nucleotide from VCF format file
     * @param format FORMAT information in VCF file
     * @param sampleInfo sample information in VCF file
     * @return int[] {ref, ref_reverse, alt, alt_reverse}
     */
    private int[] parseAd(String format, String sampleInfo) {
        String[] formatContent = format.split(":");
        String[] sampleData = sampleInfo.split(":");
        int idx = -1;
        for (int i=0; i<formatContent.length; i++) {
            if (formatContent[i].equals("AD")) {
                idx = i;
                break;
            }
        }
        if (idx < 0)
            return null;
        String[] readsRecord = sampleData[idx].split(",");
        int[] readsCount = new int[readsRecord.length];
        for (int i=0; i<readsRecord.length; i++) {
            readsCount[i] = Integer.valueOf(readsRecord[i]);
        }
        formatContent = null;
        sampleData = null;
        readsRecord = null;
        return readsCount;
    }

    private int[] findIndex(int[] readsCount, int majorCount, int minorCount) {
        int majorIdx = -1, minorIdx = -1;
        for (int i=0; i<readsCount.length; i++) {
            if (readsCount[i] == majorCount && majorIdx < 0)
                majorIdx = i;
            if (readsCount[i] == minorCount && minorIdx < 0 && majorIdx != i)
                minorIdx = i;
        }
        return new int[] {majorIdx, minorIdx};
    }

    /**
     * get reference and alternative nucleotide from VCF format file
     * @param dp4String INFO column in VCF format file
     * @return int[] {ref, ref_reverse, alt, alt_reverse}
     */
    private int[] parseDp4(String dp4String) {
        String[] info = dp4String.split(";"), readsCount = new String[4];
        for (String s: info) {
            if (s.startsWith("DP4")) {
                readsCount = s.split("=")[1].split(",");
            }
        }
        int[] reads = new int[4];
        for (int i=0; i<readsCount.length; i++) {
            reads[i] = Integer.parseInt(readsCount[i]);
        }

        return reads;
    }

    private boolean ifInExon(HashMap<String, IntervalTree> chrExon, String geneId, int position) {
        IntervalTree exons = chrExon.getOrDefault(geneId, null);
        if (exons == null)
            return false;

        IntervalTreeNode itn = exons.search(exons.root, position);

        return itn != null;
    }

    /**
     * get gene name
     * @param recordInfo GTF information
     * @return gene name
     */
    private String[] getGeneInfo(String recordInfo) {
        String[] info = recordInfo.split("; ");
        String geneName = null, geneId = null;
        for (String s: info) {
            if (s.startsWith("gene_id")) {
                String[] name = s.split(" ");
                geneId = name[1].substring(1, name[1].length() -1);
            }
            if (s.startsWith("gene_name")) {
                String[] name = s.split(" ");
                geneName = name[1].substring(1, name[1].length() -1);
            }
        }

        return new String[] {geneId, geneName};
    }
}
