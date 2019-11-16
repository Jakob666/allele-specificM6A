package heterozygoteSiteAnalysis;

import GTFComponent.GeneExonIntervalTree;
import org.apache.log4j.Logger;

import java.io.*;
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
    public PeakCoveredSNP(String gtfFile, String vcfRecordFile, String peakCallingRes, String outputFile, int readsInfimum, Logger logger) {
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
            String chrNum, refNc, altNc, majorCount, minorCount, majorAlleleStrand, majorNc, minorNc;
            int[] refAndAltCount;
            int position;
            bfw.write("#chr\tgeneId\tgeneName\tpeakStrand\tposition\tpeakStart\tpeakEnd\tmajorAlleleType\tmajorNc\tminorNc\tmajorCount\tminorCount\n");
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
                    // only take single nucleotide mutation into consideration
                    if (refNc.length() > 1 | altNc.length() > 1)
                        continue;
                    position = Integer.parseInt(info[1]);
                    refAndAltCount = this.getReadsCountViaDp4(info[7]);
                    // filter homo SNV sites
                    if (refAndAltCount[0] == 0 |refAndAltCount[1] == 0)
                        continue;
                    // filter low coverage SNP sites
                    if (Math.max(refAndAltCount[0], refAndAltCount[1]) < this.readsInfimum)
                        continue;
                    if (refAndAltCount[0] >= refAndAltCount[1]) {
                        majorCount = Integer.toString(refAndAltCount[0]);
                        minorCount = Integer.toString(refAndAltCount[1]);
                        majorAlleleStrand = "ref";
                        majorNc = refNc;
                        minorNc = altNc;
                    } else {
                        majorCount = Integer.toString(refAndAltCount[1]);
                        minorCount = Integer.toString(refAndAltCount[0]);
                        majorAlleleStrand = "alt";
                        majorNc = altNc;
                        minorNc = refNc;
                    }

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
                                                        majorAlleleStrand, majorNc, minorNc, majorCount, minorCount};
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
                                                        majorAlleleStrand, majorNc, minorNc, majorCount, minorCount};
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
     * extract reference and alternative reads using DP4 info
     * @param record VCF INFO column
     * @return [refReadsCount, altReadsCount]
     */
    private int[] getReadsCountViaDp4(String record) {
        // record DP=5;SGB=-0.379885;RPB=1;MQB=1;MQSB=1;BQB=1;MQ0F=0;ICB=1;HOB=0.5;AC=0;AN=0;DP4=1,3,0,1;MQ=60
        String[] data = record.split(";");
        String key = null, value = null;
        String[] kAndv;
        for (String kv: data) {
            kAndv = kv.split("=");
            key = kAndv[0];
            if (key.equals("DP4")) {
                value = kAndv[1];
                break;
            }
        }
        if (value == null) {
            this.logger.error("can not find DP4 field in VCF file");
            System.exit(2);
        }
        String[] reads = value.split(",");
        int refReadsCount = Integer.parseInt(reads[0]) + Integer.parseInt(reads[1]);
        int altReadsCount = Integer.parseInt(reads[2]) + Integer.parseInt(reads[3]);

        return new int[]{refReadsCount, altReadsCount};
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
