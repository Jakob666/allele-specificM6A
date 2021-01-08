package heterozygoteSiteAnalysis;

import GTFComponent.GTFIntervalTree;
import GTFComponent.GTFIntervalTreeNode;
import GTFComponent.GeneExonIntervalTree;
import org.apache.commons.cli.*;

import java.io.*;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;

public class PeakWithSNV {
    private String gtfFile, bedFile, vcfFile, outputFile;
    private HashMap<String, IntervalTree> geneIntervalTree;
    // chr -> geneId -> peakTree;  chr -> geneId -> exonIntervalTree
    private HashMap<String, HashMap<String, IntervalTree>> peakTreeMap, geneExonIntervalTree;
    private HashMap<String, String> geneNames;
    private HashMap<String, ArrayList<String>> peakGeneId;

    public static void main(String[] args) {
        Options options = new Options();
        CommandLine commandLine = setCommandLine(options, args);
        String bedFile, vcfFile, gtfFile, outputFile;

        bedFile = new File(commandLine.getOptionValue("b")).getAbsolutePath();
        vcfFile = new File(commandLine.getOptionValue("v")).getAbsolutePath();
        gtfFile = new File(commandLine.getOptionValue("g")).getAbsolutePath();
        outputFile = new File(commandLine.getOptionValue("o")).getAbsolutePath();

        PeakWithSNV pws = new PeakWithSNV(gtfFile, bedFile, vcfFile, outputFile);
        pws.buildPeakIntervalTree();
        pws.parseGTFFile();
        pws.parseVCFFile();
    }

    public PeakWithSNV(String gtfFile, String bedFile, String vcfFile, String outputFile) {
        this.gtfFile = gtfFile;
        this.bedFile = bedFile;
        this.vcfFile = vcfFile;
        this.outputFile = outputFile;
    }

    public void locateSnvInPeak() {
        this.buildGTFIntervalTree();
        this.buildPeakIntervalTree();
        this.parseGTFFile();
        this.parseVCFFile();
    }

    private void buildGTFIntervalTree() {
        GTFIntervalTree git = new GTFIntervalTree(this.gtfFile);
        git.parseGTFFile();
        this.geneIntervalTree = git.getGtfIntervalTrees();
    }

    private void buildPeakIntervalTree() {
        BufferedReader bfr = null;
        String chrNum, label, geneId, strand, blockSize, blockStart;
        int start, end, peakCenter, blockCount;
        ArrayList<String> peakCoveredGenes;
        LinkedList<GTFIntervalTreeNode> potentialLocateGene;
        this.peakTreeMap = new HashMap<>();
        this.peakGeneId = new HashMap<>();
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.bedFile))));
            String line = "";
            String[] info;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    // BED format file must contains chr, peakStart, peakEnd
                    chrNum = info[0];
                    start = Integer.parseInt(new BigDecimal(info[1]).toPlainString());
                    end = Integer.parseInt(new BigDecimal(info[2]).toPlainString());
                    peakCenter = (int) Math.round((start + end) * 0.5);
                    IntervalTree it = this.geneIntervalTree.getOrDefault(chrNum, null);
                    if (it == null)
                        continue;
                    GTFIntervalTreeNode gitn = (GTFIntervalTreeNode) it.search(it.root, peakCenter);
                    if (gitn == null)
                        continue;
                    potentialLocateGene = new LinkedList<>();
                    this.searchLocateGene(gitn, peakCenter, potentialLocateGene);

                    HashMap<String, IntervalTree> chrTree = peakTreeMap.getOrDefault(chrNum, new HashMap<>());
                    peakCoveredGenes = new ArrayList<>();
                    for (GTFIntervalTreeNode node: potentialLocateGene) {
                        strand = node.strand;
                        IntervalTreeNode newNode = new IntervalTreeNode(start, end, start, end);
                        IntervalTree tree = chrTree.getOrDefault(strand, new IntervalTree());
                        tree = tree.insertNode(tree, newNode);
                        chrTree.put(strand, tree);
                        peakCoveredGenes.add(node.geneId);
                    }
                    peakTreeMap.put(chrNum, chrTree);
                    label = String.join("->",
                                        new String[] {new BigDecimal(info[1]).toPlainString(),
                                                      new BigDecimal(info[2]).toPlainString()});
                    this.peakGeneId.put(label, peakCoveredGenes);
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

    private IntervalTreeNode[] multiBlock(int start, int end, String blockSize, String blockStarts) {
        String[] sizes = blockSize.split(",");
        String[] starts = blockStarts.split(",");

        IntervalTreeNode[] nodes = new IntervalTreeNode[starts.length];
        for (int i = 0; i < starts.length; i++) {
            int intervalStart = start + Integer.parseInt(starts[i]);
            int intervalEnd = intervalStart + Integer.parseInt(sizes[i]);
            IntervalTreeNode newNode = new IntervalTreeNode(intervalStart, intervalEnd, start, end);
            nodes[i] = newNode;
        }

        return nodes;
    }

    private void parseGTFFile() {
        this.geneIdMatchName();
        this.generateExonTree();
    }

    private void geneIdMatchName() {
        this.geneNames = new HashMap<>();
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.gtfFile))));
            String line = "", geneId, geneName;
            String[] info, geneInfo;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    if (!info[2].equals("gene"))
                        continue;
                    geneInfo = this.getGeneInfo(info[8]);
                    geneId = geneInfo[0];
                    geneName = geneInfo[1];
                    this.geneNames.put(geneId, geneName);
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

    private void generateExonTree() {
        GeneExonIntervalTree geit = new GeneExonIntervalTree(this.gtfFile);
        geit.generateExonTree();
        this.geneExonIntervalTree = geit.getGeneExonIntervalTree();
    }

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

    private void parseVCFFile() {
        BufferedReader bfr = null;
        BufferedWriter bfw = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.vcfFile))));
            bfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(this.outputFile))));
            bfw.write(String.join("\t", new String[] {"#chr", "peakStart", "peakEnd", "mutatePosition",
                    "geneId", "geneName", "qualityScore", "qualityLabel", "majorAllele", "minorAllele",
                    "majorAlleleCount", "minorAlleleCount"}));
            bfw.newLine();
            String readIn = "", writeOut, chrNum, qualityScore, qualityLabel, refNc, altNc, majorNc, minorNc,
                   geneName, peakStart, peakEnd, label;
            int position, allele1Count, allele2Count, majorAlleleCount, minorAlleleCount;
            boolean specialMutation;
            String[] info;
            int[] readsCount;
            ArrayList<String> geneIds;
            IntervalTree posStrandTree, negStrandTree;
            while (readIn != null) {
                readIn = bfr.readLine();
                if (readIn != null) {
                    if (readIn.startsWith("#"))
                        continue;
                    info = readIn.split("\t");
                    chrNum = info[0];
                    position = Integer.valueOf(info[1]);
                    refNc = info[3];
                    altNc = info[4];
                    qualityScore = info[5];
                    qualityLabel = info[6];
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
                                System.out.println("invalid VCF record: " + readIn);
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
                            System.out.println("invalid VCF record: " + readIn);
                            System.exit(2);
                        }
                        allele1Count = readsCount[0] + readsCount[1];
                        allele2Count = readsCount[2] + readsCount[3];
                        majorAlleleCount = (allele1Count >= allele2Count)? allele1Count: allele2Count;
                        minorAlleleCount = (allele1Count < allele2Count)? allele1Count: allele2Count;
                        majorNc = (allele1Count >= allele2Count)? refNc : altNc.split(",")[0];
                        minorNc = (allele1Count < allele2Count)? refNc: altNc.split(",")[0];
                    }


                    HashMap<String, IntervalTree> trees = this.peakTreeMap.getOrDefault(chrNum, null);
                    if (trees == null)
                        continue;
                    posStrandTree = trees.getOrDefault("+", null);
                    negStrandTree = trees.getOrDefault("-", null);
                    if (posStrandTree != null) {
                        IntervalTreeNode posStrandSearchResult = posStrandTree.search(posStrandTree.root, position);
                        if (posStrandSearchResult != null) {
                            peakStart = Integer.toString(posStrandSearchResult.peakStart);
                            peakEnd = Integer.toString(posStrandSearchResult.peakEnd);
                            label = String.join("->", new String[] {peakStart, peakEnd});
                            geneIds = this.peakGeneId.get(label);
                            for (String geneId: geneIds) {
                                geneName = this.geneNames.getOrDefault(geneId, null);
                                if (geneName == null)
                                    continue;
                                if (!this.ifInExon(chrNum, geneId, position))
                                    continue;
                                writeOut = String.join("\t",
                                        new String[] {chrNum, peakStart, peakEnd, info[1], geneId, geneName, qualityScore,
                                                qualityLabel, majorNc, minorNc, String.valueOf(majorAlleleCount),
                                                String.valueOf(minorAlleleCount)});
                                bfw.write(writeOut);
                                bfw.newLine();
                            }
                        }
                    }
                    if (negStrandTree != null) {
                        IntervalTreeNode negStrandSearchResult = negStrandTree.search(negStrandTree.root, position);
                        if (negStrandSearchResult != null) {
                            peakStart = Integer.toString(negStrandSearchResult.peakStart);
                            peakEnd = Integer.toString(negStrandSearchResult.peakEnd);
                            label = String.join("->", new String[] {peakStart, peakEnd});
                            geneIds = this.peakGeneId.get(label);
                            for (String geneId: geneIds) {
                                geneName = this.geneNames.getOrDefault(geneId, null);
                                if (geneName == null)
                                    continue;
                                if (!this.ifInExon(chrNum, geneId, position))
                                    continue;
                                writeOut = String.join("\t",
                                        new String[] {chrNum, peakStart, peakEnd, info[1], geneId, geneName, qualityScore,
                                                qualityLabel, majorNc, minorNc, String.valueOf(majorAlleleCount),
                                                String.valueOf(minorAlleleCount)});
                                bfw.write(writeOut);
                                bfw.newLine();
                            }
                        }
                    }
                }
                info = null;
                readsCount = null;
                writeOut = null;
            }
            this.geneNames = null;
            this.peakGeneId = null;
            this.peakTreeMap = null;
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
            if (bfw != null) {
                try {
                    bfw.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
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

    private void searchLocateGene(GTFIntervalTreeNode gitn, int position,
                                  LinkedList<GTFIntervalTreeNode> potentialLocateGene) {
        potentialLocateGene.add(gitn);

        if (gitn.rightChild != null) {
            GTFIntervalTreeNode rc = (GTFIntervalTreeNode) gitn.rightChild;
            if (rc.intervalStart <= position && position <= rc.intervalEnd)
                this.searchLocateGene(rc, position, potentialLocateGene);
        }
        if (gitn.leftChild != null) {
            GTFIntervalTreeNode lc = (GTFIntervalTreeNode) gitn.leftChild;
            if (lc.intervalStart <= position && position <= lc.intervalEnd)
                this.searchLocateGene(lc, position, potentialLocateGene);
        }
    }

    private boolean ifInExon(String chrNum, String geneId, int position) {
        HashMap<String, IntervalTree> treeMap = this.geneExonIntervalTree.getOrDefault(chrNum, null);
        if (treeMap == null)
            return false;
        IntervalTree it = treeMap.getOrDefault(geneId, null);
        if (it == null)
            return false;
        IntervalTreeNode itn = it.search(it.root, position);

        return itn != null;
    }

    private static CommandLine setCommandLine(Options options, String[] args) {
        Option option = new Option("g", "gtf", true, "GTF format file");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("b", "bed", true, "BED format file");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("v", "vcf", true, "VCF format file");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("o", "output", true, "Output file path");
        option.setRequired(true);
        options.addOption(option);

        CommandLineParser parser = new DefaultParser();
        CommandLine commandLine = null;
        try {
            commandLine = parser.parse(options, args);
        } catch (ParseException pe) {
            pe.printStackTrace();
            System.exit(2);
        }

        return commandLine;
    }
}
