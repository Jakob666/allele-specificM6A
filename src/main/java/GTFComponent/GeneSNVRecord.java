package GTFComponent;

import heterozygoteSiteAnalysis.IntervalTree;
import heterozygoteSiteAnalysis.IntervalTreeNode;
import org.apache.commons.cli.*;
import org.apache.commons.math3.exception.OutOfRangeException;

import java.io.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;

public class GeneSNVRecord {
    private String gtfFile, vcfFile, outputFile;
    private HashMap<String, IntervalTree> geneIntervalTree;
    private HashMap<String, HashMap<String, IntervalTree>> geneExonIntervalTree;

    public static void main(String[] args) {
        Options options = new Options();
        CommandLine commandLine = setCommandLine(options, args);
        String gtfFile, vcfFile, outputFile;

        gtfFile = new File(commandLine.getOptionValue("g")).getAbsolutePath();
        vcfFile = new File(commandLine.getOptionValue("v")).getAbsolutePath();
        outputFile = new File(commandLine.getOptionValue("o")).getAbsolutePath();

        GeneSNVRecord gsr = new GeneSNVRecord(gtfFile, vcfFile, outputFile);
        gsr.parseGTFFile();
        gsr.parseVCFFile();
    }

    public GeneSNVRecord(String gtfFile, String vcfFile, String outputFile) {
        this.gtfFile = gtfFile;
        this.vcfFile = vcfFile;
        this.outputFile = outputFile;
    }

    public void locateSnv() {
        this.parseGTFFile();
        this.parseVCFFile();
    }

    private void parseGTFFile() {
        GTFIntervalTree git = new GTFIntervalTree(this.gtfFile);
        git.parseGTFFile();
        this.geneIntervalTree = git.getGtfIntervalTrees();
        GeneExonIntervalTree geit = new GeneExonIntervalTree(this.gtfFile);
        geit.generateExonTree();
        this.geneExonIntervalTree = geit.getGeneExonIntervalTree();
    }

    private void parseVCFFile() {
        BufferedReader bfr = null;
        BufferedWriter bfw = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.vcfFile))));
            bfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(this.outputFile))));
            bfw.write(String.join("\t", new String[] {"#chr", "geneId", "geneName", "position",
                                                                "majorAllele", "minorAllele",
                                                                "majorAlleleCount", "minorAlleleCount", "qualityLabel"}));
            bfw.newLine();
            String readIn = "", writeOut, chrNum, geneId, geneName, refNc, altNc, majorNc, minorNc, qualityLabel;
            int position, allele1Count, allele2Count, majorAlleleCount, minorAlleleCount;
            boolean specialMutation;
            String[] info;
            int[] readsCount;
            LinkedList<IntervalTreeNode> potentialLocateGene;
            while (readIn != null) {
                readIn = bfr.readLine();
                if (readIn != null) {
                    if (readIn.startsWith("#"))
                        continue;
                    info = readIn.split("\t");
                    chrNum = info[0];
                    IntervalTree it = this.geneIntervalTree.getOrDefault(chrNum, null);
                    if (it == null)
                        continue;
                    position = Integer.valueOf(info[1]);
                    GTFIntervalTreeNode gitn = (GTFIntervalTreeNode) it.search(it.root, position);
                    if (gitn == null)
                        continue;
                    potentialLocateGene = new LinkedList<>();
                    this.searchLocateGene(gitn, chrNum, position, potentialLocateGene);
                    if (potentialLocateGene.size() == 0) {
                        potentialLocateGene = null;
                        continue;
                    }

                    refNc = info[3];
                    altNc = info[4];
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

                    for (IntervalTreeNode itn: potentialLocateGene) {
                        GTFIntervalTreeNode node = (GTFIntervalTreeNode) itn;
                        geneId = node.geneId;
                        geneName = node.geneName;
                        writeOut = String.join("\t", new String[]{chrNum, geneId, geneName, info[1], majorNc, minorNc,
                                String.valueOf(majorAlleleCount), String.valueOf(minorAlleleCount), qualityLabel});
                        bfw.write(writeOut);
                        bfw.newLine();
                    }
                }
                writeOut = null;
                info = null;
                potentialLocateGene = null;
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

    private void searchLocateGene(GTFIntervalTreeNode gitn, String chrNum, int position,
                                  LinkedList<IntervalTreeNode> potentialLocateGene) {
        String geneId = gitn.geneId;
        boolean inExon = this.ifInExon(chrNum, geneId, position);
        if (inExon)
            potentialLocateGene.add(gitn);


        if (gitn.rightChild != null) {
            GTFIntervalTreeNode rc = (GTFIntervalTreeNode) gitn.rightChild;
            if (rc.intervalStart <= position && position <= rc.intervalEnd)
                this.searchLocateGene(rc, chrNum, position, potentialLocateGene);
        }
        if (gitn.leftChild != null) {
            GTFIntervalTreeNode lc = (GTFIntervalTreeNode) gitn.leftChild;
            if (lc.intervalStart <= position && position <= lc.intervalEnd)
                this.searchLocateGene(lc, chrNum, position, potentialLocateGene);
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

    public static CommandLine setCommandLine(Options options, String[] args) {
        Option option = new Option("g", "gtf", true, "GTF format file");
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
