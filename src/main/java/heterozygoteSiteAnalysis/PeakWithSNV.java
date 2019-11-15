package heterozygoteSiteAnalysis;

import GTFComponent.GeneExonIntervalTree;
import org.apache.commons.cli.*;

import java.io.*;
import java.math.BigDecimal;
import java.util.HashMap;

public class PeakWithSNV {
    private String gtfFile, bedFile, vcfFile, outputFile;
    private boolean exonMutation;
    // chr -> geneId -> peakTree;  chr -> geneId -> exonIntervalTree
    private HashMap<String, HashMap<String, IntervalTree>> peakTreeMap, geneExonIntervalTree;
    private HashMap<String, String> geneNames, peakGeneId;

    public static void main(String[] args) {
        Options options = new Options();
        CommandLine commandLine = setCommandLine(options, args);
        String bedFile, vcfFile, gtfFile, outputFile;
        boolean inExon = false;

        bedFile = new File(commandLine.getOptionValue("b")).getAbsolutePath();
        vcfFile = new File(commandLine.getOptionValue("v")).getAbsolutePath();
        gtfFile = new File(commandLine.getOptionValue("g")).getAbsolutePath();
        outputFile = new File(commandLine.getOptionValue("o")).getAbsolutePath();

        if (commandLine.hasOption("ex")) {
            int val = Integer.valueOf(commandLine.getOptionValue("ex"));
            if (val == 1)
                inExon = true;
            else if (val != 0){
                System.out.println("invalid value of parameter ex, must be 0 or 1");
                System.exit(2);
            }
        }

        PeakWithSNV pws = new PeakWithSNV(gtfFile, bedFile, vcfFile, outputFile, inExon);
        pws.buildPeakIntervalTree();
        pws.parseGTFFile();
        pws.parseVCFFile();
    }

    public PeakWithSNV(String gtfFile, String bedFile, String vcfFile, String outputFile, boolean exonMutation) {
        this.gtfFile = gtfFile;
        this.bedFile = bedFile;
        this.vcfFile = vcfFile;
        this.outputFile = outputFile;
        this.exonMutation = exonMutation;
    }

    private void buildPeakIntervalTree() {
        BufferedReader bfr = null;
        String chrNum, label, geneId, strand, blockSize, blockStart;
        int start, end, blockCount;
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
                    chrNum = info[0];
                    start = Integer.parseInt(new BigDecimal(info[1]).toPlainString());
                    end = Integer.parseInt(new BigDecimal(info[2]).toPlainString());
                    geneId = info[3];
                    strand = info[5];
                    blockCount = Integer.valueOf(info[9]);
                    blockSize = info[10];
                    blockStart = info[11];

                    HashMap<String, IntervalTree> chrTree = peakTreeMap.getOrDefault(chrNum, new HashMap<>());
                    if (blockCount > 1) {
                        IntervalTreeNode[] nodes = this.multiBlock(start, end, blockSize, blockStart);
                        for (IntervalTreeNode node : nodes) {
                            IntervalTree it = chrTree.getOrDefault(strand, new IntervalTree());
                            it = it.insertNode(it, node);
                            chrTree.put(strand, it);
                        }
                    } else {
                        IntervalTreeNode newNode = new IntervalTreeNode(start, end, start, end);
                        IntervalTree it = chrTree.getOrDefault(strand, new IntervalTree());
                        it = it.insertNode(it, newNode);
                        chrTree.put(strand, it);
                    }
                    peakTreeMap.put(chrNum, chrTree);
                    label = String.join("->",
                                        new String[] {new BigDecimal(info[1]).toPlainString(),
                                                      new BigDecimal(info[2]).toPlainString()});
                    this.peakGeneId.put(label, geneId);
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
        if (this.exonMutation)
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
                                                                 "geneId", "geneName", "refNc", "altNc", "refCount", "altCount"}));
            bfw.newLine();
            String readIn = "", writeOut, chrNum, refNc, altNc, geneId, geneName, peakStart, peakEnd, label;
            int position, referenceCount, alternativeCount;
            String[] info;
            int[] readsCount;
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
                    readsCount = this.parseDp4(info[7]);
                    referenceCount = readsCount[0] + readsCount[1];
                    alternativeCount = readsCount[2] + readsCount[3];

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
                            geneId = this.peakGeneId.get(label);
                            geneName = this.geneNames.getOrDefault(geneId, null);
                            if (geneName == null)
                                continue;
                            if (this.exonMutation && !this.ifInExon(chrNum, geneId, position))
                                continue;
                            writeOut = String.join("\t",
                                                    new String[] {chrNum, peakStart, peakEnd, info[1], geneId, geneName,
                                                    refNc, altNc, Integer.toString(referenceCount),
                                                    Integer.toString(alternativeCount)});
                            bfw.write(writeOut);
                            bfw.newLine();
                        }
                    }
                    if (negStrandTree != null) {
                        IntervalTreeNode negStrandSearchResult = negStrandTree.search(negStrandTree.root, position);
                        if (negStrandSearchResult != null) {
                            peakStart = Integer.toString(negStrandSearchResult.peakStart);
                            peakEnd = Integer.toString(negStrandSearchResult.peakEnd);
                            label = String.join("->", new String[] {peakStart, peakEnd});
                            geneId = this.peakGeneId.get(label);
                            geneName = this.geneNames.getOrDefault(geneId, null);
                            if (geneName == null)
                                continue;
                            if (this.exonMutation && !this.ifInExon(chrNum, geneId, position))
                                continue;
                            writeOut = String.join("\t",
                                    new String[] {chrNum, peakStart, peakEnd, info[1], geneId, geneName,
                                                refNc, altNc, Integer.toString(referenceCount),
                                                Integer.toString(alternativeCount)});
                            bfw.write(writeOut);
                            bfw.newLine();
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

        option = new Option("ex", "exon", true, "If the mutation must be in exon. 1 or 0, represent must or not, respectively. Default 0");
        option.setRequired(false);
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
