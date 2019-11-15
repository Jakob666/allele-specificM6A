package GTFComponent;

import heterozygoteSiteAnalysis.IntervalTree;
import heterozygoteSiteAnalysis.IntervalTreeNode;
import org.apache.commons.cli.*;

import java.io.*;
import java.util.HashMap;
import java.util.LinkedList;

public class GeneSNVRecord {
    private String gtfFile, vcfFile, outputFile;
    private boolean exonMutation;
    private HashMap<String, IntervalTree> geneIntervalTree;
    private HashMap<String, HashMap<String, IntervalTree>> geneExonIntervalTree;

    public static void main(String[] args) {
        Options options = new Options();
        CommandLine commandLine = setCommandLine(options, args);
        String gtfFile, vcfFile, outputFile;
        boolean inExon = false;

        gtfFile = new File(commandLine.getOptionValue("g")).getAbsolutePath();
        vcfFile = new File(commandLine.getOptionValue("v")).getAbsolutePath();
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

        GeneSNVRecord gsr = new GeneSNVRecord(gtfFile, vcfFile, outputFile, inExon);
        gsr.parseGTFFile();
        gsr.parseVCFFile();
    }

    public GeneSNVRecord(String gtfFile, String vcfFile, String outputFile, boolean inExon) {
        this.gtfFile = gtfFile;
        this.vcfFile = vcfFile;
        this.outputFile = outputFile;
        this.exonMutation = inExon;
    }

    private void parseGTFFile() {
        GTFIntervalTree git = new GTFIntervalTree(this.gtfFile);
        git.parseGTFFile();
        this.geneIntervalTree = git.getGtfIntervalTrees();
        if (this.exonMutation) {
            GeneExonIntervalTree geit = new GeneExonIntervalTree(this.gtfFile);
            geit.generateExonTree();
            this.geneExonIntervalTree = geit.getGeneExonIntervalTree();
        } else
            this.geneExonIntervalTree = null;
    }

    private void parseVCFFile() {
        BufferedReader bfr = null;
        BufferedWriter bfw = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.vcfFile))));
            bfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(this.outputFile))));
            bfw.write(String.join("\t", new String[] {"#chr", "geneId", "geneName", "position",
                                                                "ref", "alt", "refCount", "altCount"}));
            bfw.newLine();
            String readIn = "", writeOut, chrNum, geneId, geneName, refNc, altNc;
            int position, referenceCount, alternativeCount;
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
                    readsCount = this.parseDp4(info[7]);
                    referenceCount = readsCount[0] + readsCount[1];
                    alternativeCount = readsCount[2] + readsCount[3];

                    for (IntervalTreeNode itn: potentialLocateGene) {
                        GTFIntervalTreeNode node = (GTFIntervalTreeNode) itn;
                        geneId = node.geneId;
                        geneName = node.geneName;
                        writeOut = String.join("\t", new String[]{chrNum, geneId, geneName, info[1], refNc, altNc,
                                               String.valueOf(referenceCount), String.valueOf(alternativeCount)});
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
        if (this.exonMutation) {
            boolean inExon = this.ifInExon(chrNum, geneId, position);
            if (inExon)
                potentialLocateGene.add(gitn);
        } else
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
