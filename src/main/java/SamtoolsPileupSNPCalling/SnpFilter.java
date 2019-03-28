package SamtoolsPileupSNPCalling;


import org.apache.log4j.Logger;

import java.io.*;
import java.util.HashMap;

public class SnpFilter {

    private File[] referenceVcfFiles;
    private File outputFile, rawVcfFile;
    private HashMap<String, Integer> cols;

    /**
     * use 1000 Genome project SNP to filter samtools pileup function output SNP
     * @param refVcfDirPath directory which stores 1000 Genome VCF file for each chromosome
     * @param rawVcfFile samtools pileup function output reads abundant file
     */
    public SnpFilter(String refVcfDirPath, String rawVcfFile, Logger log) {
        File referenceVcfDir = new File(refVcfDirPath);
        if (!referenceVcfDir.isDirectory()) {
            log.error("need reference VCF directory");
            System.exit(4);
        }
        this.referenceVcfFiles = referenceVcfDir.listFiles();
        if (this.referenceVcfFiles == null) {
            log.error("empty reference VCF directory");
            System.exit(5);
        }

        this.rawVcfFile = new File(rawVcfFile);
        String outputFileName = rawVcfFile.substring(0, rawVcfFile.lastIndexOf("_")) + "ReadsCount.txt";
        this.outputFile = new File(outputFileName);

        this.cols = new HashMap<>();
        this.cols.put("A", 4);
        this.cols.put("C", 5);
        this.cols.put("T", 6);
        this.cols.put("G", 7);
    }

    /**
     * filter SNP
     */
    public void filterVcf() {
        HashMap<String, BinarySearchTree> trees = this.buildTreeForChr();
        try {
            BufferedReader bfr = new BufferedReader(
                    new InputStreamReader(new FileInputStream(this.rawVcfFile))
            );

            BufferedWriter bwr = new BufferedWriter(
                    new OutputStreamWriter(new FileOutputStream(this.outputFile))
            );

            String readIn = "";
            String writeOut = "";
            while (readIn != null) {
                readIn = bfr.readLine();
                if (readIn != null) {
                    if (readIn.startsWith("#"))
                        continue;
                    String[] lineInfo = readIn.split("\t");
                    // skip letters "chr"
                    String chrNum = lineInfo[0];
                    int position = Integer.parseInt(lineInfo[1]);

                    BinarySearchTree targetTree = trees.get(chrNum);
                    VcfTreeNode searchRes = targetTree.search(targetTree.root, position);
                    if (searchRes != null) {
                        String id = searchRes.id;
                        String altNc = searchRes.altNucleotide;
                        if (!this.cols.containsKey(altNc))
                            continue;
                        // output information contains chr, position, id(1000Genome), reference nc, alternative nc, ref count, alt count
                        String[] usefulInfo = new String[]{lineInfo[0], lineInfo[1], id, lineInfo[2], altNc, lineInfo[this.cols.get(lineInfo[2])], lineInfo[this.cols.get(altNc)]};
                        StringBuffer stringBuf = new StringBuffer();
                        for (String str: usefulInfo) {
                            stringBuf.append(str);
                            stringBuf.append("\t");
                        }
                        writeOut = stringBuf.toString();
                        bwr.write(writeOut);
                        bwr.newLine();
                    }

                }
            }

            bwr.flush();
            bfr.close();
            bwr.close();
        } catch (FileNotFoundException fne) {
            fne.printStackTrace();
            System.exit(2);
        } catch (IOException ie) {
            ie.printStackTrace();
            System.exit(3);
        }
    }

    /**
     * build binary search tree for each chromosome SNP
     * @param referenceFile 1000 genome project SNP for a particular chromosome
     * @return binary search tree
     */
    private BinarySearchTree buildVcfTree(File referenceFile) {
        VcfSearchTree vst = new VcfSearchTree(referenceFile);

        return vst.buildTree();
    }

    /**
     * build binary search tree for all reference SNP
     * @return hashmap which contains binary search tree for each chromosome
     */
    private HashMap<String, BinarySearchTree> buildTreeForChr() {

        HashMap<String, BinarySearchTree> chrBinaryTree = new HashMap<String, BinarySearchTree>();

        for (File vcfFile : this.referenceVcfFiles) {
            if (!vcfFile.getAbsolutePath().endsWith(".vcf"))
                continue;
            String chrNum = this.getChrNum(vcfFile);
            if (chrNum == null) {
                System.out.println("Invalid VCF reference file " + vcfFile.getAbsolutePath());
                System.exit(6);
            }

            chrBinaryTree.put(chrNum, this.buildVcfTree(vcfFile));
        }

        return chrBinaryTree;
    }

    /**
     * get chromosome number of a reference VCF File
     * @param referenceVcf 1000 genome VCF file
     * @return chrNum
     */
    private String getChrNum(File referenceVcf) {
        String chrNum = null;
        try {
            BufferedReader bfr = new BufferedReader(
                    new InputStreamReader(new FileInputStream(referenceVcf))
            );

            String line;
            while ((line = bfr.readLine()) != null) {
                if (line.startsWith("#"))
                    continue;
                chrNum = line.split("\t")[0];
                break;
            }

            bfr.close();
        } catch (FileNotFoundException fne) {
            fne.printStackTrace();
            System.exit(2);
        } catch (IOException ie) {
            ie.printStackTrace();
            System.exit(3);
        }

        return chrNum;
    }
}
