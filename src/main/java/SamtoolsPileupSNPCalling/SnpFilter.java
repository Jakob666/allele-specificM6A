package SamtoolsPileupSNPCalling;


import org.apache.log4j.Logger;

import java.io.*;
import java.util.*;

public class SnpFilter {

    private File referenceVcfFile;
    private File outputFile, rawVcfFile;
    private Logger log;

    /**
     * use 1000 Genome project SNP to filter samtools pileup function output SNP
     * @param refVcfFile directory which stores 1000 Genome VCF file for each chromosome
     * @param rawVcfFile samtools pileup function output reads abundant file
     */
    public SnpFilter(String refVcfFile, String rawVcfFile, Logger log) {
        this.log = log;
        this.referenceVcfFile = new File(refVcfFile);
        if (this.referenceVcfFile.isDirectory()) {
            this.log.error("need reference VCF file, not a directory");
            System.exit(2);
        }

        this.rawVcfFile = new File(rawVcfFile);
        String outputFileName = rawVcfFile.substring(0, rawVcfFile.lastIndexOf("_")) + "_readsCount.txt";
        this.outputFile = new File(outputFileName);
    }

    public static void main(String[] args) {
        System.setProperty("log_home", "/data1/hubs/test_output/MT4_T-cells_Control/INPUT");
        Logger log = Logger.getLogger(SnpFilter.class);
        SnpFilter snpFilter = new SnpFilter("/data/hbs/dbsnp/dbsnp.vcf",
                "/data1/hubs/test_output/MT4_T-cells_Control/INPUT/SRR2648293_abundant.txt", log);
        snpFilter.filterVcf();
    }

    /**
     * filter SNP
     */
    public void filterVcf() {
        VcfSearch vs = new VcfSearch(this.referenceVcfFile, this.log);
        HashMap<String, ArrayList<Integer>> vcfPositions = vs.vcfList();
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
                    String chrNum = lineInfo[0];
                    int position = Integer.parseInt(lineInfo[1]);

                    ArrayList<Integer> positions = vcfPositions.get(chrNum);
                    boolean searchRes = vs.binarySearch(positions, position);
                    if (searchRes) {
                        HashMap<Integer, String> cols = new HashMap<Integer, String>();
                        cols.put(Integer.parseInt(lineInfo[4]), "A");
                        cols.put(Integer.parseInt(lineInfo[5]), "C");
                        cols.put(Integer.parseInt(lineInfo[6]), "T");
                        cols.put(Integer.parseInt(lineInfo[7]), "G");
                        HashMap<String, Integer> majMinHaplotype = getMajorMinorHaplotype(cols);
                        String refNc = cols.get(majMinHaplotype.get("first"));
                        String altNc = cols.get(majMinHaplotype.get("second"));
                        int refCount = majMinHaplotype.get("first");
                        int altCount = majMinHaplotype.get("second");
                        // output chr, position, ref-nucleotide, alt-nucleotide, ref-reads count, alt-reads count
                        String[] outputInfo = new String[]{chrNum, Integer.toString(position), refNc, altNc,
                                                           Integer.toString(refCount), Integer.toString(altCount)};
                        writeOut = String.join("\t", outputInfo);
                        bwr.write(writeOut);
                        bwr.newLine();
                    }
                }
            }
            vcfPositions = null;
            bwr.flush();
            bfr.close();
            bwr.close();
        } catch (IOException ie) {
            this.log.error(ie.getMessage());
            System.exit(2);
        }
    }

//    /**
//     * build binary search tree for each chromosome SNP. Deprecated, use class VcfSearch instead
//     * @param referenceFile 1000 genome project SNP for a particular chromosome
//     * @return binary search tree for each chromosome
//     */
//    private HashMap<String, BinarySearchTree> buildVcfTree(File referenceFile) {
//        VcfSearchTree vst = new VcfSearchTree(referenceFile, this.log);
//
//        return vst.buildTree();
//    }

    /**
     * get base on major haplotype and minor haplotype
     * @param readsCounts HashMap, key is the reads cover SNP, value is base
     * @return HashMap
     */
    private HashMap<String, Integer> getMajorMinorHaplotype(HashMap<Integer, String> readsCounts) {
        //get key set
        ArrayList<Integer> list= new ArrayList<>(readsCounts.keySet());
        //override sort function, descend
        Collections.sort(list, new Comparator<Integer>() {
            @Override
            public int compare(Integer o1, Integer o2) {
                return o1 < o2 ? 1:-1;
            }
        });

        int first = list.get(0);
        int second = list.get(1);
        HashMap<String, Integer> map = new HashMap<>();
        map.put("first", first);
        map.put("second", second);

        return map;
    }

}
