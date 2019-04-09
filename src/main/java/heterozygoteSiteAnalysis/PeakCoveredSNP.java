package heterozygoteSiteAnalysis;

import org.apache.log4j.Logger;

import java.io.*;
import java.util.HashMap;

public class PeakCoveredSNP {
    private File vcfFile, peakCallingRes, outputFile;
    private Logger logger;

    /**
     * Constructor
     * @param vcfRecordFile vcf record file
     * @param peakCallingRes m6a peak calling result bed file
     * @param logger Logger instance
     */
    public PeakCoveredSNP(String vcfRecordFile, String peakCallingRes, Logger logger) {
        this.logger = logger;
        String outputFileName = vcfRecordFile.substring(0, vcfRecordFile.lastIndexOf("_")) + "_peakCoveredSNP.txt";
        this.vcfFile = new File(vcfRecordFile);
        if (!vcfFile.exists()) {
            this.logger.error("vcf file not exists");
            System.exit(2);
        }
        this.peakCallingRes = new File(peakCallingRes);
        if (!this.peakCallingRes.exists()) {
            this.logger.error("peak calling result bed file not exists");
            System.exit(2);
        }
        this.outputFile = new File(outputFileName);
    }

    /**
     * filter peaks which contains SNP, and record SNP information. Output in file.
     */
    public void filterSNPAndPeak() {
        HashMap<String, HashMap<String, IntervalTree>> m6aTreeMap = this.getM6aPeakTree();
        BufferedReader bfr;
        BufferedWriter bfw;
        try {
            bfr = new BufferedReader(
                    new InputStreamReader(new FileInputStream(this.vcfFile))
            );
            bfw = new BufferedWriter(
                    new OutputStreamWriter(new FileOutputStream(this.outputFile))
            );

            String line = "", writeOut;
            String[] info;
            String chrNum, refNc, altNc, refCount, altCount;
            int position;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    info = line.split("\t");
                    chrNum = info[0];
                    refNc = info[2];
                    altNc = info[3];
                    position = Integer.parseInt(info[1]);
                    refCount = info[4];
                    altCount = info[5];

                    HashMap<String, IntervalTree> chrTree = m6aTreeMap.getOrDefault(chrNum, null);
                    if (chrTree == null)
                        continue;
                    IntervalTree posStrandTree = chrTree.get("+");
                    IntervalTree negStrandTree = chrTree.get("-");
                    IntervalTreeNode posStrandSearchResult = posStrandTree.search(posStrandTree.root, position);
                    IntervalTreeNode negStrandSearchResult = negStrandTree.search(negStrandTree.root, position);
                    if (posStrandSearchResult != null) {
                        String[] newLine = new String[]{chrNum, "+", Integer.toString(position), Integer.toString(posStrandSearchResult.peakStart),
                                                        Integer.toString(posStrandSearchResult.peakEnd), refNc, altNc, refCount, altCount};
                        writeOut = String.join("\t", newLine);
                        bfw.write(writeOut);
                        bfw.newLine();
                    }else if (negStrandSearchResult != null) {
                        String[] newLine = new String[]{chrNum, "-", Integer.toString(position), Integer.toString(negStrandSearchResult.peakStart),
                                                        Integer.toString(negStrandSearchResult.peakEnd), refNc, altNc, refCount, altCount};
                        writeOut = String.join("\t", newLine);
                        bfw.write(writeOut);
                        bfw.newLine();
                    }
                }
            }

            bfw.flush();
            bfw.close();
            bfr.close();
        } catch (IOException ie) {
            this.logger.error("load peak calling result information failed.");
            this.logger.error(ie.getMessage());
        }
    }

    /**
     * get M6a peak trees using peak calling result bed file
     * @return HashMap
     */
    private HashMap<String, HashMap<String, IntervalTree>> getM6aPeakTree() {
        PeakIntervalTree pit = new PeakIntervalTree(this.peakCallingRes.getAbsolutePath(), this.logger);
        HashMap<String, HashMap<String, IntervalTree>> treeMap = pit.getPeakTrees();
        pit = null;

        return treeMap;
    }
}
