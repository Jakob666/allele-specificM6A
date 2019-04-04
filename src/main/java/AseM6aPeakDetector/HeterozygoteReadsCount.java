package AseM6aPeakDetector;

import org.apache.log4j.Logger;

import java.io.*;
import java.util.HashMap;
import java.util.LinkedList;

public class HeterozygoteReadsCount {
    private File peakCoveredSNPFile;
    private Logger log;
    public HeterozygoteReadsCount(String peakCoveredSnpFile, Logger log) {
        this.peakCoveredSNPFile = new File(peakCoveredSnpFile);
        this.log = log;
    }

    /**
     * get major and minor haplotype SNP reads count for each m6a peak
     * @return HashMap
     */
    public HashMap<String, HashMap<String, HashMap<String, LinkedList<Integer>>>> getMajorMinorHaplotype() {
        HashMap<String, HashMap<String, HashMap<String, LinkedList<Integer>>>> majorMinorHaplotype = new HashMap<>();
        BufferedReader bfr;
        try {
            bfr = new BufferedReader(
                    new InputStreamReader(new FileInputStream(this.peakCoveredSNPFile))
            );
            String line = "";
            String[] info;
            String chr, peakStart, peakEnd, peakRange;
            int refCount, altCount;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    info = line.split("\t");
                    chr = info[0];
                    peakStart = info[3];
                    peakEnd = info[4];
                    peakRange = peakStart +":"+peakEnd;
                    refCount = Integer.parseInt(info[7]);
                    altCount = Integer.parseInt(info[8]);

                    HashMap<String, HashMap<String, LinkedList<Integer>>> chrMap = majorMinorHaplotype.getOrDefault(chr, null);
                    if (chrMap == null)
                        chrMap = new HashMap<String, HashMap<String, LinkedList<Integer>>>();
                    HashMap<String, LinkedList<Integer>> peakSnp = chrMap.getOrDefault(peakRange, null);
                    if (peakSnp == null)
                        peakSnp = new HashMap<String, LinkedList<Integer>>();
                    LinkedList<Integer> majorHaplotypeReads = peakSnp.getOrDefault("major", null);
                    if (majorHaplotypeReads == null)
                        majorHaplotypeReads = new LinkedList<Integer>();
                    LinkedList<Integer> minorHaplotypeReads = peakSnp.getOrDefault("minor", null);
                    if (minorHaplotypeReads == null)
                        minorHaplotypeReads = new LinkedList<Integer>();

                    majorHaplotypeReads.add(Math.max(refCount, altCount));
                    minorHaplotypeReads.add(Math.min(refCount, altCount));
                    peakSnp.put("major", majorHaplotypeReads);
                    peakSnp.put("minor", minorHaplotypeReads);
                    chrMap.put(peakRange, peakSnp);
                    majorMinorHaplotype.put(chr, chrMap);
                }
            }
            bfr.close();
        } catch (IOException ie) {
            this.log.error("load file failed");
            this.log.error(ie.getMessage());
        }

        return majorMinorHaplotype;
    }

}
