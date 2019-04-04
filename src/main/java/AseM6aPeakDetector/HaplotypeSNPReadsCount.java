package AseM6aPeakDetector;

import org.apache.log4j.Logger;

import java.io.*;
import java.util.HashMap;
import java.util.LinkedList;

public class HaplotypeSNPReadsCount {
    private File peakCoveredSNPFile;
    private Logger log;

    public HaplotypeSNPReadsCount(String peakCoveredSNPFile, Logger log) {
        this.peakCoveredSNPFile = new File(peakCoveredSNPFile);
        this.log = log;
    }

    /**
     * get major and minor haplotype SNP reads count of all SNP site
     * @return HashMap
     */
    public HashMap<String, LinkedList<Integer>> haplotypeSnpReadsCount() {
        HashMap<String, LinkedList<Integer>> haplotypeSNP = new HashMap<>();
        BufferedReader bfr;
        try {
            bfr = new BufferedReader(
                    new InputStreamReader(new FileInputStream(this.peakCoveredSNPFile))
            );
            String line = "";
            String[] info;
            int majorCount, minorCount;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    info = line.split("\t");
                    majorCount = Integer.parseInt(info[7]);
                    minorCount = Integer.parseInt(info[8]);

                    LinkedList<Integer> majorHaplotype = haplotypeSNP.getOrDefault("major", null);
                    LinkedList<Integer> minorHaplotype = haplotypeSNP.getOrDefault("minor", null);
                    if (majorHaplotype == null)
                        majorHaplotype = new LinkedList<>();
                    if (minorHaplotype == null)
                        minorHaplotype = new LinkedList<>();
                    majorHaplotype.add(majorCount);
                    minorHaplotype.add(minorCount);

                    haplotypeSNP.put("major", majorHaplotype);
                    haplotypeSNP.put("minor", minorHaplotype);
                }
            }
            bfr.close();
        } catch (IOException ie) {
            log.error("can not reads file");
            log.error(ie.getMessage());
        }

        return haplotypeSNP;
    }
}
