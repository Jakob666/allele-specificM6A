package AseM6aPeakDetector;

import org.apache.log4j.Logger;

import java.io.*;
import java.util.HashMap;
import java.util.LinkedList;

public class HeterozygoteReadsCount {
    private File peakCoveredSNPFile;
    private HashMap<String, String> peakMajorAlleleStrand = new HashMap<>();
    private Logger log;

    /**
     * Constructor
     * @param peakCoveredSnpFile heterozygoteSiteAnalysis.PeakCoveredSNP类生成的文件
     * @param log 日志对象
     */
    public HeterozygoteReadsCount(String peakCoveredSnpFile, Logger log) {
        this.peakCoveredSNPFile = new File(peakCoveredSnpFile);
        this.log = log;
    }

    /**
     * 得到每个m6A信号major和minor haplotype上SNP reads count
     * @return [chr1: [peak1: [major: [count1, count2,...], minor: [count1, count2,...]]], chr2:....]
     */
    public HashMap<String, HashMap<String, HashMap<String, LinkedList<Integer>>>> getMajorMinorHaplotype() {
        HashMap<String, HashMap<String, HashMap<String, LinkedList<Integer>>>> majorMinorHaplotype = new HashMap<>();
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(
                    new InputStreamReader(new FileInputStream(this.peakCoveredSNPFile))
            );
            String line = "";
            String[] info;
            String chr, peakStart, peakEnd, peakRange, majorAlleleStrand;
            int refCount, altCount;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    chr = info[0];
                    peakStart = info[3];
                    peakEnd = info[4];
                    peakRange = peakStart +":"+peakEnd;
                    majorAlleleStrand = info[5];
                    refCount = Integer.parseInt(info[6]);
                    altCount = Integer.parseInt(info[7]);

                    this.peakMajorAlleleStrand.put(chr+":"+peakRange, majorAlleleStrand);
                    HashMap<String, HashMap<String, LinkedList<Integer>>> chrMap = majorMinorHaplotype.getOrDefault(chr, new HashMap<>());
                    HashMap<String, LinkedList<Integer>> peakSnp = chrMap.getOrDefault(peakRange, new HashMap<>());
                    LinkedList<Integer> majorHaplotypeReads = peakSnp.getOrDefault("major", new LinkedList<Integer>());
                    LinkedList<Integer> minorHaplotypeReads = peakSnp.getOrDefault("minor", new LinkedList<Integer>());

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
        } finally {
            if (bfr != null) {
                try {
                    bfr.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }

        return majorMinorHaplotype;
    }

    public HashMap<String, String> getPeakMajorAlleleStrand() {
        return peakMajorAlleleStrand;
    }

}
