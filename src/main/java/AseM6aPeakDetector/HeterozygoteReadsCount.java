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
     * @return [chr1: [peak1: position1: [major: count, minor: count], position2:[major: count, minor:count]], chr2:....]
     */
    public HashMap<String, HashMap<String, HashMap<String, HashMap<String, Integer>>>> getMajorMinorHaplotype() {
        HashMap<String, HashMap<String, HashMap<String, HashMap<String, Integer>>>> majorMinorHaplotype = new HashMap<>();
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(
                    new InputStreamReader(new FileInputStream(this.peakCoveredSNPFile))
            );
            String line = "";
            String[] info;
            String chr, position, peakStart, peakEnd, peakRange, majorAlleleStrand, majorNc, minorNc;
            int majorCount, minorCount;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    chr = info[0];
                    position = info[2];
                    peakStart = info[3];
                    peakEnd = info[4];
                    peakRange = peakStart +":"+peakEnd;
                    majorAlleleStrand = info[5];
                    majorNc = info[6];
                    minorNc = info[7];
                    majorCount = Integer.parseInt(info[6]);
                    minorCount = Integer.parseInt(info[7]);

                    this.peakMajorAlleleStrand.put(chr+":"+peakRange, majorAlleleStrand);
                    HashMap<String, HashMap<String, HashMap<String, Integer>>> chrMap = majorMinorHaplotype.getOrDefault(chr, new HashMap<>());
                    HashMap<String, HashMap<String, Integer>> peakSnp = chrMap.getOrDefault(peakRange, new HashMap<>());

                    HashMap<String, Integer> nucleotideReads = new HashMap<>();
                    nucleotideReads.put(majorNc, majorCount);
                    nucleotideReads.put(minorNc, minorCount);

                    peakSnp.put(position, nucleotideReads);
                    chrMap.put(peakRange, peakSnp);
                    majorMinorHaplotype.put(chr, chrMap);
                }
            }
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
