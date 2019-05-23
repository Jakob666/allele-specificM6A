package AseM6aPeakDetector;

import org.apache.log4j.Logger;

import java.io.*;
import java.util.HashMap;
import java.util.LinkedList;

public class HaplotypeSNPReadsCount {
    private File peakCoveredSNPFile;
    private Logger log;

    /**
     * Constructor
     * @param peakCoveredSNPFile heterozygoteSiteAnalysis.PeakCoveredSNP类生成的文件
     * @param log 日志对象
     */
    public HaplotypeSNPReadsCount(String peakCoveredSNPFile, Logger log) {
        this.peakCoveredSNPFile = new File(peakCoveredSNPFile);
        this.log = log;
    }

    /**
     * 得到major haplotype和minor haplotype上每个ASE位点的reads count
     * @return [major: [count1, count2..., countN], minor: [count1, count2..., countN]]
     */
    public HashMap<String, LinkedList<Integer>> haplotypeSnpReadsCount() {
        HashMap<String, LinkedList<Integer>> haplotypeSNP = new HashMap<>();
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(
                    new InputStreamReader(new FileInputStream(this.peakCoveredSNPFile))
            );
            String line = "";
            String[] info;
            int refCount, altCount, majorCount, minorCount;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    // 获取文件第8和第9列的信息，分别是ref和alt的read数目
                    refCount = Integer.parseInt(info[7]);
                    altCount = Integer.parseInt(info[8]);
                    if (refCount > altCount) {
                        majorCount = refCount;
                        minorCount = altCount;
                    } else {
                        majorCount = altCount;
                        minorCount = refCount;
                    }

                    LinkedList<Integer> majorHaplotype = haplotypeSNP.getOrDefault("major", new LinkedList<>());
                    LinkedList<Integer> minorHaplotype = haplotypeSNP.getOrDefault("minor", new LinkedList<>());
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
        } finally {
            if (bfr != null) {
                try {
                    bfr.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }

        return haplotypeSNP;
    }
}
