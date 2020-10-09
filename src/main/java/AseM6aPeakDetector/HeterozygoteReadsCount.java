package AseM6aPeakDetector;

import org.apache.log4j.Logger;

import java.io.*;
import java.util.HashMap;
import java.util.LinkedList;

public class HeterozygoteReadsCount {
    private int readsCountThreshold;
    private File peakCoveredSNPFile;
    private HashMap<String, String> peakCoveredGene = new HashMap<>();
    private HashMap<String, LinkedList<String>> peakMajorAlleleNc = new HashMap<>();
    private Logger log;

    /**
     * Constructor
     * @param peakCoveredSnpFile file which record peak covered SNV information
     * @param readsCountThreshold reads count threshold
     * @param log log4j instance
     */
    public HeterozygoteReadsCount(String peakCoveredSnpFile, int readsCountThreshold, Logger log) {
        this.readsCountThreshold = readsCountThreshold;
        this.peakCoveredSNPFile = new File(peakCoveredSnpFile);
        this.log = log;
    }

    /**
     * allele reads count of each SNV sites under m6A signal range
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
            String chr, position, peakStart, peakEnd, peakRange, geneId, geneName, majorNc, minorNc;
            int majorCount, minorCount;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    chr = info[0];
                    peakStart = info[1];
                    peakEnd = info[2];
                    position = info[3];
                    peakRange = peakStart +":"+peakEnd;
                    geneId = info[4];
                    majorNc = info[6];
                    minorNc = info[7];
                    majorCount = Integer.parseInt(info[8]);
                    minorCount = Integer.parseInt(info[9]);
                    if (majorCount < this.readsCountThreshold)
                        continue;

                    this.peakCoveredGene.put(chr+":"+peakStart+":"+peakEnd, geneId);
                    String label = chr+":"+peakRange;
                    LinkedList<String> majorAlleleNcRecords = this.peakMajorAlleleNc.getOrDefault(label, new LinkedList<>());
                    majorAlleleNcRecords.add(String.join(":", new String[]{position, majorNc}));
                    this.peakMajorAlleleNc.put(label, majorAlleleNcRecords);
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

    public HashMap<String, LinkedList<String>> getPeakMajorAlleleNucleotides() {
        return peakMajorAlleleNc;
    }

    public HashMap<String, String> getPeakCoveredGene() {
        return this.peakCoveredGene;
    }

}
