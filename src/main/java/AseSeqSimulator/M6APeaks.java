package AseSeqSimulator;

import org.apache.commons.math3.distribution.UniformIntegerDistribution;

import java.util.HashMap;

public class M6APeaks {
    private int m6aPeakInterval = 500;
    private M6AGenerator m6AGenerator;
    private HashMap<String, HashMap<String, HashMap<Integer, Integer>>> geneM6aSites = new HashMap<>();

    private M6APeaks(M6AGenerator m6AGenerator) {
        this.m6AGenerator = m6AGenerator;
    }

    public static M6APeaks getInstance(M6AGenerator m6AGenerator) {
        return new M6APeaks(m6AGenerator);
    }

    public HashMap<String, HashMap<String, HashMap<Integer, Integer>>> getGeneM6aSites() {
        return geneM6aSites;
    }

    public void geneM6aPeakRange(Gene gene, int m6aPeakLength, int readLength, int peakNumSupremum) {
        int exonSeqLength = gene.getExonSeq().length();
        String geneId = gene.getGeneId(), chrNum = gene.getChr(), strand = gene.getStrand();
        String label = String.join(":", new String[]{chrNum, geneId, strand});
        // gene上最多有多少个peak，并确定peak center
        int maxPeakNum = exonSeqLength / (2 * m6aPeakLength + this.m6aPeakInterval);
        int peakNum;
        if (maxPeakNum != 0 && peakNumSupremum <= maxPeakNum) {
            if (peakNumSupremum <0)
                peakNum = new UniformIntegerDistribution(1, maxPeakNum).sample();
            else
                peakNum = new UniformIntegerDistribution(1, peakNumSupremum).sample();
        } else
            peakNum = 0;

        if (peakNum == 0) {
            int start = 1;
            int end = (start + m6aPeakLength < exonSeqLength)? (start + m6aPeakLength): exonSeqLength;
            int[] m6aSite = this.geneM6aSites(gene, start, end, readLength);
            int peakCenter = m6aSite[0];
            int peakStart = (peakCenter - m6aPeakLength/2 > 0)? peakCenter-m6aPeakLength/2:1;
            int peakEnd = (peakCenter + m6aPeakLength/2 < exonSeqLength)? peakCenter+m6aPeakLength/2:exonSeqLength;
            int genomePeakStart = this.m6AGenerator.m6aGenomePosition(gene, peakStart);
            int genomePeakEnd = this.m6AGenerator.m6aGenomePosition(gene, peakEnd);
            String[] peakGenomeRange = (genomePeakStart < genomePeakEnd)?
                    new String[]{Integer.toString(genomePeakStart), Integer.toString(genomePeakEnd),
                                 Integer.toString(peakStart), Integer.toString(peakEnd)}:
                    new String[]{Integer.toString(genomePeakEnd), Integer.toString(genomePeakStart),
                                 Integer.toString(peakStart), Integer.toString(peakEnd)};
            HashMap<Integer, Integer> modifySite = new HashMap<>();
            modifySite.put(m6aSite[0], m6aSite[1]);
            HashMap<String, HashMap<Integer, Integer>> peakCoverSites = new HashMap<>();
            peakCoverSites.put(String.join(":", peakGenomeRange), modifySite);
            this.geneM6aSites.put(label, peakCoverSites);
        } else {
            for (int i = 0; i < peakNum; i++) {
                double randNum = Math.random();
                if (randNum < 0.5 && this.geneM6aSites.keySet().contains(label))
                    continue;
                int start = i * (2 * m6aPeakLength + this.m6aPeakInterval) + m6aPeakLength / 2;
                int end = start + m6aPeakLength;
                // 获取每个peak下覆盖的m6A修饰位点
                int[] m6aSite = this.geneM6aSites(gene, start, end, readLength);
                int peakStart = m6aSite[0] - m6aPeakLength/2;
                int peakEnd = m6aSite[0] + m6aPeakLength/2;
                int genomePeakStart = this.m6AGenerator.m6aGenomePosition(gene, peakStart);
                int genomePeakEnd = this.m6AGenerator.m6aGenomePosition(gene, peakEnd);

                String[] peakGenomeRange = (genomePeakStart < genomePeakEnd)?
                                            new String[]{Integer.toString(genomePeakStart), Integer.toString(genomePeakEnd),
                                                         Integer.toString(peakStart), Integer.toString(peakEnd)}:
                                            new String[]{Integer.toString(genomePeakEnd), Integer.toString(genomePeakStart),
                                                         Integer.toString(peakStart), Integer.toString(peakEnd)};
                HashMap<Integer, Integer> modifySite = new HashMap<>();
                modifySite.put(m6aSite[0], m6aSite[1]);
                HashMap<String, HashMap<Integer, Integer>> peakCoveredSites = this.geneM6aSites.getOrDefault(label, new HashMap<>());
                peakCoveredSites.put(String.join(":", peakGenomeRange), modifySite);
                this.geneM6aSites.put(label, peakCoveredSites);
            }
        }
    }

    public int[] geneM6aSites(Gene gene, int start, int end, int readLength) {
        return this.m6AGenerator.generateM6aSites(gene, start, end, readLength);
    }
}
