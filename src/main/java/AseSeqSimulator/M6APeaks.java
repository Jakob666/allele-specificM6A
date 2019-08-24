package AseSeqSimulator;

import GTFComponent.ElementRecord;
import org.apache.commons.math3.distribution.UniformIntegerDistribution;

import java.util.HashMap;

public class M6APeaks {
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

    public void geneM6aPeakRange(Gene gene, int m6aPeakLength, int readLength) {
        int exonSeqLength = gene.getExonSeq().length();
        String geneId = gene.getGeneId(), chrNum = gene.getChr(), strand = gene.getStrand();
        String label = String.join(":", new String[]{chrNum, geneId, strand});
        // gene上最多有多少个peak，并确定peak center
        int maxPeakNum = exonSeqLength / (2 * m6aPeakLength);
        int peakNum = new UniformIntegerDistribution(0, maxPeakNum).sample();
        if (peakNum == 0) {
            int start = 1;
            int end = (start + m6aPeakLength < exonSeqLength)? (start + m6aPeakLength): exonSeqLength;
            int genomePeakStart = this.getGenomePosition(gene, start);
            int genomePeakEnd = this.getGenomePosition(gene, end);
            String[] peakGenomeRange = (genomePeakStart < genomePeakEnd)?
                    new String[]{Integer.toString(genomePeakStart), Integer.toString(genomePeakEnd),
                                 Integer.toString(start), Integer.toString(end)}:
                    new String[]{Integer.toString(genomePeakEnd), Integer.toString(genomePeakStart),
                                 Integer.toString(start), Integer.toString(end)};
            HashMap<Integer, Integer> m6aSites = this.geneM6aSites(gene, start, end, readLength);
            HashMap<String, HashMap<Integer, Integer>> peakCoverSites = new HashMap<>();
            peakCoverSites.put(String.join(":", peakGenomeRange), m6aSites);
            this.geneM6aSites.put(label, peakCoverSites);
        } else {
            for (int i = 0; i < peakNum; i++) {
                int start = i * (2 * m6aPeakLength) + m6aPeakLength / 2;
                int end = (i+1) * (2 * m6aPeakLength) - m6aPeakLength / 2;
                int peakCenter = new UniformIntegerDistribution(start, end).sample();
                int exonPeakStart = peakCenter - m6aPeakLength / 2;
                int exonPeakEnd = peakCenter + m6aPeakLength / 2;

                int genomePeakStart = this.getGenomePosition(gene, exonPeakStart);
                int genomePeakEnd = this.getGenomePosition(gene, exonPeakEnd);
                String[] peakGenomeRange = (genomePeakStart < genomePeakEnd)?
                                            new String[]{Integer.toString(genomePeakStart), Integer.toString(genomePeakEnd),
                                                         Integer.toString(start), Integer.toString(end)}:
                                            new String[]{Integer.toString(genomePeakEnd), Integer.toString(genomePeakStart),
                                                         Integer.toString(start), Integer.toString(end)};
                // 获取每个peak下覆盖的m6A修饰位点
                HashMap<Integer, Integer> m6aSites = this.geneM6aSites(gene, exonPeakStart, exonPeakEnd, readLength);
                HashMap<String, HashMap<Integer, Integer>> peakCoveredSites = this.geneM6aSites.getOrDefault(label, new HashMap<>());
                peakCoveredSites.put(String.join(":", peakGenomeRange), m6aSites);
                this.geneM6aSites.put(label, peakCoveredSites);
            }
        }
    }

    public HashMap<Integer, Integer> geneM6aSites(Gene gene, int start, int end, int readLength) {
        return this.m6AGenerator.generateM6aSites(gene, start, end, readLength);
    }

    /**
     * 将外显子序列的修饰位点转换为基因组上对应的位置
     * @param gene Gene对象
     * @param m6APos 外显子序列上的修饰位点
     * @return 基因组对应位置
     */
    private int getGenomePosition(Gene gene, int m6APos) {
        ElementRecord exon = gene.getExonList();
        int genomePosition;
        int length = 0, distance = m6APos;
        while (exon != null) {
            length = length + exon.getElementEnd() - exon.getElementStart() + 1;
            if (length >= m6APos)
                break;
            exon = exon.getNextElement();
            distance = m6APos - length;
        }
        if (gene.getStrand().equals("+")) {
            int start = exon.getElementStart();
            genomePosition = start + distance;
        } else {
            int end = exon.getElementEnd();
            genomePosition = end - distance;
        }

        return genomePosition;
    }
}
