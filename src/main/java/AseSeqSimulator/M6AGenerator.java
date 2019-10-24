package AseSeqSimulator;

import GTFComponent.ElementRecord;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.UniformIntegerDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;

import java.io.*;
import java.util.*;

/**
 * randomly generate m6A sites for each mRNA
 */
public class M6AGenerator {
    private NormalDistribution m6aFrequency = new NormalDistribution(3, 1);
    private HashMap<String, HashMap<Integer, Double>> asmRatio = new HashMap<>();
    private HashMap<String, HashMap<Integer, Boolean>> asmBias = new HashMap<>();

    /**
     * Constructor
     */
    public M6AGenerator() {}

    /**
     * randomly generate m6A sites
     * @param m6aModifyGene Gene instance
     * @return m6A modification sites
     */
    public int[] generateM6aSites(Gene m6aModifyGene, int peakStart, int peakEnd, int readLength) {
        UniformIntegerDistribution uid;
        int modifyPosition;

        uid = new UniformIntegerDistribution(peakStart + readLength/2, peakEnd - readLength/2);

        // randomly select m6A modification sites on exon
        modifyPosition = uid.sample();

        int genomePosition;
        genomePosition = this.m6aGenomePosition(m6aModifyGene, modifyPosition);
        // release
        uid = null;

        return new int[] {modifyPosition, genomePosition};
    }

    /**
     * write the generated m6A sites into file
     * @param simulatedM6aSites HashMap<String, HashMap<Integer, Integer>> geneId: exon position: genome position
     * @param outputFile output file path
     */
    public void storeGeneM6aSites(HashMap<String, HashMap<String, HashMap<Integer, Integer>>> simulatedM6aSites, File outputFile,
                                  HashMap<String, ArrayList<String>> peakRanges, HashMap<String, ArrayList<Boolean>> asmPeaks) {
        UniformRealDistribution urd = new UniformRealDistribution(0.6, 0.9);
        BufferedWriter bfw = null;
        try {
            bfw = new BufferedWriter(
                    new OutputStreamWriter(new FileOutputStream(outputFile))
            );
            bfw.write("geneId\texonPosition\tgenomePosition\tmajorASMRatio\tMajorAlleleSpecific\n");

            for (String label: simulatedM6aSites.keySet()) {
                String[] info = label.split(":");
                String geneId = info[1];
                ArrayList<String> genePeaks = peakRanges.get(geneId);
                ArrayList<Boolean> asms = asmPeaks.get(geneId);
                HashMap<String, HashMap<Integer, Integer>> peakCoveredM6aSites = simulatedM6aSites.get(label);
                assert genePeaks.size() == asms.size();
                // get m6A ASM ratio
                for (int i = 0; i < genePeaks.size(); i++) {
                    boolean asm = asms.get(i);
                    String peak = genePeaks.get(i);
                    double asmRatio = 0.5;
                    if (asm)
                        asmRatio = urd.sample();
                    HashMap<Integer, Integer> m6aSites = peakCoveredM6aSites.get(peak);
                    for (Integer exonPosition: m6aSites.keySet()) {
                        Integer genomePosition = m6aSites.get(exonPosition);
                        bfw.write(geneId + "\t" + exonPosition + "\t" + genomePosition + "\t" + asmRatio + "\t" + asm);
                        bfw.newLine();

                        HashMap<Integer, Double> geneAsmRatio = this.asmRatio.getOrDefault(geneId, new HashMap<>());
                        HashMap<Integer, Boolean> geneAsmBias = this.asmBias.getOrDefault(geneId, new HashMap<>());
                        geneAsmRatio.put(exonPosition, asmRatio);
                        geneAsmBias.put(exonPosition, asm);
                        this.asmRatio.put(geneId, geneAsmRatio);
                        this.asmBias.put(geneId, geneAsmBias);
                    }
                }
            }
        } catch (IOException ie) {
            ie.printStackTrace();
        } finally {
            if (bfw != null) {
                try {
                    bfw.close();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
    }

    public HashMap<String, HashMap<Integer, Double>> getAseRatio() {return asmRatio;}

    public HashMap<String, HashMap<Integer, Boolean>> getAseBias() {return asmBias;}

    /**
     * get m6A site genome location via exon location
     * @param gene Gene instance
     * @param exonM6aPos m6A exon modification position
     * @return genome position
     */
    public int m6aGenomePosition(Gene gene, int exonM6aPos) {
        ElementRecord exon = gene.getExonList();
        int genomePosition;
        int length = 0, distance = exonM6aPos;
        while (exon != null) {
            length = length + exon.getElementEnd() - exon.getElementStart() + 1;
            if (length >= exonM6aPos)
                break;
            exon = exon.getNextElement();
            distance = exonM6aPos - length;
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
