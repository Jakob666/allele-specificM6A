package AseM6aPeakDetector;

import betaBinomialMetaAnalysis.RhoEstimator;
import betaBinomialMetaAnalysis.SignificantTest;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.*;

public class AseM6aPeakDetector {
    private String peakCoveredSNPFile;
    private File aseM6aPeakFile;
    private double initialRho, learningRate, unImproveThreshold;
    private Logger log;

    public AseM6aPeakDetector(String peakCoveredSNPFile, double initialRho, double learningRate,
                              double unImproveThreshold, Logger log) {
        this.peakCoveredSNPFile = peakCoveredSNPFile;
        this.initialRho = initialRho;
        this.learningRate = learningRate;
        this.unImproveThreshold = unImproveThreshold;
        this.log = log;
        String outputFileName = peakCoveredSNPFile.substring(0, peakCoveredSNPFile.lastIndexOf("_")) + "_asePeak.txt";
        this.aseM6aPeakFile = new File(outputFileName);
    }

    /**
     * detect ASE m6A peaks and output to file
     */
    public void detectAsePeak() {
        // get global rho
        double rho = this.getGlobalRho();
        double epsilon = 0.00001;
        BufferedWriter bfw;
        try {
            bfw = new BufferedWriter(
                    new OutputStreamWriter(new FileOutputStream(this.aseM6aPeakFile))
            );
            HashMap<String, Double> peakPValues = this.getPValueOfEachPeak(rho);
            Set<Map.Entry<String, Double>> pVals = peakPValues.entrySet();
            for (Map.Entry<String, Double> pv: pVals) {
                System.out.println(pv.getKey() + "=>" + pv.getValue());
            }
            HashMap<String, Double> aseM6aPeaks = this.bhRecalibrationOfEachPeak(peakPValues);
            Set<Map.Entry<String, Double>> qValues = aseM6aPeaks.entrySet();
            double qVal;
            String peakLabel, writeOut;
            String[] info, outputLine;
            for (Map.Entry<String, Double> qValue: qValues) {
                peakLabel = qValue.getKey();
                info = peakLabel.split(":");
                String chrNum = info[0];
                String peakStart = info[1];
                String peakEnd = info[2];
                qVal = qValue.getValue();
                outputLine = new String[]{chrNum, peakStart, peakEnd, Double.toString(qVal)};
                writeOut = String.join("\t", outputLine);
                bfw.write(writeOut);
                bfw.newLine();
            }
            bfw.close();
        } catch (IOException ie) {
            this.log.error("can not write ASE m6a peak result file");
            this.log.error(ie.getMessage());
        }
    }

    /**
     * get major and minor haplotype SNP site reads count for each m6a peak
     * @return HashMap
     */
    private HashMap<String, HashMap<String, HashMap<String, LinkedList<Integer>>>> getPeakSNPReadsCount() {
        HeterozygoteReadsCount hrc = new HeterozygoteReadsCount(this.peakCoveredSNPFile, this.log);
        // chr: peakRange: maor/minor: readsCounts
        HashMap<String, HashMap<String, HashMap<String, LinkedList<Integer>>>> majorMinorHaplotype = hrc.getMajorMinorHaplotype();
        hrc = null;

        return majorMinorHaplotype;
    }

    /**
     * get major and minor haplotype SNP reads count of all SNP
     * @return HashMap
     */
    private HashMap<String, LinkedList<Integer>> getHaplotypeSNPReadsCount() {
        HaplotypeSNPReadsCount hsrc = new HaplotypeSNPReadsCount(this.peakCoveredSNPFile, log);
        HashMap<String, LinkedList<Integer>> haplotypeSNPReadsCount = hsrc.haplotypeSnpReadsCount();
        hsrc = null;

        return haplotypeSNPReadsCount;
    }

    /**
     * get the overdispersion value rho of the beta binomial distribution
     * @return rho
     */
    private double getGlobalRho() {
        HashMap<String, LinkedList<Integer>> haplotypeSNPReadsCount = this.getHaplotypeSNPReadsCount();
        LinkedList<Integer> majorHaplotype = haplotypeSNPReadsCount.get("major");
        int[] majorSNPReadsCount = new int[majorHaplotype.size()];
        for (int i = 0 ; i < majorHaplotype.size(); i++ ) {
            majorSNPReadsCount[i] = majorHaplotype.get(i);
        }
        majorHaplotype = null;
        LinkedList<Integer> minorHaplotype = haplotypeSNPReadsCount.get("minor");
        int[] minorSNPReadsCount = new int[minorHaplotype.size()];
        for (int i = 0 ; i < minorHaplotype.size(); i++ ) {
            minorSNPReadsCount[i] = minorHaplotype.get(i);
        }
        minorHaplotype = null;
        RhoEstimator re = new RhoEstimator(majorSNPReadsCount, minorSNPReadsCount, this.initialRho,
                                           this.learningRate, this.unImproveThreshold);
        double rho = re.gradientAscend();
        majorSNPReadsCount = null;
        minorSNPReadsCount = null;

        return rho;
    }

    /**
     * get p values for each m6a peak
     * @param rho beta binomial distribution parameter
     * @return HashMap
     */
    private HashMap<String, Double> getPValueOfEachPeak(double rho) {

        double pVal;
        SignificantTest st;
        HashMap<String, Double> pValues = new HashMap<>();

        HashMap<String, HashMap<String, HashMap<String, LinkedList<Integer>>>> peakSnpReadsCount = this.getPeakSNPReadsCount();
        Set<Map.Entry<String, HashMap<String, HashMap<String, LinkedList<Integer>>>>> chrPeaks = peakSnpReadsCount.entrySet();
        for (Map.Entry<String, HashMap<String, HashMap<String, LinkedList<Integer>>>> chrPeak: chrPeaks) {
            String chrNum = chrPeak.getKey();
            // peaks on a particular chromosome
            HashMap<String, HashMap<String, LinkedList<Integer>>> peaksOnChr = chrPeak.getValue();

            Set<Map.Entry<String, HashMap<String, LinkedList<Integer>>>> peakRanges = peaksOnChr.entrySet();
            for (Map.Entry<String, HashMap<String, LinkedList<Integer>>> peakRange: peakRanges) {
                // the format of key => start:end
                String peakStartToEnd = peakRange.getKey();
                // major and minor haplotype SNP reads in a certain peak range
                HashMap<String, LinkedList<Integer>> haplotypeSnpReads = peakRange.getValue();

                LinkedList<Integer> majorSNPReadsCount = haplotypeSnpReads.get("major");
                int[] majorCount = new int[majorSNPReadsCount.size()];
                for (int i = 0 ; i < majorSNPReadsCount.size(); i++ ) {
                    majorCount[i] = majorSNPReadsCount.get(i);
                }
                majorSNPReadsCount = null;
                LinkedList<Integer> minorSNPReadsCount = haplotypeSnpReads.get("minor");
                int[] minorCount = new int[minorSNPReadsCount.size()];
                for (int i = 0 ; i < minorSNPReadsCount.size(); i++ ) {
                    minorCount[i] = minorSNPReadsCount.get(i);
                }
                minorSNPReadsCount = null;
                // meta analysis for all SNP under a peak, get p value
                st = new SignificantTest(majorCount, minorCount, rho);
                pVal = st.testSignificant();
                pValues.put(chrNum+":"+peakStartToEnd, pVal);
            }
        }

        return pValues;
    }

    /**
     * recalibrate peak p values using BH method, output Q value. Select Q values < 0.05 as ASE m6a peak
     * @param peakPValues p values for each m6a peak
     * @return HashMap
     */
    private HashMap<String, Double> bhRecalibrationOfEachPeak(HashMap<String, Double> peakPValues) {
        List<Map.Entry<String, Double>> sortedByValue = new ArrayList<Map.Entry<String, Double>>(peakPValues.entrySet());
        // sort peak p values from small to large
        Collections.sort(sortedByValue, new Comparator<Map.Entry<String, Double>>() {
            public int compare(Map.Entry<String, Double> o1,
                               Map.Entry<String, Double> o2) {
                return (o1.getValue()).compareTo(o2.getValue());
            }
        });

        HashMap<String, Double> aseM6aPeaks = new HashMap<>();
        int rankage = 0;
        int totalPeak = sortedByValue.size();
        double peakPValue;
        String peakLabel;
        for (int i = 0; i < totalPeak; i++) {
            rankage += 1;
            peakLabel = sortedByValue.get(i).getKey();
            peakPValue = sortedByValue.get(i).getValue();
            double qValue = SignificantTest.BHRecalibration(peakPValue, rankage, totalPeak);
            if ((qValue - 0.05) > 0.00001)
                break;
            aseM6aPeaks.put(peakLabel, qValue);
        }

        return aseM6aPeaks;
    }
}
