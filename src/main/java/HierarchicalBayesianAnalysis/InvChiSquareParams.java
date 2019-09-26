package HierarchicalBayesianAnalysis;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.special.Gamma;

import java.io.*;
import java.util.ArrayList;
import java.util.TreeMap;

public class InvChiSquareParams {
    private String vcfFile;
    private int readsCoverageThreshold;
    public ArrayList<Double> lorList = new ArrayList<>();
    public Double lorMean, lorStd, infimum, supremum;
    private double degreeOfFreedom;
    private ChiSquaredDistribution csd;

    public InvChiSquareParams(String vcfFile, int readsCoverageThreshold, double degreeOfFreedom) {
        this.vcfFile = vcfFile;
        this.readsCoverageThreshold = readsCoverageThreshold;
        this.degreeOfFreedom = degreeOfFreedom;
        this.csd = new ChiSquaredDistribution(degreeOfFreedom);
    }

    public InvChiSquareParams(double lorStd, double degreeOfFreedom) {
        this.lorStd = lorStd;
        this.degreeOfFreedom = degreeOfFreedom;
        this.csd = new ChiSquaredDistribution(degreeOfFreedom);
    }

    public void paramEstimate() {
        this.parseVcfFile();
        this.calcLorMean();
        this.calcLorStd();
        this.calcInfimumAndSupremum();
    }

    public void parseVcfFile() {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.vcfFile))));
            String line = "";
            String[] info, record;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    Double lor = null;
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    record = info[7].split(";");
                    for (String re: record) {
                        if (re.startsWith("DP4"))
                            lor = this.calcLOR(re);
                    }
                    if (lor != null)
                        this.lorList.add(lor);
                }
            }
        } catch (IOException ie) {
            ie.printStackTrace();
        } finally {
            if (bfr != null) {
                try {
                    bfr.close();
                } catch (IOException ie) {
                    ie.printStackTrace();
                }
            }
        }
    }

    public void calcLorMean() {
        Double sum = 0.0;
        for (Double lor: this.lorList) {
            sum += lor;
        }

        this.lorMean = sum / this.lorList.size();
    }

    public void calcLorStd() {
        double variance = 0.0;
        for (Double lor: this.lorList) {
            variance += Math.pow((lor - this.lorMean), 2);
        }

        this.lorStd = Math.sqrt(variance / this.lorList.size());
    }

    public void calcInfimumAndSupremum() {
        this.infimum = Double.MAX_VALUE;
        this.supremum = Double.MIN_VALUE;
        for (Double lor: this.lorList) {
            double distance = Math.abs((lor - this.lorMean));
            if (distance < this.infimum)
                this.infimum = distance;
            if (distance > this.supremum)
                this.supremum = distance;
        }
    }

    private Double calcLOR(String dp4) {
        String[] readsCount = dp4.split("=")[1].split(",");
        int refCount = Integer.valueOf(readsCount[0]) + Integer.valueOf(readsCount[1]);
        int altCount = Integer.valueOf(readsCount[2]) + Integer.valueOf(readsCount[3]);
        if (refCount == 0 || altCount == 0)
            return null;
        if (Math.max(refCount, altCount) < this.readsCoverageThreshold)
            return null;
        int majorAlleleCount = (refCount >= altCount)? refCount: altCount;
        int minorAlleleCount = (refCount < altCount)? refCount: altCount;

        return Math.log((double) majorAlleleCount/(double) minorAlleleCount);
    }

    public TreeMap<Double, Integer> invChiSquareSampling(int samplingTime) {
        TreeMap<Double, Integer> sampleFrequency = new TreeMap<>();

        for (int i=0; i<samplingTime; i++) {
            double val = this.degreeOfFreedom * Math.pow(this.lorStd, 2) / this.csd.sample();
            int frequency = sampleFrequency.getOrDefault(val, 0);
            frequency += 1;
            sampleFrequency.put(val, frequency);
        }

        return sampleFrequency;
    }

    public double sample() {
        return this.degreeOfFreedom * Math.pow(this.lorStd, 2) / this.csd.sample();
    }

    public boolean isGoodParam(TreeMap<Double, Integer> sampleFrequency, int samplingTime) {
        int confidenceIntervalInfimum = (int) (samplingTime * 0.025);
        int confidenceIntervalSupremum = (int) (samplingTime * 0.975);
        int cum = 0;
        double sampleLorSupremum = Double.MAX_VALUE, sampleLorInfimum = Double.MIN_VALUE;
        for (Object lor: sampleFrequency.descendingKeySet()) {
            cum += sampleFrequency.get(lor);
            if (cum > confidenceIntervalInfimum)
                sampleLorSupremum = (double) lor;
            if (cum > confidenceIntervalSupremum)
                sampleLorInfimum = (double) lor;
        }
        boolean b1 = false, b2 = false, b3 = false;
        if (sampleLorInfimum >= this.infimum * 0.95)
            b1 = true;
        if (sampleLorSupremum <= this.supremum * 1.05)
            b2 = true;
        double sampleMedian = this.getSampleMedian(sampleFrequency, samplingTime);
        if (sampleMedian >= this.lorStd*0.95 && sampleMedian <= this.lorStd*1.05)
            b3 = true;

        return (b1 && b2) && b3;
    }

    public Double getSampleMedian(TreeMap<Double, Integer> sampleFreq, int samplingTime) {
        double cum = 0;
        for (Double val: sampleFreq.keySet()) {
            cum += sampleFreq.get(val);
            if (cum >= samplingTime / 2)
                return val;
        }

        return null;
    }

    public void getIntervalBoundary(TreeMap<Double, Integer> sampleFreq, int samplingTime) {
        int start = (int) (samplingTime * 0.025), end = (int) (samplingTime*0.975), cum = 0;
        double low = -1, high = -1;
        for (Double val: sampleFreq.descendingKeySet()) {
            cum += sampleFreq.get(val);
            if (cum >= start) {
                high = val;
                start = Integer.MAX_VALUE;
            }
            if (cum >= end) {
                low = val;
                System.out.println("low: " + low + "\thigh: " + high);
                break;
            }
        }
    }

    public double density(double x) {
        return Math.pow(2, this.degreeOfFreedom/2) / Gamma.gamma(this.degreeOfFreedom/2) * Math.pow(x, -this.degreeOfFreedom/2-1) * Math.exp(-1/(2*x));
    }

    public void output(File outputFile, TreeMap<Double, Integer> sampleFreq) {
        BufferedWriter bfw = null;
        try {
            bfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outputFile)));
            bfw.write("LOR\tfrequency\n");
            for (Object obj: sampleFreq.descendingKeySet()) {
                double lor = (double) obj;
                String[] record = new String[] {String.valueOf(lor), String.valueOf(sampleFreq.get(lor))};
                bfw.write(String.join("\t", record));
                bfw.newLine();
            }
        } catch (IOException ie) {
            ie.printStackTrace();
        } finally {
            if (bfw != null) {
                try {
                    bfw.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }
}
