package HierarchicalBayesianAnalysis;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.special.Gamma;

import java.io.*;
import java.util.TreeMap;

public class InvChiSquareParams {
    public double lorStd, degreeOfFreedom;
    private ChiSquaredDistribution csd;

    public InvChiSquareParams(double lorStd, double degreeOfFreedom) {
        this.lorStd = lorStd;
        this.degreeOfFreedom = degreeOfFreedom;
        this.csd = new ChiSquaredDistribution(degreeOfFreedom);
    }

    public double sample() {
        return this.degreeOfFreedom * Math.pow(this.lorStd, 2) / this.csd.sample();
    }

    public double density(double x) {
        return Math.pow(2, this.degreeOfFreedom/2) / Gamma.gamma(this.degreeOfFreedom/2) * Math.pow(x, -this.degreeOfFreedom/2-1) * Math.exp(-1/(2*x));
    }

    public double logDensity(double x) {
        return this.degreeOfFreedom/2 * Math.log(2) - Gamma.logGamma(this.degreeOfFreedom/2) + (-this.degreeOfFreedom/2-1) * Math.log(x) + (-1 / (2*x));
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
