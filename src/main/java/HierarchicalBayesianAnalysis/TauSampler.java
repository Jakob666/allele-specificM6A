package HierarchicalBayesianAnalysis;


import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.Arrays;

public class TauSampler extends MHSampling {
    private ScaledInvChiSquareDistribution priorTau;
    private NormalDistribution nd;

    /**
     * Constructor
     * @param df number of chi-squared degrees of freedom
     * @param scaleParam scaling parameter
     */
    public TauSampler(double df, double scaleParam) {
        super();
        this.priorTau = new ScaledInvChiSquareDistribution(df, scaleParam);
        this.nd = new NormalDistribution(0, 0.1);
    }

    /**
     * sampling a new tau from its priority distribution
     * @return randomly sample tau
     */
    public double randomInit() {
        return this.priorTau.sample();
    }

    /**
     * sampling a new tau from proposal distribution
     * @return randomly sample tau
     */
    public double randomTau(double prevTau) {
        return Math.abs(prevTau + this.nd.sample());
    }

    /**
     * logarithm tau posterior density, more robust than posteriorTau method
     * @param curSamplingTau tau
     * @param logOddRatios observed LOR, correspond to vector y
     * @param variances observed variance, correspond to vector sigma^2
     * @param globalLORMean global LOR mean, correspond to u^
     * @param globalLORVar global LOR Var, correspond to Vu
     * @return logarithm tau posterior density
     */
    public double logPosteriorTau(double curSamplingTau, double[] logOddRatios, double[] variances, double globalLORMean, double globalLORVar) {
        double sum = 0, std, frac;
        for (int i=0; i<variances.length; i++) {
            std = variances[i] + Math.pow(curSamplingTau, 2);
            frac = Math.pow(logOddRatios[i] - globalLORMean, 2) / 2 / std;
            sum -= 0.5 * Math.log(std);
            sum -= frac;
        }

        return sum + this.priorTau.logDensity(1/curSamplingTau) + 0.5 * Math.log(globalLORVar);
    }
}
