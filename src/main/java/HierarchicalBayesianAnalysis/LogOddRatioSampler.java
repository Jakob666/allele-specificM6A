package HierarchicalBayesianAnalysis;

import org.apache.commons.math3.distribution.NormalDistribution;

public class LogOddRatioSampler extends MHSampling {
    private NormalDistribution nd;

    public LogOddRatioSampler() {
        this.nd = new NormalDistribution(0, 0.1);
    }

    public double randomU(double curU) {
        return curU + this.nd.sample();
    }

    public double logPosteriorObserveLOR(double observedLOR, double expectedLOR, double variance) {
        NormalDistribution nd = new NormalDistribution(expectedLOR, Math.sqrt(variance));
        return nd.logDensity(observedLOR);
    }

    public double logPosteriorExpectedLOR(double observedLOR, double expectedLOR, double variance, double t, double u) {
        double numerator = 1 / variance * observedLOR + 1 / Math.pow(t, 2) * u;
        double denominator = 1 / variance + 1 / Math.pow(t, 2);
        double thetaMean = numerator / denominator;
        double thetaVar = 1 / denominator;
        NormalDistribution nd = new NormalDistribution(thetaMean, Math.sqrt(thetaVar));
        return nd.logDensity(expectedLOR);
    }
}
