package HierarchicalBayesianAnalysis;

import org.apache.commons.math3.distribution.GammaDistribution;

public class ScaledInvChiSquareDistribution {
    private GammaDistribution gd;

    /**
     * Constructor
     * if X ~ scaled-inv-Chi(v, t^2), then X ~ inv-gamma(v/2, vt^2/2)
     * if X ~ gamma(k, theta), then 1/X ~ inv-gamma(k, 1 / theta) => if X ~ gamma(v/2, 2/vt^2), then 1/X ~ inv-gamma(v/2, vt^2/2)
     * @param v number of chi-squared degrees of freedom, large than 0
     * @param tauSquare scaling parameter, large than 0
     */
    public ScaledInvChiSquareDistribution(double v, double tauSquare) {
        double gammaShape = v * 0.5;
        double gammaScale = 2 / v / tauSquare;
        this.gd = new GammaDistribution(gammaShape, gammaScale);
    }

    public double sample() {
        return Math.pow(this.gd.sample(), -1);
    }

    public double logDensity(double x) {
        return this.gd.logDensity(1/x);
    }
}
