package HierarchicalBayesianAnalysis;

import org.apache.commons.math3.distribution.UniformRealDistribution;

public class MHSampling {
    private UniformRealDistribution urd;

    public MHSampling() {
        urd = new UniformRealDistribution(0, 1);
    }

    /**
     * get new round sampling result via random u and receptance
     * @param curSamplingVal sampling result of new round
     * @param curSamplingDensity posterior density of sampling result of new round
     * @param prevSamplingVal sampling result of last round
     * @param prevSamplingDensity posterior density of sampling result of last round
     * @return sampling result
     */
    public double[] getSamplingRes(double curSamplingVal, double curSamplingDensity,
                                 double prevSamplingVal, double prevSamplingDensity) {
        double u, receptance;
        u = this.getRandomVal();
        receptance = this.getReceptance(curSamplingDensity, prevSamplingDensity);
        if (u < receptance)
            return new double[]{curSamplingVal, curSamplingDensity};
        else
            return new double[]{prevSamplingVal, prevSamplingDensity};
    }

    /**
     * randomly sample from uniform(0, 1)
     * @return sampling result
     */
    protected double getRandomVal() {
        return urd.sample();
    }

    /**
     * calculate MH sampling receptance
     * @return receptance
     */
    protected double getReceptance(double curSamplingDensity, double prevSamplingDensity) {
        return Math.min(1, curSamplingDensity/prevSamplingDensity);
    }
}
