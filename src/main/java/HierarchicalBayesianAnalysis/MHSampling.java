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
     * @param curSamplingLogDensity posterior density of sampling result of new round
     * @param prevSamplingVal sampling result of last round
     * @param prevSamplingLogDensity posterior density of sampling result of last round
     * @return sampling result
     */
    public double[] getSamplingRes(double curSamplingVal, double curSamplingLogDensity,
                                   double prevSamplingVal, double prevSamplingLogDensity) {
        double u, receptance;
        u = this.getRandomVal();
        receptance = this.getReceptance(curSamplingLogDensity, prevSamplingLogDensity);
        if (u < receptance)
            return new double[]{curSamplingVal, curSamplingLogDensity};
        else
            return new double[]{prevSamplingVal, prevSamplingLogDensity};
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
    protected double getReceptance(double curSamplingLogDensity, double prevSamplingLogDensity) {
        return Math.min(1, Math.exp(curSamplingLogDensity - prevSamplingLogDensity));
    }
}
