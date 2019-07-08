package HierarchicalBayesianAnalysis;

import org.apache.commons.math3.distribution.UniformRealDistribution;

public class MHSampling {
    private UniformRealDistribution urd;

    public MHSampling() {
        urd = new UniformRealDistribution(0, 1);
    }

    /**
     * 根据随机抽取的u和接受率的值决定本轮采样结果
     * @param curSamplingVal 本轮采样得到的样本值
     * @param curSamplingDensity 本轮采样样本值对应的概率密度
     * @param prevSamplingVal 上一轮采样得到的样本值
     * @param prevSamplingDensity 上一轮采样样本值对应的概率密度
     * @return 本轮采样结果
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
     * 随机在均匀分布(0, 1)中抽取一个值
     * @return uniform(0, 1)中的值
     */
    protected double getRandomVal() {
        return urd.sample();
    }

    /**
     * 计算MH采样的接受率
     * @return 接受率
     */
    protected double getReceptance(double curSamplingDensity, double prevSamplingDensity) {
        return Math.min(1, curSamplingDensity/prevSamplingDensity);
    }
}
