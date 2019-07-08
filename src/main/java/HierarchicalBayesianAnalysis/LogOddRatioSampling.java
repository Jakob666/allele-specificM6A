package HierarchicalBayesianAnalysis;

import org.apache.commons.math3.distribution.NormalDistribution;

public class LogOddRatioSampling {
    public LogOddRatioSampling() {}

    /**
     * 通过新一轮采样得到的tau，获取新一轮的全局对数优势比采样值
     * @param curTau 新一轮采样得到的tau
     * @param logOddRatios 现有的ASE位点的对数优势比
     * @param variances 现有的ASE位点的对数优势比的方差
     * @return 新一轮采样的对数优势比
     */
    public double globalLogOddRatioSampling(double curTau, double[] logOddRatios, double[] variances) {
        double miu = 0, v = 0;
        for (int i=0; i<logOddRatios.length; i++) {
            miu += 1.0 / (variances[i] + Math.pow(curTau, 2)) * logOddRatios[i];
            v += 1.0 / (variances[i] + Math.pow(curTau, 2));
        }
        miu = miu * v;
        v = 1.0 / v;

        NormalDistribution nd = new NormalDistribution(miu, v);

        return nd.sample();
    }

    /**
     * 通过新一轮采样得到的tau和全局对数优势比采样值，采样每个ASE位点的对数优势比
     * @param curTau 新一轮采样得到的tau
     * @param curGlobalOddRatio 新一轮采样得到的全局对数优势比
     * @param logOddRatios 现有的ASE位点的对数优势比
     * @param variances 现有的ASE位点的对数优势比的方差
     * @return ASE位点新一轮对数优势比的采样值
     */
    public double[] singleAseOddRatioSampling(double curTau, double curGlobalOddRatio, double[] logOddRatios,
                                            double[] variances) {
        double miu = 0, v = 0;
        NormalDistribution nd;
        double[] newAseLogOddRatio = new double[logOddRatios.length];
        for (int i=0; i<logOddRatios.length; i++) {
            miu = (1.0/variances[i] * logOddRatios[i] + 1.0/Math.pow(curTau, 2) * curGlobalOddRatio) / (1.0/variances[i] + 1.0/Math.pow(curTau, 2));
            v = 1.0 / (1.0/variances[i] + 1.0/Math.pow(curTau, 2));
            nd = new NormalDistribution(miu, v);
            double LORMean = nd.sample();
            nd = new NormalDistribution(LORMean, variances[i]);
            newAseLogOddRatio[i] = nd.sample();
        }

        return newAseLogOddRatio;
    }
}
