package HierarchicalBayesianAnalysis;

import org.apache.commons.math3.distribution.UniformRealDistribution;

public class TauSampling extends MHSampling {
    private UniformRealDistribution priorTau;

    /**
     * Constructor
     * @param low Tau先验分布服从的均匀分布的下界
     * @param high Tau先验分布服从的均匀分布的上界
     */
    public TauSampling(double low, double high) {
        super();
        this.priorTau = new UniformRealDistribution(low, high);
    }

    /**
     * 获取新一轮采样最终结果
     * @param prevTau 当前的tau
     * @param prevTauDensity 当前tau对应的后验概率
     * @param logOddRatios 当前各个ASE位点major allele reads对数优势比
     * @param variances 当前各个ASE位点major allele reads对数优势比的方差
     * @param miu 当前全局对数优势比的均值
     * @param sigma 当前全局对数优势比的方差
     * @return 采样结果及该结果对应的后验概率
     */
    public double[] sampling(double prevTau, double prevTauDensity, double[] logOddRatios, double[] variances,
                             double miu, double sigma) {
        double curTau = this.randomTau();
        double curTauPosteriorDensity = this.posteriorTau(curTau, logOddRatios, variances, miu, sigma);
        System.out.println("random Tau: " + curTau + "\tposterior density: " + curTauPosteriorDensity);

        return this.getSamplingRes(curTau, curTauPosteriorDensity, prevTau, prevTauDensity);
    }

    /**
     * 从先验分布中随机抽取一个值作为tau的值作为新一轮采样值
     * @return 随机抽取的tau的值
     */
    public double randomTau() {
        return this.priorTau.sample();
    }

    /**
     * 随机采样的新一轮的tau在先验分布中对应的概率
     * @param prevSamplingTau 上一轮采样得到的tau的值
     * @return 对应的先验概率密度
     */
    public double priorTauDensity(double prevSamplingTau) {
        return this.priorTau.density(prevSamplingTau);
    }

    /**
     * 随机抽取的新一轮采样tau的值对应的后验概率近似值
     * @param curSamplingTau 新一轮采样的tau值
     * @param logOddRatios 当前各个ASE位点major allele reads对数优势比
     * @param variances 当前各个ASE位点major allele reads对数优势比的方差
     * @param miu 当前全局对数优势比的均值
     * @param sigma 当前全局对数优势比的方差
     * @return 后验概率近似值
     */
    public double posteriorTau(double curSamplingTau, double[] logOddRatios, double[] variances, double miu, double sigma) {
        double cumProd = 1;

        for (int i=0; i<variances.length; i++) {
            cumProd = cumProd * Math.pow(Math.pow(curSamplingTau, 2) + variances[i], -0.5) * Math.exp(-1 * Math.pow(logOddRatios[i]-miu, 2)/(2*(variances[i]+Math.pow(curSamplingTau, 2))));
        }

        double curTauPriorDensity = this.priorTauDensity(curSamplingTau);

        return curTauPriorDensity * Math.pow(sigma, 0.5) * cumProd;
    }
}
