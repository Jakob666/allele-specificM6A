package HierarchicalBayesianAnalysis;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;

import java.io.*;
import java.util.HashMap;
import java.util.HashSet;

public class LogOddRatioSampling {
    private UniformRealDistribution urd;

    public LogOddRatioSampling() {
        this.urd = new UniformRealDistribution(0, 1.0);
    }

    /**
     * 通过新一轮采样得到的tau，获取新一轮的全局对数优势比均值、方差及采样值
     * @param curTau 新一轮采样得到的tau
     * @param logOddRatios 现有的ASE位点的对数优势比
     * @param variances 现有的ASE位点的对数优势比的方差
     * @param prevGlobalLOR 上一轮抽样的全局对数优势比的值
     * @param prevGlobalLORPosteriorDensity 上一轮抽样的全局对数优势比的后验概率
     * @return 新一轮采样的对数优势比均值、方差及采样值
     */
    public double[] globalLogOddRatioSampling(double curTau, double[] logOddRatios, double[] variances,
                                              double prevGlobalLOR, double prevGlobalLORPosteriorDensity) {
        double miu = 0, v = 0;
        for (int i=0; i<logOddRatios.length; i++) {
            miu += 1.0 / (variances[i] + Math.pow(curTau, 2)) * logOddRatios[i];
            v += 1.0 / (variances[i] + Math.pow(curTau, 2));
        }
        miu = miu / v;
        v = 1.0 / v;

        NormalDistribution nd = new NormalDistribution(miu, v);
        double sampleGlobalLOR = nd.sample();
        double sampleGlobalLORDensity = nd.density(sampleGlobalLOR);
        double receptance = Math.min(1.0, sampleGlobalLORDensity/prevGlobalLORPosteriorDensity);

        double randNum = this.urd.sample();
        if (randNum < receptance)
            return new double[]{miu, v, sampleGlobalLOR, sampleGlobalLORDensity};
        else
            return new double[]{miu, v, prevGlobalLOR, prevGlobalLORPosteriorDensity};
    }

    /**
     * 通过新一轮采样得到的tau和全局对数优势比采样值，采样每个ASE位点的对数优势比均值
     * @param curTau 新一轮采样得到的tau
     * @param curGlobalOddRatio 新一轮采样得到的全局对数优势比
     * @param logOddRatios 现有的ASE位点的对数优势比
     * @param variances 现有的ASE位点的对数优势比的方差
     * @return ASE位点新一轮对数优势比的采样值
     */
    public HashMap<String, double[]> singleAseOddRatioMeanSampling(double curTau, double curGlobalOddRatio, double[] logOddRatios,
                                            double[] variances, double[] prevLOR, double[] prevLORPosteriorDensity) {
        double miu = 0, v = 0, prevSampleLOR, prevSampleLORPosteriorDensity, sampleLOR, sampleLORPosteriorDensity,
                randNum, receptance;
        NormalDistribution nd;
        double[] newAseLogOddRatio = new double[logOddRatios.length], newAseLORPosteriorDensity = new double[logOddRatios.length];
        for (int i=0; i<logOddRatios.length; i++) {
            miu = (1.0/variances[i] * logOddRatios[i] + 1.0/Math.pow(curTau, 2) * curGlobalOddRatio) / (1.0/variances[i] + 1.0/Math.pow(curTau, 2));
            v = 1.0 / (1.0/variances[i] + 1.0/Math.pow(curTau, 2));
            nd = new NormalDistribution(miu, v);

            prevSampleLOR = prevLOR[i];
            prevSampleLORPosteriorDensity = prevLORPosteriorDensity[i];
            sampleLOR = nd.sample();
            sampleLORPosteriorDensity = nd.density(sampleLOR);

            randNum = this.urd.sample();
            receptance = Math.min(1.0, sampleLORPosteriorDensity / prevSampleLORPosteriorDensity);
            if (randNum < receptance) {
                newAseLogOddRatio[i] = sampleLOR;
                newAseLORPosteriorDensity[i] = sampleLORPosteriorDensity;
            } else {
                newAseLogOddRatio[i] = prevSampleLOR;
                newAseLORPosteriorDensity[i] = prevSampleLORPosteriorDensity;
            }
        }

        HashMap<String, double[]> sampleRes = new HashMap<>();
        sampleRes.put("LOR", newAseLogOddRatio);
        sampleRes.put("density", newAseLORPosteriorDensity);

        return sampleRes;
    }

    public HashMap<String, double[]> singleAseOddRatioSampling(double[] singleAseLORMean, double[] variance,
                                                               double[] prevSingleAseLOR, double[] prevSingleAseLORPosteriorDensity) {
        double[] newSingleAseLOR = new double[singleAseLORMean.length];
        double[] newSingleAseLORPosteriorDensity = new double[singleAseLORMean.length];
        for (int j=0; j < singleAseLORMean.length; j++) {
            NormalDistribution nd = new NormalDistribution(singleAseLORMean[j], variance[j]);
            double sampleLor = nd.sample();
            double sampleLorPosteriorDensity = nd.density(sampleLor);
            double receptance = Math.min(1.0, sampleLorPosteriorDensity/prevSingleAseLORPosteriorDensity[j]);
            double randNum = this.urd.sample();
            if (randNum < receptance) {
                newSingleAseLOR[j] = sampleLor;
                newSingleAseLORPosteriorDensity[j] = nd.density(sampleLor);
            } else {
                newSingleAseLOR[j] = prevSingleAseLOR[j];
                newSingleAseLORPosteriorDensity[j] = prevSingleAseLORPosteriorDensity[j];
            }
        }
        HashMap<String, double[]> sampleRes = new HashMap<>();
        sampleRes.put("LOR", newSingleAseLOR);
        sampleRes.put("density", newSingleAseLORPosteriorDensity);

        return sampleRes;
    }
}
