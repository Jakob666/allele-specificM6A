package HierarchicalBayesianAnalysis;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;

import java.util.HashMap;

public class LogOddRatioSampler {
    private UniformRealDistribution urd;

    public LogOddRatioSampler() {
        this.urd = new UniformRealDistribution(0, 1.0);
    }

    /**
     * start a new round and sample for tau, globalLOR
     * @param curTau current tau
     * @param logOddRatios current LOR for SNV sites
     * @param variances current variance for SNV sites
     * @param prevGlobalLOR globalLOR in last round
     * @param prevGlobalLORPosteriorDensity the posterior probability of last round globalLOR
     * @return sampling globalLOR, variance value in the new round
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
     * start a new round and sample for SNV sites expected LOR
     * @param curTau current tau
     * @param curGlobalOddRatio current globalLOR
     * @param logOddRatios LOR in last round
     * @param variances LOR variance in last round
     * @return sampling result in new round
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
