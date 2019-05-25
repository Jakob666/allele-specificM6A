package betaBinomialMetaAnalysis;

import org.apache.commons.math3.distribution.NormalDistribution;

public class SignificantTest {

    private MetaOddRatio mor;

    public SignificantTest(int[] majorPeakReads, int[] minorPeakReads, double rho) {
        this.mor = new MetaOddRatio(majorPeakReads, minorPeakReads, rho);
    }

    /**
     * calculate p value with the result of meta analysis
     * @return significant test p value
     */
    public double testSignificant() {
        double[] metaValues = this.mor.metaAnalysisOddRatio();
        double logMetaOddRatio = metaValues[1];
        double variantSum = metaValues[2];

        return cumulativeProbability(logMetaOddRatio, 1.0 / variantSum);
    }

    /**
     * BH method for recalibrating p values of significant test and get q value
     * @param pValues p value from significant test
     * @param rank the rank of the values corresponding to the peak, sorted from smallest to larges
     * @param peakNum total number of m6A peaks on haplotype
     * @return recalibrated q value
     */
    public static double BHRecalibration(double pValues, int rank, int peakNum) {
        return Math.min(pValues * peakNum / rank, 1.0);
    }

    /**
     * calculate cumulative probability
     * @param miu mean value of normal distribution
     * @param sigma standard deviation of normal distribution
     * @return cumulative probability
     */
    private double cumulativeProbability(double miu, double sigma) {
        double cumProba;
        NormalDistribution normalDistribution = new NormalDistribution(miu, sigma);

        if (miu > 0) {
            cumProba = 2 * normalDistribution.cumulativeProbability(0);
        } else {
            cumProba = 2 * (1 - normalDistribution.cumulativeProbability(0));
        }

        return cumProba;
    }

}
