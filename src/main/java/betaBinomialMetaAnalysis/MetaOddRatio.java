package betaBinomialMetaAnalysis;


public class MetaOddRatio {

    private int[] majorSNPReads, minorSNPReads;
    private double rho;

    /**
     * Constructor
     * @param majorSNPReads number of reads cover major SNP sites under a m6A peak
     * @param minorSNPReads number of reads cover minor SNP sites under a m6A peak
     * @param rho parameter of beta-binomial distribution
     */
    public MetaOddRatio(int[] majorSNPReads, int[] minorSNPReads, double rho) {
        this.majorSNPReads = majorSNPReads;
        this.minorSNPReads = minorSNPReads;
        this.rho = rho;
    }

    /**
     * combine all the major SNP under a m6A peak to calculate the combined odd ratio
     * @return double array contains meta odd ratio、log meta odd ratio and variance of log meta odd ratio
     */
    public double[] metaAnalysisOddRatio() {

        double[] calibratedSingleSNPOddRatio = new double[this.majorSNPReads.length];
        double[] singleSNPWeight = new double[this.majorSNPReads.length];
        double[] singleSNPVariant = new double[this.majorSNPReads.length];
        double singleOR, singleWeight, singleVariant, calibratedOddRatio;

        for (int i = 0; i < this.majorSNPReads.length; i++) {
            singleOR = oddRatioSingleSNP(this.majorSNPReads[i], this.minorSNPReads[i]);
            singleWeight = oddRatioWeightSingleSNP(this.majorSNPReads[i], this.minorSNPReads[i]);
            singleVariant = oddRatioVarianceSingleSNP(this.majorSNPReads[i], this.minorSNPReads[i]);

            calibratedOddRatio = singleWeight * singleOR;
            calibratedSingleSNPOddRatio[i] = calibratedOddRatio;
            singleSNPWeight[i] = singleWeight;
            singleSNPVariant[i] = singleVariant;
        }

        double weightSum = getSum(singleSNPWeight);
        double calibratedORSum = getSum(calibratedSingleSNPOddRatio);
        double variantSum = getSum(singleSNPVariant);

        double metaOddRatio = calibratedORSum / weightSum;
        double logMetaOddRatio = Math.log(metaOddRatio);

        return new double[]{metaOddRatio, logMetaOddRatio, variantSum};
    }

    /**
     * odd ratio for a single SNP site
     * @param majorSNPReads number of reads covered by a single major SNP site
     * @param minorSNPReads number of reads covered by a single minor SNP site
     * @return odd ratio for a single SNP site
     */
    private double oddRatioSingleSNP(int majorSNPReads, int minorSNPReads) {
        double majorSNPProbability = (double)majorSNPReads / (double)(majorSNPReads + minorSNPReads);
        return majorSNPProbability / (1 - majorSNPProbability);
    }

    /**
     * variance of odd ratio for a single SNP site
     * @param majorSNPReads number of reads covered by a single major SNP site
     * @param minorSNPReads number of reads covered by a single minor SNP site
     * @return variance of odd ratio
     */
    private double oddRatioVarianceSingleSNP(int majorSNPReads, int minorSNPReads) {
        double majorSNPProbability = (double)majorSNPReads / (double)(majorSNPReads + minorSNPReads);
        double minorSNPProbability = 1 - majorSNPProbability;
        int totalReads = majorSNPReads + minorSNPReads;
        return (1 + (totalReads - 1) * this.rho) / (totalReads * majorSNPProbability * minorSNPProbability) + (1 + (totalReads - 1) * this.rho) / (0.25 * totalReads);
    }

    /**
     * weight of odd ratio for a single SNP site
     * @param majorSNPReads number of reads covered by a single major SNP site
     * @param minorSNPReads number of reads covered by a single minor SNP site
     * @return weight of odd ratio
     */
    private double oddRatioWeightSingleSNP(int majorSNPReads, int minorSNPReads) {
        double majorSNPProbability = (double)majorSNPReads / (double)(majorSNPReads + minorSNPReads);
        int totalReads = majorSNPReads + minorSNPReads;
        double calibrator = 1 + (totalReads) * this.rho;
        return 1 / (2 * calibrator / totalReads) * 0.5 * (1 - majorSNPProbability);
    }

    /**
     * get sum of the SNP data
     * @param SNPData Array contains SNP data
     * @return sum result
     */
    private double getSum(double[] SNPData) {
        double sum = 0;
        for (double data : SNPData) {
            sum += data;
        }
        return sum;
    }

}
