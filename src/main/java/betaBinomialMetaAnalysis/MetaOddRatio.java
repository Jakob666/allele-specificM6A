package betaBinomialMetaAnalysis;


public class MetaOddRatio {

    private int[] majorSNPReads, minorSNPReads;
    private double rho;

    /**
     * Constructor
     * @param majorSNPReads 某个m6A信号下各 major SNP位点上reads的数目
     * @param minorSNPReads 某个m6A信号下各 minor SNP位点上reads的数目
     * @param rho beta-binomial二项分布全局离散度
     */
    public MetaOddRatio(int[] majorSNPReads, int[] minorSNPReads, double rho) {
        this.majorSNPReads = majorSNPReads;
        this.minorSNPReads = minorSNPReads;
        this.rho = rho;
    }

    /**
     * 将一个m6A peak下所有的major SNP的优势比整合计算整个peak的优势比
     * @return [优势比，对数优势比，各位点方差]
     */
    public double[] metaAnalysisOddRatio() {

        // 3个数组分别记录每个位点的优势比、权重及方差
        double[] calibratedSingleSNPOddRatio = new double[this.majorSNPReads.length];
        double[] singleSNPWeight = new double[this.majorSNPReads.length];
        double[] singleSNPVariant = new double[this.majorSNPReads.length];
        double singleOR, singleWeight, singleVariant, calibratedOddRatio;

        for (int i = 0; i < this.majorSNPReads.length; i++) {
            singleOR = oddRatioSingleSNP(this.majorSNPReads[i], this.minorSNPReads[i]);
            singleWeight = oddRatioWeightSingleSNP(this.majorSNPReads[i], this.minorSNPReads[i]);
            singleVariant = oddRatioVarianceSingleSNP(this.majorSNPReads[i], this.minorSNPReads[i]);

            // 每个位点校正后的优势比(加权优势比)
            calibratedOddRatio = singleWeight * singleOR;
            calibratedSingleSNPOddRatio[i] = calibratedOddRatio;
            singleSNPWeight[i] = singleWeight;
            singleSNPVariant[i] = singleVariant;
        }

        double weightSum = getSum(singleSNPWeight);
        double calibratedORSum = getSum(calibratedSingleSNPOddRatio);
        double variantSum = getSum(singleSNPVariant);

        // 元分析，计算各个SNP位点整合优势比和对数整合优势比
        double metaOddRatio = calibratedORSum / weightSum;
        double logMetaOddRatio = Math.log(metaOddRatio);

        return new double[]{metaOddRatio, logMetaOddRatio, variantSum};
    }

    /**
     * 单个SNP位点的优势比
     * @param majorSNPReads 该位点 major haplotype上覆盖的reads数目
     * @param minorSNPReads 该位点 minor haplotype上覆盖的reads数目
     * @return 该位点处的优势比
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
        return (1 + (totalReads - 1) * this.rho) / (totalReads * majorSNPProbability * minorSNPProbability) + (1 + (totalReads - 1) * this.rho) * 4 / totalReads;
    }

    /**
     * 每个SNP位点的权重
     * @param majorSNPReads 该位点 major haplotype上覆盖的reads数目
     * @param minorSNPReads 该位点 minor haplotype上覆盖的reads数目
     * @return 位点在元分析中所占权重
     */
    private double oddRatioWeightSingleSNP(int majorSNPReads, int minorSNPReads) {
        double majorSNPProbability = (double)majorSNPReads / (double)(majorSNPReads + minorSNPReads);
        int totalReads = majorSNPReads + minorSNPReads;
        double calibrator = 1 + (totalReads - 1) * this.rho;
        return totalReads / (2 * calibrator) * 0.5 * (1 - majorSNPProbability);
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
