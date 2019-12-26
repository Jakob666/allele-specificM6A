package HierarchicalBayesianAnalysis;

import java.util.HashMap;

/**
 * calculate the LOR and its variance of SNV sites locate in the range of a gene or m6A signal
 */
public class OddRatioCalc {
    private int[] majorSNPReads, minorSNPReads, majorAlleleBackground, minorAlleleBackground;

    /**
     * Constructor
     * @param majorSNPReads MeRIP-seq INPUT data major allele reads count
     * @param minorSNPReads MeRIP-seq INPUT data minor allele reads count
     * @param majorAlleleBackground WES data major allele reads count
     * @param minorAlleleBackground WES data minor allele reads count
     */
    public OddRatioCalc(int[] majorSNPReads, int[] minorSNPReads, int[] majorAlleleBackground, int[] minorAlleleBackground) {
        this.majorSNPReads = majorSNPReads;
        this.minorSNPReads = minorSNPReads;
        this.majorAlleleBackground = majorAlleleBackground;
        this.minorAlleleBackground = minorAlleleBackground;
    }

    /**
     * calculate major allele LOR and the variance
     * @return {"LOR": [log orr ratios], "VAR": [variances]}
     */
    public HashMap<String, double[]> getLogOddRatio() {
        double[] logOR = new double[this.majorSNPReads.length];
        double[] variance = new double[this.majorSNPReads.length];
        int majorReads, minorReads, majorBackground, minorBackground;
        double logOddratio, var;
        for (int i=0; i<this.majorSNPReads.length; i++) {
            majorReads = this.majorSNPReads[i];
            minorReads = this.minorSNPReads[i];
            if (this.majorAlleleBackground == null)
                majorBackground = (majorReads + minorReads) / 2;
            else
                majorBackground = (this.majorAlleleBackground[i] == 0)? (majorReads + minorReads) / 2: this.majorAlleleBackground[i];
            if (this.minorAlleleBackground == null)
                minorBackground = (majorReads + minorReads) / 2;
            else
                minorBackground = (this.minorAlleleBackground[i] == 0)? (majorReads + minorReads) / 2: this.minorAlleleBackground[i];
            logOddratio = this.calculateLogOddRatio(majorReads, minorReads, majorBackground, minorBackground);
            var = this.calculateVariance(majorReads, minorReads, majorBackground, minorBackground);
            logOR[i] = logOddratio;
            variance[i] = var;
        }
        HashMap<String, double[]> calcRes = new HashMap<>();
        calcRes.put("LOR", logOR);
        calcRes.put("VAR", variance);

        return calcRes;
    }

    /**
     * formula of LOR calculation
     *      yi = ln[y_major/(total-y_major)] - ln[0.5*total/(total-0.5*total)]
     *         = ln[y_major/(total-y_major)]
     * @param majorAlleleReads MeRIP-seq INPUT data major allele reads count
     * @param minorAlleleReads MeRIP-seq INPUT data minor allele reads count
     * @param majorBackground WES data major allele reads count
     * @param minorBackground WES data minor allele reads count
     * @return LOR of a SNV site
     */
    private double calculateLogOddRatio(int majorAlleleReads, int minorAlleleReads, int majorBackground, int minorBackground) {
        double oddRatio;
        if (minorAlleleReads == 0 | minorBackground == 0)
            oddRatio = ((majorAlleleReads + 0.5) / (minorAlleleReads + 0.5)) / ((majorBackground + 0.5) / (minorBackground + 0.5));
        else
            oddRatio = (double) majorAlleleReads / (double) (minorAlleleReads) / ((double) (majorBackground) / (double) (minorBackground));
        return Math.log(oddRatio);
    }

    /**
     * formula of LOR variance calculation
     *      s = 1/y_major + 1/(total-y_major) + 1/majorBackground + 1/minorBackground
     * @param majorAlleleReads MeRIP-seq INPUT data major allele reads count
     * @param minorAlleleReads MeRIP-seq INPUT data minor allele reads count
     * @param majorBackground WES data major allele reads count
     * @param minorBackground WES data minor allele reads count
     * @return LOR variance of a SNV site
     */
    private double calculateVariance(int majorAlleleReads, int minorAlleleReads, int majorBackground, int minorBackground) {
        if (minorAlleleReads == 0 | minorBackground == 0)
            return  1 / (majorAlleleReads + 0.5) + 1 / (minorAlleleReads + 0.5) + 1 / (majorBackground + 0.5) + 1 / (minorBackground + 0.5);
        else
            return 1/(double) majorAlleleReads + 1/(double)(minorAlleleReads) + 1/(double) (majorBackground) + 1/(double)(minorBackground);
    }
}
