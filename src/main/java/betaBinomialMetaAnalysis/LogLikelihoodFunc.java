package betaBinomialMetaAnalysis;


public class LogLikelihoodFunc {

    private int totalSNP;
    private double majorProbability, minorProbability;
    private int[] majorSNPReadsCount, minorSNPReadsCount;

    /**
     * Constructor
     * @param majorSNPReadsCount reads count of each major SNP site
     * @param minorSNPReadsCount reads count of each minor SNP site
     */
    public LogLikelihoodFunc(int[] majorSNPReadsCount, int[] minorSNPReadsCount) {
        this.majorSNPReadsCount = majorSNPReadsCount;
        this.minorSNPReadsCount = minorSNPReadsCount;

        this.totalSNP = majorSNPReadsCount.length;

        int totalMajorSNPReads = getSum(majorSNPReadsCount);
        int totalMinotSNPReads = getSum(minorSNPReadsCount);
        int totalSNPReads = totalMajorSNPReads + totalMinotSNPReads;
        this.majorProbability = (double)totalMajorSNPReads / (double)totalSNPReads;
        this.minorProbability = 1 - majorProbability;
    }

    /**
     * log likelihood function of beta binomial distribution
     * @param rho overdispersion value
     * @return the likelihood function result
     */
    public double logLikelihoodFunc(double rho) {

        double logLikelihoodValue = 0.0;

        // log-likelihood function
        for (int i = 0; i < this.totalSNP; i++) {
            int siteiTotalReadsSNP = this.majorSNPReadsCount[i] + this.minorSNPReadsCount[i];

            // a constant during the calculation procedure, skip it can greatly improve the optimization speed
//             for (int j = 0; j < this.totalSNP; j++) {
//                 logLikelihoodValue += CombinatoricsUtils.binomialCoefficientLog(siteiTotalReadsSNP, this.majorSNPReadsCount[i]);
//             }

            for (int k = 0; k < this.majorSNPReadsCount[i]; k++) {
                logLikelihoodValue = logLikelihoodValue + Math.log(rho * k + this.majorProbability);
            }

            for (int l = 0; l < this.minorSNPReadsCount[i]; l++) {
                logLikelihoodValue = logLikelihoodValue + Math.log(this.minorProbability + l * rho);
            }

            for (int m = 0; m < siteiTotalReadsSNP; m++) {
                logLikelihoodValue = logLikelihoodValue - Math.log(1 + m * rho);
            }
        }
        return logLikelihoodValue;
    }


    /**
     * calculate reads counts
     * @param SNPReads int array records reads number of SNP sites
     * @return sum of reads on SNP sites
     */
    private int getSum(int[] SNPReads) {
        int sum = 0;
        for (int i = 0; i < SNPReads.length; i++) {
            sum += SNPReads[i];
        }
        return sum;
    }
}
