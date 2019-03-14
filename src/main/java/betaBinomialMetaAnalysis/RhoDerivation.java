package betaBinomialMetaAnalysis;

public class RhoDerivation {

    private int totalSNP;
    private double majorProbability, minorProbability;
    private int[] majorSNPReadsCount, minorSNPReadsCount;

    /**
     * Constructor
     * @param majorSNPReadsCount reads count of each major SNP site
     * @param minorSNPReadsCount reads count of each minor SNP site
     */
    RhoDerivation(int[] majorSNPReadsCount, int[] minorSNPReadsCount) {
        this.majorSNPReadsCount = majorSNPReadsCount;
        this.minorSNPReadsCount = minorSNPReadsCount;

        // total SNP number on haplotype
        this.totalSNP = majorSNPReadsCount.length;

        // calculate the observe probability of major and minor SNP reads
        int totalMajorSNPReads = getSum(majorSNPReadsCount);
        int totalMinorSNPReads = getSum(minorSNPReadsCount);
        int totalReads = totalMajorSNPReads + totalMinorSNPReads;
        this.majorProbability = (double)totalMajorSNPReads / (double)totalReads;
        this.minorProbability = 1 - this.majorProbability;
    }

    /**
     * calculate the first order derivation value for a particular rho
     * @param rho the overdispersion value
     * @return gradient(first order derivation) value
     */
    public double firstOrderDerivation(double rho) {
        double firstOrderDerivationVal = (double)0;

        // first order derivation
        for (int k = 0; k < this.totalSNP; k++) {
            int siteiTotalReadsSNP = this.majorSNPReadsCount[k] + this.minorSNPReadsCount[k];

            for (int j = 0; j < this.majorSNPReadsCount[k]; j++) {
                firstOrderDerivationVal += (double)j / (this.majorProbability + j * rho);
            }

            for (int j = 0; j < this.minorSNPReadsCount[k]; j++) {
                firstOrderDerivationVal += (double)j / (this.minorProbability + j * rho);
            }

            for (int j = 0; j < siteiTotalReadsSNP; j++) {
                firstOrderDerivationVal -= (double)j / (1 + j * rho);
            }
        }
        return  firstOrderDerivationVal;
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
