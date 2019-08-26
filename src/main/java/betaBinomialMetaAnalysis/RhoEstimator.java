package betaBinomialMetaAnalysis;


public class RhoEstimator {

    private int[] majorSNPReadsCount, minorSNPReadsCount;
    private double learningRate, unImproveThreshold;
    private double r, delta, initialRho;

    /**
     * Constructor
     * @param majorSNPReadsCount reads count of each major SNP site
     * @param minorSNPReadsCount reads count of each minor SNP site
     * @param learningRate learning rate, recommend to set a small value such as 0.001. However, cautious about setting
     *                     an extremely tiny learning rate may train for a long time.
     * @param unImproveThreshold the threshold for early stop
     */
    public RhoEstimator(int[] majorSNPReadsCount, int[] minorSNPReadsCount, double initialRho, double learningRate, double unImproveThreshold){
        this.majorSNPReadsCount = majorSNPReadsCount;
        this.minorSNPReadsCount = minorSNPReadsCount;
        this.initialRho = initialRho;
        this.learningRate = learningRate;
        this.unImproveThreshold = unImproveThreshold;
        this.r = 0.0;
        this.delta = 0.0000001;
    }

    /**
     * estimate best rho parameter of beta binomial distribution with gradient ascend, update parameter rho with AdaGrad
     * optimizer
     * @return The value of rho which corresponding to the maximum value of Log likelihood function
     */
    public double gradientAscend() {
        int iter = 0;

        LogLikelihoodFunc llf = new LogLikelihoodFunc(this.majorSNPReadsCount, this.minorSNPReadsCount);
        RhoDerivation gradient = new RhoDerivation(this.majorSNPReadsCount, this.minorSNPReadsCount);

        double curTargetFuncVal = llf.logLikelihoodFunc(this.initialRho);
        double bestTargetFuncVal = curTargetFuncVal;
        double preTargetFuncVal = curTargetFuncVal;

        // these 2 values for early stop
        final double epsilon = 0.00001;
        int unImproveTime = 0;

        double curGradient;
        double curRho = this.initialRho;
        double bestRho = this.initialRho;

        while (true) {
//            if ((iter+1) % 100 == 0)
//                System.out.println(iter+1 +" times, curRho=" + curRho + ", bestRho=" + bestRho + ", LogLikelihoodFunc: " + curTargetFuncVal);
            // use first order derivation calculate gradient and update value of rho
            curGradient = gradient.firstOrderDerivation(curRho);

            // cumulative gradient square update
            this.r = this.r + Math.pow(curGradient, 2);

            curRho = curRho + this.learningRate / (this.delta + Math.sqrt(this.r)) * curGradient;

            // re-calculate the target function (log likelihood function)
            curTargetFuncVal = llf.logLikelihoodFunc(curRho);

            // renew best rho value
            if (curTargetFuncVal - bestTargetFuncVal > epsilon && curRho > 0) {
                bestRho = curRho;
                bestTargetFuncVal = curTargetFuncVal;
            }

            // early stop if there is no progress
            if (Math.abs(curTargetFuncVal - preTargetFuncVal) < this.unImproveThreshold) {
                unImproveTime++;
            } else {
                unImproveTime = 0;
            }

            if (unImproveTime >= 500) break;

            preTargetFuncVal = curTargetFuncVal;
            iter ++;
        }

        return bestRho;
    }
}
