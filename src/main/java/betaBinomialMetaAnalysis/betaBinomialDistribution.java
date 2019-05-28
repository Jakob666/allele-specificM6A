package betaBinomialMetaAnalysis;

import org.apache.commons.math3.distribution.UniformRealDistribution;
import org.apache.commons.math3.special.Beta;
import org.apache.commons.math3.util.CombinatoricsUtils;

public class betaBinomialDistribution {

    /**
     * simulation of random sampling process from beta-binomial distribution
     * @param sampleSize sample number
     * @param n parameter of beta binomial distribution
     * @param theta parameter of beta binomial distribution
     * @param tao parameter of beta binomial distribution
     * @return random sampling result
     */
    public int[] sampling(int sampleSize, int n, double theta, double tao) {
        int[] randomSample = new int[sampleSize];

        UniformRealDistribution ufd = new UniformRealDistribution();
        double[] sampleProba = ufd.sample(sampleSize);
        for (int i = 0; i < sampleSize; i++) {
            randomSample[i] = sampleBinarySearch(n, theta, tao, sampleProba[i]);
        }
        return randomSample;
    }

    /**
     * implement beta-binomial distribution probability
     *      P(Y=y|n,theta,tao) = C(n, y) * Beta(y + theta/tao, n-y + (1-theta)/tao) / Beta(theta/tao, (1-theta)/tao)
     * @param n parameter of beta binomial distribution
     * @param y parameter of beta binomial distribution
     * @param theta parameter of beta binomial distribution
     * @param tao parameter of beta binomial distribution
     * @return probability value for a certain group of parameters
     */
    public double betaBinomialDistributionProbability(int n, int y, double theta, double tao) {
        double alpha1 = y+theta/tao;
        double beta1 = n-y+(1-theta)/tao;
        double alpha2 = theta/tao;
        double beta2 = (1-theta)/tao;

        double logBetaBinomial = logCombination(n, y) + logBetaFunction(alpha1, beta1) - logBetaFunction(alpha2, beta2);
        return Math.exp(logBetaBinomial);
    }

    /**
     * implement beta-binomial distribution cumulative probability
     * @param lowerBound start value of cdf
     * @param upperBound end value of cdf
     * @param n parameter of beta binomial distribution
     * @param theta parameter of beta binomial distribution
     * @param tao parameter of beta binomial distribution
     * @return cdf value cumulates from lower bound to upper bound
     */
    public double betaBinomialCdf(int lowerBound, int upperBound, int n, double theta, double tao) {
        assert (lowerBound < upperBound): "invalid input";
        double cdf = 0.0;
        while (lowerBound <= upperBound) {
            cdf = cdf + this.betaBinomialDistributionProbability(n, lowerBound, theta, tao);
            lowerBound++;
        }
        return cdf;
    }

    /**
     * implement beta distribution computational formula,
     *      beta(a, b) = gamma(a) * gamma(b) / gamma(a+b)
     * using logistic beta value to avoid arithmetic overflow
     * @return beta function log value for a particular alpha and beta parameter
     */
    private double logBetaFunction(double alpha, double beta){
        return Beta.logBeta(alpha, beta);
    }

    /**
     * calculate logistic combination value
     *      C(total, part) = total! / (part! * (total-part)!)
     *      log C(total, part) = log total! - log part! - log (total-part)!
     *                         = sum(log 1 to log total) - sum(log 1 to log part) - sum(log 1 to log total-part)
     * use logistic value aims to avoid arithmetic overflow
     * @param total total number
     * @param part sampling number
     * @return logistic combination value
     */
    private double logCombination(int total, int part) {
        return CombinatoricsUtils.binomialCoefficientLog(total, part);
    }

    /**
     * binary search of value n which satisfies betaBinomialCDF(n-1) < uniformProba <=  betaBinomialCDF(n)
     * the value n denotes as the sampling number which submit to beta binomial distribution
     * @param n parameter of beta binomial distribution
     * @param theta parameter of beta binomial distribution
     * @param tao parameter of beta binomial distribution
     * @param uniformProba parameter of beta binomial distribution
     * @return binary search value satisfies condition "betaBinomialCDF(n-1) < uniformProba <=  betaBinomialCDF(n)"
     */
    private int sampleBinarySearch(int n, double theta, double tao, double uniformProba) {
        int tmp;
        int lowerBound = 0;
        int upperBound = n;
        int mid = (lowerBound + n) / 2;
        double cdfProbaNMinus1, cdfProbaN;
        final double epsilon = 0.000001;
        while (mid > lowerBound) {
            cdfProbaNMinus1 = betaBinomialCdf(0, mid-1, n, theta, tao);
            cdfProbaN = betaBinomialCdf(0, mid, n, theta, tao);
            if (cdfProbaN-uniformProba >= epsilon && uniformProba-cdfProbaNMinus1 > epsilon) {
                return mid;
            } else if (cdfProbaN-uniformProba < epsilon) {
                tmp = mid;
                mid =  (upperBound + mid) / 2;
                lowerBound = tmp;
            } else if (uniformProba-cdfProbaNMinus1 < epsilon) {
                tmp = mid;
                mid = (lowerBound + mid) / 2;
                upperBound = tmp;
            }
        }
        return mid;
    }
}
