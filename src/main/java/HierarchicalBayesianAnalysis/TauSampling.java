package HierarchicalBayesianAnalysis;


public class TauSampling extends MHSampling {
    private InvChiSquareParams priorTau;

    /**
     * Constructor
     * @param lorStd the standard deviation of LOR
     * @param df Inv-Chi-square degree of freedom
     */
    public TauSampling(double lorStd, double df) {
        super();
        this.priorTau = new InvChiSquareParams(lorStd, df);
    }

    /**
     * sampling tau in new round
     * @param prevTau current tau
     * @param prevTauDensity the posterior density of current tau
     * @param logOddRatios LOR for each SNV sites
     * @param variances LOR variance for each SNV sites
     * @param miu expectation of globalLOR
     * @param sigma variance of globalLOR
     * @return sampling result and its posterior density
     */
    public double[] sampling(double prevTau, double prevTauDensity, double[] logOddRatios, double[] variances,
                             double miu, double sigma) {
        double curTau = this.randomTau();
        double curTauPosteriorDensity = this.posteriorTau(curTau, logOddRatios, variances, miu, sigma);

        return this.getSamplingRes(curTau, curTauPosteriorDensity, prevTau, prevTauDensity);
    }

    /**
     * sampling a new tau from its priority distribution
     * @return randomly sample tau
     */
    public double randomTau() {
        return this.priorTau.sample();
    }

    /**
     * the priority density of new round sampling tau
     * @param prevSamplingTau last round sampling tau
     * @return corresponding priority density
     */
    public double priorTauDensity(double prevSamplingTau) {
        return this.priorTau.density(prevSamplingTau);
    }

    /**
     * approximate posterior probability of sampling tau in new round
     * @param curSamplingTau new round sampling tau
     * @param logOddRatios LOR for each SNV sites
     * @param variances LOR variance for each SNV sites
     * @param miu expectation of globalLOR
     * @param sigma variance of globalLOR
     * @return approximate posterior density
     */
    public double posteriorTau(double curSamplingTau, double[] logOddRatios, double[] variances, double miu, double sigma) {
        double cumProd = 1;

        for (int i=0; i<variances.length; i++) {
            cumProd = cumProd * Math.pow(Math.pow(curSamplingTau, 2) + variances[i], -0.5) * Math.exp(-1 * Math.pow(logOddRatios[i]-miu, 2)/(2*(variances[i]+Math.pow(curSamplingTau, 2))));
        }

        double curTauPriorDensity = this.priorTauDensity(curSamplingTau);

        return curTauPriorDensity * Math.pow(sigma, 0.5) * cumProd;
    }
}
