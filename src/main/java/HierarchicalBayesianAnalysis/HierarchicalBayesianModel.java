package HierarchicalBayesianAnalysis;

import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.HashMap;

public class HierarchicalBayesianModel {
    private TauSampling ts;
    private LogOddRatioSampling lors;
    private double curTau, curGlobalLOR, curTauPosteriorDensity, curGlobalLORPosteriorDensity, globalLORMean, globalLORSigma;
    private int samplingTime, burnIn;
    private double[] observeLogOddRatio, variances, singleASELORMean, singleASELORMeanPosteriorDensity, singleASELOR,
                     samplingGlobalLORs;
    private int[] majorAlleleReads, minorAlleleReads, majorAlleleBackground, minorAlleleBackground;

    /**
     * Constructor
     */
    public HierarchicalBayesianModel(double lorStd, double df, int samplingTime, int burnIn,
                                     int[] majorAlleleReads, int[] minorAlleleReads,
                                     int[] majorAlleleBackground, int[] minorAlleleBackground) {
        this.ts = new TauSampling(lorStd, df);
        this.lors = new LogOddRatioSampling();
        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
        this.majorAlleleReads = majorAlleleReads;
        this.minorAlleleReads = minorAlleleReads;
        this.samplingGlobalLORs = new double[samplingTime];
        this.singleASELORMean = new double[minorAlleleReads.length];
        this.singleASELORMeanPosteriorDensity = new double[minorAlleleReads.length];
        this.singleASELOR = new double[minorAlleleReads.length];
        this.majorAlleleBackground = majorAlleleBackground;
        this.minorAlleleBackground = minorAlleleBackground;
    }

    /**
     * test for ASE or ASM significance. Return significant p value
     * @return p value
     */
    public double testSignificant() {
        this.initializer();
        this.sampling();
        double positiveLOR = 0, negativeLOR = 0;
        for (double globalLOR: this.samplingGlobalLORs) {
            if (globalLOR <= 0)
                negativeLOR++;
            else
                positiveLOR++;
        }

        return 2 * Math.min(positiveLOR, negativeLOR) / this.samplingTime;
    }

    /**
     * hierarchical model parameters initialization
     */
    private void initializer() {
        // calculate the LOR(y) and variance(var) for each bi-allele sites
        HashMap<String, double[]> initLORAndVar = this.getInitLORAndVar();
        this.observeLogOddRatio = initLORAndVar.get("LOR");
        this.variances = initLORAndVar.get("VAR");

        // randomly sample a value tau from priority distribution as initial value
        this.curTau = this.ts.randomInit();

        // when tau、y and var are known, the global LOR(average LOR) posterior distribution of gene or m6A signal
        //          P(globalLOR|tau,y)~N(Theta, V) can be calculate via y, var and tau
        // initial the global LOR from prior distribution
        double miu_denominator = 0, miu_numernator = 0;
        for (int i = 0; i < this.observeLogOddRatio.length; i++) {
            miu_denominator += 1.0 / (this.variances[i] + Math.pow(this.curTau, 2));
            miu_numernator += this.observeLogOddRatio[i] / (1.0 / (this.variances[i] + Math.pow(this.curTau, 2)));
        }
        this.globalLORMean = miu_numernator / miu_denominator;
        this.globalLORSigma = 1.0 / miu_denominator;
        NormalDistribution nd = new NormalDistribution(this.globalLORMean, this.globalLORSigma);
        this.curGlobalLOR = nd.sample();

        // the approximate posterior distribution of hyper-parameter tau
        this.curTauPosteriorDensity = this.ts.posteriorTau(this.curTau, this.observeLogOddRatio, this.variances,
                                                           this.globalLORMean, this.globalLORSigma);
        // the approximate posterior distribution of average LOR
        this.curGlobalLORPosteriorDensity = nd.density(this.curGlobalLOR);

        // when y, var, tau and globalLOR are known, the expectation LOR of each SNV sites(LORMean_i) satisfy normal distribution
        //               P(LORMean_i|y,Tau,globalLOR)~N(Theta_i, Vi)
        for (int i = 0; i < this.observeLogOddRatio.length; i++) {
            double mean = (1.0 / this.variances[i] * this.observeLogOddRatio[i] + 1.0 / Math.pow(this.curTau, 2) * this.curGlobalLOR) / (1.0 / this.variances[i] + 1.0 / Math.pow(this.curTau, 2));
            double sigma = 1.0 / (1.0 / this.variances[i] + 1.0 / Math.pow(this.curTau, 2));
            nd = new NormalDistribution(mean, sigma);
            double lor = nd.sample();
            // initial the expected LOR for each SNV sites and record the corresponding posterior probability
            this.singleASELORMean[i] = lor;
            this.singleASELORMeanPosteriorDensity[i] = nd.density(lor);
        }

        // when var, LORMean_i are known，LOR for each SNV sites satisfy normal distribution
        //          P(LOR_i|LORMean_i, var_i)~N(LORMean_i, var_i)
        for (int i = 0; i < this.singleASELORMean.length; i++) {
            double lorMean = this.singleASELORMean[i];
            double var = this.variances[i];
            nd = new NormalDistribution(lorMean, var);
            double lor = nd.sample();
            // initial LOR
            this.singleASELOR[i] = lor;
        }
    }

    /**
     * calculate LOR and variance via observed data
     * @return {"LOR": [log odd ratios], "VAR": [variances]}
     */
    private HashMap<String, double[]> getInitLORAndVar() {
        OddRatioCalc orc = new OddRatioCalc(this.majorAlleleReads, this.minorAlleleReads,
                                            this.majorAlleleBackground, this.minorAlleleBackground);
        return orc.getLogOddRatio();
    }

    /**
     * sampling the average LOR of gene or m6A peak，times = samplingTimes + burnIn
     */
    private void sampling() {
        int totalTimes = this.samplingTime + this.burnIn;

        for (int i=0; i<totalTimes; i++) {
            // sampling for tau
            double prevTau = this.curTau;
            double prevTauPosteriorDensity = this.curTauPosteriorDensity;
            double[] samplingRes = this.ts.sampling(prevTau, prevTauPosteriorDensity, this.singleASELOR,
                    this.variances, this.globalLORMean, this.globalLORSigma);
            this.curTau = samplingRes[0];
            this.curTauPosteriorDensity = samplingRes[1];

            // sampling for average LOR
            double prevGlobalLOR = this.curGlobalLOR;
            double prevGlobalLORPosteriorDensity = this.curGlobalLORPosteriorDensity;
            double[] globalLORSummary = this.lors.globalLogOddRatioSampling(this.curTau, this.observeLogOddRatio,
                                                                            this.variances, prevGlobalLOR, prevGlobalLORPosteriorDensity);
            this.globalLORMean = globalLORSummary[0];
            this.globalLORSigma = globalLORSummary[1];
            this.curGlobalLOR = globalLORSummary[2];
            this.curGlobalLORPosteriorDensity = globalLORSummary[3];
            if (i > this.burnIn)
                this.samplingGlobalLORs[i-this.burnIn] = this.curGlobalLOR;

            // sampling for the expected LOR for each SNV site
            double[] prevSingleAseLORMean = this.singleASELORMean;
            double[] prevSingleAseLORMeanPosteriorDensity = this.singleASELORMeanPosteriorDensity;
            HashMap<String, double[]> singleAseMeanSampleRes = this.lors.singleAseOddRatioMeanSampling(this.curTau, this.curGlobalLOR, this.observeLogOddRatio,
                                                                                      this.variances, prevSingleAseLORMean, prevSingleAseLORMeanPosteriorDensity);
            this.singleASELORMean = singleAseMeanSampleRes.get("LOR");
            this.singleASELORMeanPosteriorDensity = singleAseMeanSampleRes.get("density");
        }
    }
}
