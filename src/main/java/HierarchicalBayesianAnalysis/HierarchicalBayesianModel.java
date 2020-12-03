package HierarchicalBayesianAnalysis;

import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.Arrays;
import java.util.HashMap;

public class HierarchicalBayesianModel {
    private TauSampler ts;
    private LogOddRatioSampler lors;
    private double curTau, curGlobalLOR, curPosteriorDensity;
    private int samplingTime, burnIn;
    private boolean allZero = true;
    // vector of y, sigma, p(y_j|theta_j), theta
    private double[] observeLogOddRatio, variances, singleASELORMean;
    // u sampling list
    private double[] samplingGlobalLORs = null;
    private int[] majorAlleleReads, minorAlleleReads, majorAlleleBackground, minorAlleleBackground;

    /**
     * Constructor
     */
    public HierarchicalBayesianModel(double df, double scaleParam, int samplingTime, int burnIn,
                                     int[] majorAlleleReads, int[] minorAlleleReads,
                                     int[] majorAlleleBackground, int[] minorAlleleBackground) {
        this.ts = new TauSampler(df, scaleParam);
        this.lors = new LogOddRatioSampler();
        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
        this.majorAlleleReads = majorAlleleReads;
        this.minorAlleleReads = minorAlleleReads;
        this.singleASELORMean = new double[minorAlleleReads.length];
        this.majorAlleleBackground = majorAlleleBackground;
        this.minorAlleleBackground = minorAlleleBackground;
    }

    /**
     * test for ASE or ASM significance. Return significant p value
     * @return p value
     */
    public double testSignificant() {
        this.checkLORVectors();
        if (this.allZero)
            return 1.0;
        this.initializer();
        this.sampling();
        double positiveLOR = 0, negativeLOR = 0;
        for (double globalLOR: this.samplingGlobalLORs) {
            if (globalLOR <= 0.0000001)
                negativeLOR++;
            else
                positiveLOR++;
        }

        return 2 * Math.min(positiveLOR, negativeLOR) / this.samplingTime;
    }

    private void checkLORVectors() {
        int idx = 0;
        while (idx < this.majorAlleleReads.length && this.allZero) {
            if (this.majorAlleleReads[idx] != this.minorAlleleReads[idx])
                this.allZero = false;
            idx++;
        }
    }

    /**
     * hierarchical model parameters initialization
     */
    private void initializer() {
        // init LOR(vector y) and variance(vector sigma)
        HashMap<String, double[]> initLORAndVar = this.getObserveData();
        this.observeLogOddRatio = initLORAndVar.get("LOR");
        this.variances = initLORAndVar.get("VAR");

        // randomly init value tau from priority distribution
        this.curTau = this.ts.randomInit();

        // init u with known tau, y and sigma
        double std, miuDenominator = 0, miuNumerator = 0;
        for (int i = 0; i < this.observeLogOddRatio.length; i++) {
            std = this.variances[i] + Math.pow(this.curTau, 2);
            miuDenominator += 1.0 / std;
            miuNumerator += 1.0 / std * this.observeLogOddRatio[i];
        }
        double globalLORMean = miuNumerator / miuDenominator; // u^
        double globalLORSigma = 1.0 / miuDenominator;         // Vu
        NormalDistribution nd = new NormalDistribution(globalLORMean, Math.sqrt(globalLORSigma));
        this.curGlobalLOR = nd.sample();

        // init vector theta with known u, tau, y and sigma
        for (int i = 0; i < this.observeLogOddRatio.length; i++) {
            double mean = (1.0 / this.variances[i] * this.observeLogOddRatio[i] + 1.0 / Math.pow(this.curTau, 2) * this.curGlobalLOR) / (1.0 / this.variances[i] + 1.0 / Math.pow(this.curTau, 2)); // theta^_j
            double var = 1.0 / (1.0 / this.variances[i] + 1.0 / Math.pow(this.curTau, 2));  // V_j
            nd = new NormalDistribution(mean, Math.sqrt(var));
            double lor = nd.sample();
            this.singleASELORMean[i] = lor;
        }

        // logP(t | y) + logP(u | t, y) + logP(theta | u, t, y) + logP(y | theta, u, t)
        this.curPosteriorDensity = this.calcPosteriorProb(this.curTau, this.curGlobalLOR, this.singleASELORMean);
    }

    /**
     * calculate LOR and variance via observed data
     * @return {"LOR": [log odd ratios], "VAR": [variances]}
     */
    private HashMap<String, double[]> getObserveData() {
        OddRatioCalc orc = new OddRatioCalc(this.majorAlleleReads, this.minorAlleleReads,
                                            this.majorAlleleBackground, this.minorAlleleBackground);
        return orc.getLogOddRatio();
    }

    /**
     * sampling the average LOR of gene or m6A peak，times = samplingTimes + burnIn
     */
    private void sampling() {
        int totalTimes = this.samplingTime + this.burnIn;

        double newTau, newGlobalLOR, prevPosteriorDensity, newPosteriorDensity;
        double[] samplingRes;
        for (int i=0; i<totalTimes; i++) {
            // sampling for tau
            newTau = this.ts.randomTau(this.curTau);
            newPosteriorDensity = this.calcPosteriorProb(newTau, this.curGlobalLOR, this.singleASELORMean);
            samplingRes = this.ts.getSamplingRes(newTau, newPosteriorDensity, this.curTau, this.curPosteriorDensity);
            this.curTau = samplingRes[0];
            this.curPosteriorDensity = samplingRes[1];

            // sampling for global LOR (u)
            newGlobalLOR = this.lors.randomU(this.curGlobalLOR);
            newPosteriorDensity = this.calcPosteriorProb(this.curTau, newGlobalLOR, this.singleASELORMean);
            samplingRes = this.lors.getSamplingRes(newGlobalLOR, newPosteriorDensity, this.curGlobalLOR, this.curPosteriorDensity);
            this.curGlobalLOR = samplingRes[0];
            this.curPosteriorDensity = samplingRes[1];
            if (i > this.burnIn) {
                if (this.samplingGlobalLORs == null)
                    samplingGlobalLORs = new double[this.samplingTime];
                this.samplingGlobalLORs[i-this.burnIn] = this.curGlobalLOR;
            }

            // sampling for the expected LOR for each SNV site (theta)
            double prevLORMean, newLORMean;
            for (int j=0; j<this.observeLogOddRatio.length; j++) {
                prevLORMean = this.singleASELORMean[j];
                newLORMean = this.lors.randomU(prevLORMean);
                prevPosteriorDensity = this.lors.logPosteriorObserveLOR(this.observeLogOddRatio[j], prevLORMean, this.variances[j]) +
                                       this.lors.logPosteriorExpectedLOR(this.observeLogOddRatio[j], prevLORMean, this.variances[j], this.curTau, this.curGlobalLOR);
                newPosteriorDensity = this.lors.logPosteriorObserveLOR(this.observeLogOddRatio[j], newLORMean, this.variances[j]) +
                                      this.lors.logPosteriorExpectedLOR(this.observeLogOddRatio[j], newLORMean, this.variances[j], this.curTau, this.curGlobalLOR);
                samplingRes = this.lors.getSamplingRes(newLORMean, newPosteriorDensity, prevLORMean, prevPosteriorDensity);
                this.singleASELORMean[j] = samplingRes[0];
            }
            this.curPosteriorDensity = this.calcPosteriorProb(this.curTau, this.curGlobalLOR, this.singleASELORMean);
        }
    }

    /**
     * calculate p(t | y)
     */
    private double calcTauLogPosteriorProb(double tau) {
        // random sampling, renew u^, Vu, theta^ and V
        double[] param = this.calcUNormDistribParam(tau);
        double globalLORVar = param[0];    // Vu
        double globalLORMean = param[1];   // u^

        return this.ts.logPosteriorTau(tau, this.observeLogOddRatio, this.variances, globalLORMean, globalLORVar);
    }

    /**
     * calculate p(u | t, y)
     */
    private double calcULogPosteriorProb(double tau, double u) {
        double[] param = this.calcUNormDistribParam(tau);
        double globalLORMean = param[0];   // u^
        double globalLORVar = param[1];    // Vu

        NormalDistribution nd = new NormalDistribution(globalLORMean, Math.sqrt(globalLORVar));

        return nd.logDensity(u);
    }

    /**
     * calculate sum of p(theta_j | u, t, y)
     */
    private double calcThetaLogPosteriorProb(double tau, double u, double[] theta) {
        double[] logPosteriorProb = new double[theta.length];
        for (int i=0; i<this.variances.length; i++) {
            logPosteriorProb[i] = this.lors.logPosteriorExpectedLOR(this.observeLogOddRatio[i], theta[i], this.variances[i], tau, u);
        }

        return Arrays.stream(logPosteriorProb).sum();
    }

    /**
     * calculate sum of p(y_j | theta_j, u, t) => p(y_j | theta_j)
     */
    @Deprecated
    private double calcYLogPosteriorProb(double[] theta) {
        double[] logPosteriorProb = new double[this.observeLogOddRatio.length];
        for (int i=0; i<this.observeLogOddRatio.length; i++) {
            logPosteriorProb[i] = this.lors.logPosteriorObserveLOR(this.observeLogOddRatio[i], theta[i], this.variances[i]);
        }

        return Arrays.stream(logPosteriorProb).sum();
    }

    private double calcPosteriorProb(double tau, double u, double[] theta) {
        return this.calcTauLogPosteriorProb(tau) + this.calcULogPosteriorProb(tau, u)
               + this.calcThetaLogPosteriorProb(tau, u, theta); // + this.calcYLogPosteriorProb(theta);
    }

    private double[] calcUNormDistribParam(double tau) {
        double std, miuDenominator = 0, miuNumerator = 0;

        for (int i=0; i<this.variances.length; i++) {
            std = this.variances[i] + Math.pow(tau, 2);
            miuDenominator += 1 / std;
            miuNumerator += 1 / std * this.observeLogOddRatio[i];
        }
        double globalLORVar = 1 / miuDenominator;    // Vu
        double globalLORMean = miuNumerator / miuDenominator;    // u^

        return new double[] {globalLORMean, globalLORVar};
    }

    /**
     * quantify the global log odd ratio with sampling median
     * @return gene logarithm odd ratio
     */
    public double quantifyGeneLOR() {
        double median;
        if (this.allZero)
            median = 0;
        else {
            double[] sortedSamplingValue = Arrays.stream(this.samplingGlobalLORs).sorted().toArray();
            int medianIdx = sortedSamplingValue.length / 2;
            if (sortedSamplingValue.length%2==0)
                median = (sortedSamplingValue[medianIdx] + sortedSamplingValue[medianIdx+1]) / 2;
            else
                median = sortedSamplingValue[medianIdx];

            this.samplingGlobalLORs = null;
        }

        return median;
    }

    public boolean isAllZero() {
        return this.allZero;
    }
}
