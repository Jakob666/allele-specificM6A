package DifferentiationAnalysis;

import HierarchicalBayesianAnalysis.LogOddRatioSampler;
import HierarchicalBayesianAnalysis.TauSampler;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.Arrays;


public class ModelSelection {
    private int samplingTime, burnIn;
    private boolean sameLOR;
    private double gammaShape, gammaScale, curTau, curS1GlobalLOR, curS2GlobalLOR, curPosteriorDensity, maxPosteriorDensity=-1*Double.MAX_VALUE;
    private double bestS1LOR, bestS2LOR, bestTau, quantifiedS1LOR, quantifiedS2LOR;
    private double[] bestS1ExpectedLOR, bestS2ExpectedLOR;
    // samplingLOR, shape 1 × samplingTime
    private double[] s1GlobalLORList, s2GlobalLORList;
    // observed and expected LOR, shape 1 × allele number
    private double[] s1ObservedLOR, s2ObservedLOR, s1Variance, s2Variance, s1ExpectedLOR, s2ExpectedLOR;
    // parameter sampler
    private TauSampler ts;
    private LogOddRatioSampler lorSampler;

    public ModelSelection(double[] s1ObservedLOR, double[] s1Variance, double[] s2ObservedLOR, double[] s2Variance,
                          boolean sameLOR, int samplingTime, int burnIn, double df, double scaleParam) {
        assert s2Variance.length == s2ObservedLOR.length;
        assert s1Variance.length == s1ObservedLOR.length;

        this.s1ObservedLOR = s1ObservedLOR;
        this.s2ObservedLOR = s2ObservedLOR;
        this.s1Variance = s1Variance;
        this.s2Variance = s2Variance;
        this.s1ExpectedLOR = new double[s1ObservedLOR.length];
        this.s2ExpectedLOR = new double[s2ObservedLOR.length];
        this.sameLOR = sameLOR;
        this.samplingTime = samplingTime;
        this.burnIn = burnIn;

        this.gammaShape = df * 0.5;
        this.gammaScale = 2 / df / scaleParam;
        this.ts = new TauSampler(df, scaleParam);
        this.lorSampler = new LogOddRatioSampler();
    }

    public void initModel() {
        // randomly init value tau from priority distribution
        this.curTau = this.ts.randomInit();
        // init u with known tau, y and sigma
        double std, miuDenominator = 0, miuNumerator = 0, globalLORMean, globalLORSigma;
        NormalDistribution nd;
        if (this.sameLOR) {
            for (int i = 0; i < this.s1ObservedLOR.length; i++) {
                std = this.s1Variance[i] + Math.pow(this.curTau, 2);
                miuDenominator += 1.0 / std;
                miuNumerator += 1.0 / std * this.s1ObservedLOR[i];
            }
            globalLORMean = miuNumerator / miuDenominator; // u^
            globalLORSigma = 1.0 / miuDenominator;         // Vu
            nd = new NormalDistribution(globalLORMean, Math.sqrt(globalLORSigma));
            this.curS1GlobalLOR = nd.sample();
            this.curS2GlobalLOR = this.curS1GlobalLOR;
        } else {
            for (int i = 0; i < this.s1ObservedLOR.length; i++) {
                std = this.s1Variance[i] + Math.pow(this.curTau, 2);
                miuDenominator += 1.0 / std;
                miuNumerator += 1.0 / std * this.s1ObservedLOR[i];
            }
            globalLORMean = miuNumerator / miuDenominator; // u^
            globalLORSigma = 1.0 / miuDenominator;         // Vu
            nd = new NormalDistribution(globalLORMean, Math.sqrt(globalLORSigma));
            this.curS1GlobalLOR = nd.sample();
            for (int i = 0; i < this.s2ObservedLOR.length; i++) {
                std = this.s2Variance[i] + Math.pow(this.curTau, 2);
                miuDenominator += 1.0 / std;
                miuNumerator += 1.0 / std * this.s2ObservedLOR[i];
            }
            globalLORMean = miuNumerator / miuDenominator; // u^
            globalLORSigma = 1.0 / miuDenominator;         // Vu
            nd = new NormalDistribution(globalLORMean, Math.sqrt(globalLORSigma));
            this.curS2GlobalLOR = nd.sample();
        }

        // init vector theta with known u, tau, y and sigma
        double mean, var, lor;
        for (int i = 0; i < this.s1ObservedLOR.length; i++) {
            mean = (1.0 / this.s1Variance[i] * this.s1ObservedLOR[i] + 1.0 / Math.pow(this.curTau, 2) * this.curS1GlobalLOR) / (1.0 / this.s1Variance[i] + 1.0 / Math.pow(this.curTau, 2)); // theta^_j
            var = 1.0 / (1.0 / this.s1Variance[i] + 1.0 / Math.pow(this.curTau, 2));  // V_j
            nd = new NormalDistribution(mean, Math.sqrt(var));
            lor = nd.sample();
            this.s1ExpectedLOR[i] = lor;
        }

        for (int i = 0; i < this.s2ObservedLOR.length; i++) {
            mean = (1.0 / this.s2Variance[i] * this.s2ObservedLOR[i] + 1.0 / Math.pow(this.curTau, 2) * this.curS2GlobalLOR) / (1.0 / this.s2Variance[i] + 1.0 / Math.pow(this.curTau, 2)); // theta^_j
            var = 1.0 / (1.0 / this.s2Variance[i] + 1.0 / Math.pow(this.curTau, 2));  // V_j
            nd = new NormalDistribution(mean, Math.sqrt(var));
            lor = nd.sample();
            this.s2ExpectedLOR[i] = lor;
        }

        // logP(t | y) + logP(u | t, y) + logP(theta | u, t, y) + logP(y | theta, u, t)
        this.curPosteriorDensity = this.calcPosteriorProb(this.curTau, this.curS1GlobalLOR, this.curS2GlobalLOR, this.s1ExpectedLOR, this.s2ExpectedLOR);
        if (this.curPosteriorDensity - this.maxPosteriorDensity > 0.0000001) {
            this.maxPosteriorDensity = this.curPosteriorDensity;
            this.recordBestParams();
        }
    }

    public void sampling() {
        int totalTimes = this.samplingTime + this.burnIn;
        for (int i=0; i<totalTimes; i++) {
            this.tauSampling();
            this.uSampling(i);
            this.thetaSampling();
            this.curPosteriorDensity = this.calcPosteriorProb(this.curTau, this.curS1GlobalLOR, this.curS2GlobalLOR, this.s1ExpectedLOR, this.s2ExpectedLOR);
            if (this.curPosteriorDensity - this.maxPosteriorDensity > 0.0000001) {
                this.maxPosteriorDensity = this.curPosteriorDensity;
                this.recordBestParams();
            }
        }
    }

    private void tauSampling() {
        double newTau = this.ts.randomTau(this.curTau);
        double newPosteriorDensity = this.calcPosteriorProb(newTau, this.curS1GlobalLOR, this.curS2GlobalLOR, this.s1ExpectedLOR, this.s2ExpectedLOR);
        double[] samplingRes = this.ts.getSamplingRes(newTau, newPosteriorDensity, this.curTau, this.curPosteriorDensity);
        this.curTau = samplingRes[0];
        this.curPosteriorDensity = samplingRes[1];
    }

    private void uSampling(int time) {
        if (this.sameLOR) {
            double newGlobalLOR = this.lorSampler.randomU(this.curS1GlobalLOR);
            double newPosteriorDensity = this.calcPosteriorProb(this.curTau, newGlobalLOR, newGlobalLOR, this.s1ExpectedLOR, this.s2ExpectedLOR);
            double[] samplingRes = this.lorSampler.getSamplingRes(newGlobalLOR, newPosteriorDensity, this.curS1GlobalLOR, this.curPosteriorDensity);
            this.curS1GlobalLOR = samplingRes[0];
            this.curS2GlobalLOR = samplingRes[0];
            this.curPosteriorDensity = samplingRes[1];
        } else {
            double newGlobalLOR = this.lorSampler.randomU(this.curS1GlobalLOR);
            double newPosteriorDensity = this.calcPosteriorProb(this.curTau, newGlobalLOR, this.curS2GlobalLOR, this.s1ExpectedLOR, this.s2ExpectedLOR);
            double[] samplingRes = this.lorSampler.getSamplingRes(newGlobalLOR, newPosteriorDensity, this.curS1GlobalLOR, this.curPosteriorDensity);
            this.curS1GlobalLOR = samplingRes[0];
            this.curPosteriorDensity = samplingRes[1];

            newGlobalLOR = this.lorSampler.randomU(this.curS2GlobalLOR);
            newPosteriorDensity = this.calcPosteriorProb(this.curTau, this.curS1GlobalLOR, newGlobalLOR, this.s1ExpectedLOR, this.s2ExpectedLOR);
            samplingRes = this.lorSampler.getSamplingRes(newGlobalLOR, newPosteriorDensity, this.curS2GlobalLOR, this.curPosteriorDensity);
            this.curS2GlobalLOR = samplingRes[0];
            this.curPosteriorDensity = samplingRes[1];
        }
        if (time > this.burnIn) {
            if (this.s1GlobalLORList == null)
                this.s1GlobalLORList = new double[this.samplingTime];
            if (this.s2GlobalLORList == null)
                this.s2GlobalLORList = new double[this.samplingTime];
            this.s1GlobalLORList[time-this.burnIn] = this.curS1GlobalLOR;
            this.s2GlobalLORList[time-this.burnIn] = this.curS2GlobalLOR;
        }
    }

    private void thetaSampling() {
        double prevLORMean, newLORMean, prevPosteriorDensity, newPosteriorDensity;
        double[] samplingRes;
        for (int j=0; j<this.s1ObservedLOR.length; j++) {
            prevLORMean = this.s1ExpectedLOR[j];
            newLORMean = this.lorSampler.randomU(prevLORMean);
            prevPosteriorDensity = this.lorSampler.logPosteriorObserveLOR(this.s1ObservedLOR[j], prevLORMean, this.s1Variance[j]) +
                                   this.lorSampler.logPosteriorExpectedLOR(this.s1ObservedLOR[j], prevLORMean, this.s1Variance[j], this.curTau, this.curS1GlobalLOR);
            newPosteriorDensity = this.lorSampler.logPosteriorObserveLOR(this.s1ObservedLOR[j], newLORMean, this.s1Variance[j]) +
                                  this.lorSampler.logPosteriorExpectedLOR(this.s1ObservedLOR[j], newLORMean, this.s1Variance[j], this.curTau, this.curS1GlobalLOR);
            samplingRes = this.lorSampler.getSamplingRes(newLORMean, newPosteriorDensity, prevLORMean, prevPosteriorDensity);
            this.s1ExpectedLOR[j] = samplingRes[0];
        }
        for (int j=0; j<this.s2ObservedLOR.length; j++) {
            prevLORMean = this.s2ExpectedLOR[j];
            newLORMean = this.lorSampler.randomU(prevLORMean);
            prevPosteriorDensity = this.lorSampler.logPosteriorObserveLOR(this.s2ObservedLOR[j], prevLORMean, this.s2Variance[j]) +
                                   this.lorSampler.logPosteriorExpectedLOR(this.s2ObservedLOR[j], prevLORMean, this.s2Variance[j], this.curTau, this.curS2GlobalLOR);
            newPosteriorDensity = this.lorSampler.logPosteriorObserveLOR(this.s2ObservedLOR[j], newLORMean, this.s2Variance[j]) +
                                  this.lorSampler.logPosteriorExpectedLOR(this.s2ObservedLOR[j], newLORMean, this.s2Variance[j], this.curTau, this.curS2GlobalLOR);
            samplingRes = this.lorSampler.getSamplingRes(newLORMean, newPosteriorDensity, prevLORMean, prevPosteriorDensity);
            this.s2ExpectedLOR[j] = samplingRes[0];
        }
    }

    private double calcPosteriorProb(double tau, double s1GlobalLOR, double s2GlobalLOR, double[] s1ExpectedLOR, double[] s2ExpectedLOR) {
        return this.calcTauLogPosteriorProb(tau, this.s1Variance, this.s1ObservedLOR) +
               this.calcTauLogPosteriorProb(tau, this.s2Variance, this.s2ObservedLOR) +
               this.calcULogPosteriorProb(tau, s1GlobalLOR, this.s1Variance, this.s1ObservedLOR) +
               this.calcULogPosteriorProb(tau, s2GlobalLOR, this.s2Variance, this.s2ObservedLOR) +
               this.calcThetaLogPosteriorProb(tau, s1GlobalLOR, s1ExpectedLOR, this.s1Variance, this.s1ObservedLOR) +
               this.calcThetaLogPosteriorProb(tau, s2GlobalLOR, s2ExpectedLOR, this.s2Variance, this.s2ObservedLOR);
    }

    /**
     * calculate p(t | y)
     */
    private double calcTauLogPosteriorProb(double tau, double[] variances, double[] observedLOR) {
        // random sampling, renew u^, Vu, theta^ and V
        double[] param = this.calcUNormDistribParam(tau, variances, observedLOR);
        double globalLORVar = param[0];    // Vu
        double globalLORMean = param[1];   // u^

        return this.ts.logPosteriorTau(tau, observedLOR, variances, globalLORMean, globalLORVar);
    }

    /**
     * calculate p(u | t, y)
     */
    private double calcULogPosteriorProb(double tau, double u, double[] variances, double[] observedLOR) {
        double[] param = this.calcUNormDistribParam(tau, variances, observedLOR);
        double globalLORMean = param[0];   // u^
        double globalLORVar = param[1];    // Vu

        NormalDistribution nd = new NormalDistribution(globalLORMean, Math.sqrt(globalLORVar));

        return nd.logDensity(u);
    }

    /**
     * calculate sum of p(theta_j | u, t, y)
     */
    private double calcThetaLogPosteriorProb(double tau, double u, double[] theta, double[] variances, double[] observedLOR) {
        double[] logPosteriorProb = new double[theta.length];
        for (int i=0; i<variances.length; i++) {
            logPosteriorProb[i] = this.lorSampler.logPosteriorExpectedLOR(observedLOR[i], theta[i], variances[i], tau, u);
        }

        return Arrays.stream(logPosteriorProb).sum();
    }

    private double[] calcUNormDistribParam(double tau, double[] variances, double[] observedLOR) {
        double std, miuDenominator = 0, miuNumerator = 0;

        for (int i=0; i<variances.length; i++) {
            std = variances[i] + Math.pow(tau, 2);
            miuDenominator += 1 / std;
            miuNumerator += 1 / std * observedLOR[i];
        }
        double globalLORVar = 1 / miuDenominator;    // Vu
        double globalLORMean = miuNumerator / miuDenominator;    // u^

        return new double[] {globalLORMean, globalLORVar};
    }

    private void recordBestParams() {
        this.bestS1LOR = this.curS1GlobalLOR;
        this.bestS2LOR = this.curS2GlobalLOR;
        this.bestTau = this.curTau;
        this.bestS1ExpectedLOR = this.s1ExpectedLOR;
        this.bestS2ExpectedLOR = this.s2ExpectedLOR;
    }

    public double calcModelMarginProb() {
        ModelMarginProba mmp = new ModelMarginProba(this.s1ObservedLOR, this.s2ObservedLOR, this.s1Variance, this.s2Variance,
                                                    this.bestS1ExpectedLOR, this.bestS2ExpectedLOR, this.bestS1LOR, this.bestS2LOR,
                                                    this.maxPosteriorDensity, this.bestTau, this.gammaShape, this.gammaScale, this.sameLOR);

        return mmp.calcMarginProb();
    }

    public void quantify() {
        double[] sortedS1Value = Arrays.stream(this.s1GlobalLORList).sorted().toArray(),
                 sortedS2Value = Arrays.stream(this.s2GlobalLORList).sorted().toArray();
        int medianIdx = sortedS1Value.length / 2;
        if (sortedS1Value.length%2==0) {
            this.quantifiedS1LOR = (sortedS1Value[medianIdx] + sortedS1Value[medianIdx+1]) / 2;
            this.quantifiedS2LOR = (sortedS2Value[medianIdx] + sortedS2Value[medianIdx+1]) / 2;
        } else {
            this.quantifiedS1LOR = sortedS1Value[medianIdx];
            this.quantifiedS2LOR = sortedS2Value[medianIdx];
        }

        this.s1GlobalLORList = null;
        this.s2GlobalLORList = null;
    }

    public double getQuantifiedS1LOR() {
        return this.quantifiedS1LOR;
    }

    public double getQuantifiedS2LOR() {
        return this.quantifiedS2LOR;
    }
}
