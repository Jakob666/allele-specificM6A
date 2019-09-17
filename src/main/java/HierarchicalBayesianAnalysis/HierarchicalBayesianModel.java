package HierarchicalBayesianAnalysis;

import org.apache.commons.math3.distribution.NormalDistribution;

//import java.io.*;
import java.text.DecimalFormat;
import java.util.HashMap;

public class HierarchicalBayesianModel {
    private TauSampling ts;
    private LogOddRatioSampling lors;
    private double curTau, curGlobalLOR, curTauPosteriorDensity, curGlobalLORPosteriorDensity, globalLORMean, globalLORSigma;
    private int samplingTime, burnIn;
    private double[] observeLogOddRatio, variances, singleASELORMean, singleASELORMeanPosteriorDensity, singleASELOR,
                     singleASELORPosteriorDensity, samplingGlobalLORs;
    private int[] majorAlleleReads, minorAlleleReads, majorAlleleBackground, minorAlleleBackground;
    private DecimalFormat df = new DecimalFormat("0.00");

    /**
     * Constructor
     */
    public HierarchicalBayesianModel(double minTau, double maxTau, int samplingTime, int burnIn,
                                     int[] majorAlleleReads, int[] minorAlleleReads,
                                     int[] majorAlleleBackground, int[] minorAlleleBackground) {
        this.ts = new TauSampling(minTau, maxTau);
        this.lors = new LogOddRatioSampling();
        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
        this.majorAlleleReads = majorAlleleReads;
        this.minorAlleleReads = minorAlleleReads;
        this.samplingGlobalLORs = new double[samplingTime];
        this.singleASELORMean = new double[minorAlleleReads.length];
        this.singleASELORMeanPosteriorDensity = new double[minorAlleleReads.length];
        this.singleASELOR = new double[minorAlleleReads.length];
        this.singleASELORPosteriorDensity = new double[minorAlleleReads.length];
        this.majorAlleleBackground = majorAlleleBackground;
        this.minorAlleleBackground = minorAlleleBackground;
    }

    /**
     * 依据采样结果返回该peak的ASE显著性水平
     * @return p值
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
     * 初始化参数tau等参数
     */
    private void initializer() {
        // 计算观测到的ASE位点的对数优势比y和方差var
        HashMap<String, double[]> initLORAndVar = this.getInitLORAndVar();
        this.observeLogOddRatio = initLORAndVar.get("LOR");
        this.variances = initLORAndVar.get("VAR");
//        System.out.println("observeLOR: [" + this.getString(this.observeLogOddRatio) + "]");
//        System.out.println("observeVAR: [" + this.getString(this.variances) + "]");

        // 从tau的先验分布中初始化一个tau值
        this.curTau = this.ts.randomTau();

        // 已知tau、y和var，全局对数优势比的后验概率满足 P(globalLOR|tau,y)~N(Theta, V)。正态分布的两个变量可由tau、y和var求出
        // 初始化得到全局对数优势比的值
        double miu_denominator = 0, miu_numernator = 0;
        for (int i = 0; i < this.observeLogOddRatio.length; i++) {
            miu_denominator += 1.0 / (this.variances[i] + Math.pow(this.curTau, 2));
            miu_numernator += this.observeLogOddRatio[i] / (1.0 / (this.variances[i] + Math.pow(this.curTau, 2)));
        }
        this.globalLORMean = miu_numernator / miu_denominator;
        this.globalLORSigma = 1.0 / miu_denominator;
        NormalDistribution nd = new NormalDistribution(this.globalLORMean, this.globalLORSigma);
        this.curGlobalLOR = nd.sample();

        // 当前超参数Tau的近似后验概率
        this.curTauPosteriorDensity = this.ts.posteriorTau(this.curTau, this.observeLogOddRatio, this.variances,
                                                           this.globalLORMean, this.globalLORSigma);
        // 当前全局对数优势比的后验概率
        this.curGlobalLORPosteriorDensity = nd.density(this.curGlobalLOR);
//        System.out.println("initial tau: " + this.curTau+"\t initial density: " + this.curTauPosteriorDensity);

        // 已知y、var、Tau和globalLOR，各点对数优势比期望值P(LORMean_i|y,Tau,globalLOR)~N(Theta_i, Vi)。正态分布参数都可由已知量求出
        for (int i = 0; i < this.observeLogOddRatio.length; i++) {
            double mean = (1.0 / this.variances[i] * this.observeLogOddRatio[i] + 1.0 / Math.pow(this.curTau, 2) * this.curGlobalLOR) / (1.0 / this.variances[i] + 1.0 / Math.pow(this.curTau, 2));
            double sigma = 1.0 / (1.0 / this.variances[i] + 1.0 / Math.pow(this.curTau, 2));
            nd = new NormalDistribution(mean, sigma);
            double lor = nd.sample();
            // 初始化各位点的对数优势比期望值，并记录初始化优势比的后验概率
            this.singleASELORMean[i] = lor;
            this.singleASELORMeanPosteriorDensity[i] = nd.density(lor);
        }
//        System.out.println("initial single ASE site LOR mean: [" + this.getString(this.singleASELORMean) + "]");

        // 已知var、LORMean_i，各位点的对数优势比的后验概率 P(LOR_i|LORMean_i, var_i)~N(LORMean_i, var_i)
        for (int i = 0; i < this.singleASELORMean.length; i++) {
            double lorMean = this.singleASELORMean[i];
            double var = this.variances[i];
            nd = new NormalDistribution(lorMean, var);
            double lor = nd.sample();
            // 初始化一组LOR值并记录其后验概率
            this.singleASELOR[i] = lor;
            this.singleASELORPosteriorDensity[i] = nd.density(lor);
        }
//        System.out.println("initial single ASE site LOR: [" + this.getString(this.singleASELOR) + "]");
    }

    /**
     * 第一轮采样前根据已有数据计算得到每个ASE位点的对数优势比及对数优势比的方差，相当于初始化值
     * @return {"LOR": [log odd ratios], "VAR": [variances]}
     */
    private HashMap<String, double[]> getInitLORAndVar() {
        OddRatioCalc orc = new OddRatioCalc(this.majorAlleleReads, this.minorAlleleReads,
                                            this.majorAlleleBackground, this.minorAlleleBackground);
        return orc.getLogOddRatio();
    }

    /**
     * 进行采样操作，得到 samplingTimes + burnIn个采样值
     */
    private void sampling() {
        int totalTimes = this.samplingTime + this.burnIn;

        for (int i=0; i<totalTimes; i++) {
//            System.out.println("iteration " + i);
            // 首先对tau进行采样
            double prevTau = this.curTau;
            double prevTauPosteriorDensity = this.curTauPosteriorDensity;
            double[] samplingRes = this.ts.sampling(prevTau, prevTauPosteriorDensity, this.singleASELOR,
                    this.variances, this.globalLORMean, this.globalLORSigma);
            this.curTau = samplingRes[0];
            this.curTauPosteriorDensity = samplingRes[1];
//            System.out.println("new Tau value: " + this.curTau);

            // 对全局对数优势比进行采样
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
//            System.out.println("new global LOR value: " + this.curGlobalLOR);

            // 对各个ASE位点的对数优势比均值进行采样
            double[] prevSingleAseLORMean = this.singleASELORMean;
            double[] prevSingleAseLORMeanPosteriorDensity = this.singleASELORMeanPosteriorDensity;
            HashMap<String, double[]> singleAseMeanSampleRes = this.lors.singleAseOddRatioMeanSampling(this.curTau, this.curGlobalLOR, this.observeLogOddRatio,
                                                                                      this.variances, prevSingleAseLORMean, prevSingleAseLORMeanPosteriorDensity);
            this.singleASELORMean = singleAseMeanSampleRes.get("LOR");
            this.singleASELORMeanPosteriorDensity = singleAseMeanSampleRes.get("density");
//            System.out.println("new single site LOR mean value: [" + this.getString(this.singleASELORMean) + "]");

            // 对各个ASE位点的对数优势比进行采样
            double[] prevSingleAseLOR = this.singleASELOR;
            double[] prevSingleAseLORPosteriorDensity = this.singleASELORPosteriorDensity;
            HashMap<String, double[]> singleAseSampleRes = this.lors.singleAseOddRatioSampling(this.singleASELORMean, this.variances,
                                                                                               prevSingleAseLOR, prevSingleAseLORPosteriorDensity);
            this.singleASELOR = singleAseSampleRes.get("LOR");
            this.singleASELORPosteriorDensity = singleAseSampleRes.get("density");
//            System.out.println("new single site LOR value: [" + this.getString(this.singleASELOR) + "]");
//            System.out.println("==========================");
        }
    }

    private int getSum(int[] reads) {
        int total = 0;
        for (int r: reads) {
            total += r;
        }

        return total;
    }
}
