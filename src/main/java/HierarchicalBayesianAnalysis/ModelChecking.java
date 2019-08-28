package HierarchicalBayesianAnalysis;

import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.*;
import java.text.DecimalFormat;
import java.util.HashMap;

public class ModelChecking {
    private TauSampling ts;
    private LogOddRatioSampling lors;
    private double curTau, curTauPosteriorDensity, globalLORMean, globalLORSigma;
    private int samplingTime, burnIn;
    private double[] observeLogOddRatio, singleASELORMean, variances, minReplicatedData;
    private int[] majorAlleleReads, minorAlleleReads;
    private double curGlobalLOR, curGlobalLORPosteriorDensity;
    private double[] singleASELORMeanPosteriorDensity, singleASELOR,
            singleASELORPosteriorDensity, samplingGlobalLORs, samplingTaus;
    private int[] majorAlleleBackground, minorAlleleBackground;
    private double discrepancyMeasureRatio = 0;
    private DecimalFormat df = new DecimalFormat("0.0000");

    /**
     * Constructor
     */
    public ModelChecking(double minTau, double maxTau, int samplingTime, int burnIn,
                         int[] majorAlleleReads, int[] minorAlleleReads, int[] majorBackground, int[] minorBackground) {
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
        this.majorAlleleBackground = majorBackground;
        this.minorAlleleBackground = minorBackground;
        this.minReplicatedData = new double[samplingTime];
        this.samplingTaus = new double[samplingTime];
    }

    /**
     * 依据采样结果返回该peak的ASM显著性水平
     */
    public void testSignificant() {
        this.initializer();
        this.sampling();
        double pValue = this.discrepancyMeasureRatio / this.samplingTime;
        if (pValue < 0.05)
            System.out.println("extreme p value: " + pValue +", model need to revise");
        else
            System.out.println("model captures observe data according to posterior predictive result. p value = " + pValue);

        BufferedWriter tauBfw = null, globalLORBfw = null, minReplicatedLORBfw = null;
        String tauOutputFile = "C:\\Users\\hbs\\Desktop\\等位基因特异的m6A修饰位点分析平台\\model_checking\\tauRecord.txt";
        String globalLORFile = "C:\\Users\\hbs\\Desktop\\等位基因特异的m6A修饰位点分析平台\\model_checking\\globalLORRecord.txt";
        String replicateFile = "C:\\Users\\hbs\\Desktop\\等位基因特异的m6A修饰位点分析平台\\model_checking\\minLORRecord.txt";
        try {
            tauBfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(tauOutputFile))));
            globalLORBfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(globalLORFile))));
            minReplicatedLORBfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(replicateFile))));
            double minimumObserve = this.getMinLOR(this.observeLogOddRatio);
            minReplicatedLORBfw.write(df.format(minimumObserve));
            minReplicatedLORBfw.newLine();
            for (double lor: this.minReplicatedData) {
                minReplicatedLORBfw.write(df.format(lor));
                minReplicatedLORBfw.newLine();
            }
            for (double globalLOR: this.samplingGlobalLORs) {
                globalLORBfw.write(df.format(globalLOR));
                globalLORBfw.newLine();
            }
            for (double tau: this.samplingTaus) {
                tauBfw.write(this.df.format(tau));
                tauBfw.newLine();
            }
        } catch (IOException ie) {
            ie.printStackTrace();
        } finally {
            if (minReplicatedLORBfw != null) {
                try {
                    minReplicatedLORBfw.close();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
            if (tauBfw != null) {
                try {
                    tauBfw.close();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
            if (globalLORBfw != null) {
                try {
                    globalLORBfw.close();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
    }

    /**
     * 初始化参数tau等参数
     */
    private void initializer() {
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
    }

    /**
     * 第一轮采样前根据已有数据计算得到每个ASE位点的对数优势比及对数优势比的方差，相当于初始化值
     * @return {"LOR": [log odd ratios], "VAR": [variances]}
     */
    private HashMap<String, double[]> getInitLORAndVar() {
        OddRatioCalc orc = new OddRatioCalc(this.majorAlleleReads, this.minorAlleleReads, this.majorAlleleBackground, this.minorAlleleBackground);
        return orc.getLogOddRatio();
    }

    /**
     * 进行采样操作，得到 samplingTimes + burnIn个采样值
     */
    private void sampling() {
        int totalTimes = this.samplingTime + this.burnIn;
        for (int i=0; i<totalTimes; i++) {
            double prevTau = this.curTau;
            double prevTauPosteriorDensity = this.curTauPosteriorDensity;
            double[] samplingRes = this.ts.sampling(prevTau, prevTauPosteriorDensity, this.singleASELOR,
                    this.variances, this.globalLORMean, this.globalLORSigma);
            this.curTau = samplingRes[0];
            this.curTauPosteriorDensity = samplingRes[1];
            if (i > this.burnIn)
                this.samplingTaus[i-this.burnIn] = this.curTau;
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

            if (i > this.burnIn) {
                double observeDis = this.discrepancyMeasure(this.observeLogOddRatio, this.singleASELORMean);
                double replicationDis = this.discrepancyMeasure(this.singleASELOR, this.singleASELORMean);
                if (replicationDis >= observeDis) {
                    this.discrepancyMeasureRatio++;
                }
                this.minReplicatedData[i-this.burnIn] = this.getMinLOR(this.singleASELOR);
            }
        }
    }

    /**
     * 定义discrepancy measure用于posterior predictive p value求取过程中衡量replication data与observe data之间差异
     * @param logOddRatios 各个ASE位点上major allele reads的对数优势比
     * @param logOddRatioMeans 各个ASE位点上major allele reads的对数优势比的期望值
     * @return discrepancy measure
     */
    private double discrepancyMeasure(double[] logOddRatios, double[] logOddRatioMeans) {
        assert logOddRatios.length == logOddRatioMeans.length;
        double discrepancy = 0;
        for (int i=0; i<logOddRatios.length; i++) {
            discrepancy += (logOddRatios[i] - logOddRatioMeans[i]) / logOddRatioMeans[i];
        }

        return discrepancy;
    }

    /**
     * 获取生成数据中的最小值
     * @return 最小值
     */
    private double getMinLOR(double[] logOddRatioList) {
        double minimum = logOddRatioList[0];
        for (double lor: logOddRatioList) {
            if (lor < minimum)
                minimum = lor;
        }

        return minimum;
    }

    private String getString(double[] lorList) {
        String[] sb = new String[lorList.length];
        for (int i = 0; i < lorList.length; i++) {
            sb[i] = df.format(lorList[i]);
        }
        return String.join(", ", sb);
    }
}
