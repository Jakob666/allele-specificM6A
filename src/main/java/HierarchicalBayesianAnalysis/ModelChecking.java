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
    private double[] observeLogOddRatio, singleASELORMean, variances, replicatedData, minReplicatedData;
    private int[] majorAlleleReads, minorAlleleReads;
    private double discrepancyMeasureRatio = 0;
    private DecimalFormat df = new DecimalFormat("0.00");

    /**
     * Constructor
     */
    public ModelChecking(double minTau, double maxTau, int samplingTime, int burnIn,
                         int[] majorAlleleReads, int[] minorAlleleReads) {
        this.ts = new TauSampling(minTau, maxTau);
        this.lors = new LogOddRatioSampling();
        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
        this.majorAlleleReads = majorAlleleReads;
        this.minorAlleleReads = minorAlleleReads;
        this.singleASELORMean = new double[minorAlleleReads.length];
        this.replicatedData = new double[majorAlleleReads.length];
        this.minReplicatedData = new double[this.samplingTime];
    }

    /**
     * 依据采样结果返回该peak的ASE显著性水平
     */
    public void testSignificant() {
        this.initializer();
        this.sampling();
        double pValue = this.discrepancyMeasureRatio / this.samplingTime;
        if (pValue < 0.05)
            System.out.println("extreme p value: " + pValue +", model need to revise");
        else
            System.out.println("model captures observe data according to posterior predictive result. p value = " + pValue);

        BufferedWriter bfw = null;
        String outputFile = "C:\\Users\\hbs\\Desktop\\等位基因特异的m6A修饰位点分析平台\\model_checking\\minLORRecord.txt";
        try {
            bfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(outputFile))));
            double minimumObserve = this.getMinLOR(this.observeLogOddRatio);
            bfw.write(df.format(minimumObserve));
            bfw.newLine();
            for (double lor: this.minReplicatedData) {
                bfw.write(df.format(lor));
                bfw.newLine();
            }
            bfw.close();
        } catch (IOException ie) {
            ie.printStackTrace();
        } finally {
            if (bfw != null) {
                try {
                    bfw.close();
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
        // 计算观测到的ASE位点的对数优势比和方差
        HashMap<String, double[]> initLORAndVar = this.getInitLORAndVar();
        this.observeLogOddRatio = initLORAndVar.get("LOR");
        this.variances = initLORAndVar.get("VAR");
        System.out.println("observeLOR: [" + this.getString(this.observeLogOddRatio) + "]");
        System.out.println("observeVAR: [" + this.getString(this.variances) + "]");

        // 从tau的先验分布中初始化一个tau值
        this.curTau = this.ts.randomTau();
        // 依据初始的tau计算初始化的全局对数优势比的均值、方差和采样值
        double miu_denominator = 0, miu_numernator = 0;
        for (int i = 0; i < this.observeLogOddRatio.length; i++) {
            miu_denominator += 1.0 / (this.variances[i] + Math.pow(this.curTau, 2));
            miu_numernator += this.observeLogOddRatio[i] / (1.0 / (this.variances[i] + Math.pow(this.curTau, 2)));
        }
        this.globalLORMean = miu_numernator / miu_denominator;
        this.globalLORSigma = 1.0 / miu_denominator;
        NormalDistribution nd = new NormalDistribution(this.globalLORMean, this.globalLORSigma);
        double globalLOR = nd.sample();
        // 依据全局对数优势比的均值计算得到当前Tau的近似后验概率
        this.curTauPosteriorDensity = this.ts.posteriorTau(this.curTau, this.observeLogOddRatio, this.variances,
                this.globalLORMean, this.globalLORSigma);
        System.out.println("initial tau: " + this.curTau+"\t initial density: " + this.curTauPosteriorDensity);

        // 根据全局优势比得到各ASE位点优势比均值
        for (int i = 0; i < this.observeLogOddRatio.length; i++) {
            double mean = (1.0 / this.variances[i] * this.observeLogOddRatio[i] + 1.0 / Math.pow(this.curTau, 2) * globalLOR) / (1.0 / this.variances[i] + 1.0 / Math.pow(this.curTau, 2));
            double sigma = 1.0 / (1.0 / this.variances[i] + 1.0 / Math.pow(this.curTau, 2));
            nd = new NormalDistribution(mean, sigma);
            this.singleASELORMean[i] = nd.sample();
        }
        System.out.println("initial single ASE site LOR: [" + this.getString(this.singleASELORMean) + "]");
    }

    /**
     * 第一轮采样前根据已有数据计算得到每个ASE位点的对数优势比及对数优势比的方差，相当于初始化值
     * @return {"LOR": [log odd ratios], "VAR": [variances]}
     */
    private HashMap<String, double[]> getInitLORAndVar() {
        OddRatioCalc orc = new OddRatioCalc(this.majorAlleleReads, this.minorAlleleReads);
        return orc.getLogOddRatio();
    }

    /**
     * 进行采样操作，得到 samplingTimes + burnIn个采样值
     */
    private void sampling() {
        int totalTimes = this.samplingTime + this.burnIn;
        for (int i=0; i<totalTimes; i++) {
            // 首先对tau进行采样
            double prevTau = this.curTau;
            double prevTauPosteriorDensity = this.curTauPosteriorDensity;
            double[] samplingRes = this.ts.sampling(prevTau, prevTauPosteriorDensity, this.observeLogOddRatio,
                    this.variances, this.globalLORMean, this.globalLORSigma);
            this.curTau = samplingRes[0];
            this.curTauPosteriorDensity = samplingRes[1];

            // 对全局对数优势比进行采样
            double[] globalLORSummary = this.lors.globalLogOddRatioSampling(this.curTau, this.observeLogOddRatio, this.variances);
            this.globalLORMean = globalLORSummary[0];
            this.globalLORSigma = globalLORSummary[1];
            double globalLOR = globalLORSummary[2];
            System.out.println("step " + i + " current theta: " + df.format(globalLOR));

            // 对各个ASE位点的对数优势比均值进行采样
            this.singleASELORMean = this.lors.singleAseOddRatioSampling(this.curTau, globalLOR, this.observeLogOddRatio, this.variances);
            System.out.println("single ASE Site LOR: [" + this.getString(this.singleASELORMean) + "]");
            for (int j=0; j<this.singleASELORMean.length; j++) {
                this.replicatedData[j] = new NormalDistribution(this.singleASELORMean[j], this.variances[j]).sample();
            }
            System.out.println("replicate data: [" + this.getString(this.replicatedData) + "]");
            if (i > this.burnIn) {
                double observeDis = this.discrepancyMeasure(this.observeLogOddRatio, this.singleASELORMean);
                double replicationDis = this.discrepancyMeasure(this.replicatedData, this.singleASELORMean);
                if (replicationDis >= observeDis) {
                    this.discrepancyMeasureRatio++;
                }
                this.minReplicatedData[i-this.burnIn] = this.getMinLOR(this.replicatedData);
            }
            System.out.println("------------------------------------");
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
