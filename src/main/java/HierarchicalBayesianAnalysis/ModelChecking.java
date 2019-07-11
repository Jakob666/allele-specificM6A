package HierarchicalBayesianAnalysis;

import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class ModelChecking {
    private String discrepancyMethod;
    private double[] observedLogOddRatio, replicationLogOddRatio, variances, replicationDataParams;
    private double curTau, curTauPosteriorDensity, globalLORMean, globalLORSigma;
    private TauSampling ts;
    private int samplingTime, burnIn;
    private ArrayList<Double> minimumReplicationLOR;

    public ModelChecking(double minTau, double maxTau, int samplingTime, int burnIn,
                         int[] majorAlleleReads, int[] minorAlleleReads, String method) {
        this.ts = new TauSampling(minTau, maxTau);
        this.burnIn = burnIn;
        this.samplingTime = samplingTime;
        this.discrepancyMethod = method;
        this.replicationLogOddRatio = new double[majorAlleleReads.length];
        this.replicationDataParams = new double[majorAlleleReads.length];
        this.minimumReplicationLOR = new ArrayList<>();

        // 通过数据计算得到观察值对应的对数优势比、方差及discrepancy measure
        OddRatioCalc orc = new OddRatioCalc(majorAlleleReads, minorAlleleReads);
        HashMap<String, double[]> observeData = orc.getLogOddRatio();
        this.observedLogOddRatio = observeData.get("LOR");
        this.variances = observeData.get("VAR");
    }

    /**
     * 检验模型
     */
    public void checkModel() {
        double pValue = this.getPosteriorPredictivePValue();
        if (pValue < 0.05)
            System.out.println("extreme p value: " + pValue +", model need to revise");
        else
            System.out.println("model captures observe data according to posterior predictive result");
        minimumComparison();
    }

    /**
     * 返回posterior predictive p value
     * @return posterior predictive p value
     */
    public double getPosteriorPredictivePValue() {
        this.initial();

        int totalTime = this.samplingTime + this.burnIn;
        double observeDis, replicationDis, ratio = 0;
        for (int i=0; i<totalTime; i++) {
            // 首先对tau进行采样
            double prevTau = this.curTau;
            double prevTauPosteriorDensity = this.curTauPosteriorDensity;
            double[] samplingRes = this.ts.sampling(prevTau, prevTauPosteriorDensity, this.replicationLogOddRatio,
                                                    this.variances, this.globalLORMean, this.globalLORSigma);
            this.curTau = samplingRes[0];
            this.curTauPosteriorDensity = samplingRes[1];
            // 对全局对数优势比进行采样
            double globalLOR = this.globalLORSampling();
            // 对每个ASE位点的对数优势比期望进行采样，之后采样得到replication data
            this.singleAseSiteLORSampling(globalLOR);
            // 计算当前参数下replication data和observation data的discrepancy measure
            if (i >= this.burnIn) {
                double minimum = this.getMinimumLogOddRatio(this.replicationLogOddRatio);
                this.minimumReplicationLOR.add(minimum);
                observeDis = this.discrepancyMeasure(this.observedLogOddRatio, this.replicationDataParams, this.discrepancyMethod);
                replicationDis = this.discrepancyMeasure(this.replicationLogOddRatio, this.replicationDataParams, this.discrepancyMethod);
                if (replicationDis > observeDis)
                    ratio += 1;
            }
        }

        return ratio / this.samplingTime;
    }

    /**
     * 检验各组replication data得到的对数优势比的最小值是否能够包含观测数据对数优势比的最小值，并记录到文件以后可视化
     */
    public void minimumComparison() {
        // 观测数据中对数优势比最小值
        double minimumObserveLOR = this.getMinimumLogOddRatio(this.observedLogOddRatio);
        // 判断观测最小值是否在replication data最小值的分布中
        double minReplicateLORLowerBound = Collections.min(this.minimumReplicationLOR);
        double minReplicateLORUpperBound = Collections.max(this.minimumReplicationLOR);
        if (minimumObserveLOR >= minReplicateLORLowerBound && minimumObserveLOR <= minReplicateLORUpperBound)
            System.out.println("The hierarchical model clearly captured the observed variation");
        else
            System.out.println("The hierarchical model need to be revised");
        // 将数据写入到文件用于绘图
        BufferedWriter bfw = null;
        String outputFileName = "minimumLORDistribution.tsv";
        try {
            bfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(outputFileName))));
            bfw.write("ObserveMinimumLOR\t" + minimumObserveLOR);
            bfw.newLine();
            DecimalFormat df = new DecimalFormat("0.00");
            String[] minLOR = new String[this.minimumReplicationLOR.size()];
            for (int i=0; i<this.minimumReplicationLOR.size(); i++) {
                minLOR[i] = df.format(this.minimumReplicationLOR.get(i));
            }
            bfw.write("replicateMinimumLOR\t" + String.join(",", minLOR));
            bfw.newLine();
            bfw.close();
        } catch (Exception e) {
            System.out.println("Failed to write in file " + outputFileName);
        } finally {
            if (bfw != null) {
                try {
                    bfw.close();
                } catch (IOException ie) {
                    ie.printStackTrace();
                }
            }
        }
    }

    /**
     * 初始化计算过程中需要的值
     */
    private void initial() {
        this.replicationLogOddRatio = this.observedLogOddRatio;
        // 依据初始的tau
        this.curTau = this.ts.randomTau();
        this.curTauPosteriorDensity = this.ts.priorTauDensity(this.curTau);
        // 计算初始化的全局对数优势比的值
        double globalLOR = this.globalLORSampling();
        // 根据全局优势比初始化各个位点优势比
        this.singleAseSiteLORSampling(globalLOR);
    }

    /**
     * 定义discrepancy measure用于posterior predictive p value求取过程中衡量replication data与observe data之间差异
     * @param logOddRatios 各个ASE位点上major allele reads的对数优势比
     * @param logOddRatioMeans 各个ASE位点上major allele reads的对数优势比的期望值
     * @param methodName discrepancy measure的方法名
     * @return discrepancy measure
     */
    private double discrepancyMeasure(double[] logOddRatios, double[] logOddRatioMeans, String methodName) {
        assert logOddRatios.length == logOddRatioMeans.length;
        double discrepancy = 0;
        String method = methodName.toLowerCase();
        if (method.equals("gelman")) {
            for (int i=0; i<logOddRatios.length; i++) {
                discrepancy += (logOddRatios[i] - logOddRatioMeans[i]) / logOddRatioMeans[i];
            }
        }
        else if (method.equals("likelihood")) {
            for (int i=0; i< logOddRatios.length; i++) {
                discrepancy += logOddRatios[i] * Math.log(logOddRatios[i] / logOddRatioMeans[i]);
            }
            discrepancy = discrepancy * 2;
        }
        else if (method.equals("freeman")) {
            for (int i=0; i<logOddRatios.length; i++) {
                discrepancy += Math.pow(Math.sqrt(logOddRatios[i])-Math.sqrt(logOddRatioMeans[i]), 2);
            }
        }

        return discrepancy;
    }

    /**
     * 对全局对数优势比进行采样
     * @return 全局对数优势比
     */
    private double globalLORSampling() {
        double miu_denominator = 0, miu_numernator = 0;
        for (int i = 0; i < this.replicationLogOddRatio.length; i++) {
            miu_denominator += 1.0 / (this.variances[i] + Math.pow(this.curTau, 2));
            miu_numernator += this.replicationLogOddRatio[i] / (1.0 / (this.variances[i] + Math.pow(this.curTau, 2)));
        }
        this.globalLORMean = miu_numernator / miu_denominator;
        this.globalLORSigma = 1.0 / miu_denominator;
        NormalDistribution nd = new NormalDistribution(this.globalLORMean, this.globalLORSigma);
        return nd.sample();
    }

    /**
     * 生成各个ASE位点的对数优势比
     * @param globalLOR 当前全局优势比
     */
    private void singleAseSiteLORSampling(double globalLOR) {
        NormalDistribution nd;
        for (int i = 0; i < this.replicationLogOddRatio.length; i++) {
            double mean = (1.0 / this.variances[i] * this.replicationLogOddRatio[i] + 1.0 / Math.pow(this.curTau, 2) * globalLOR) / (1.0 / this.variances[i] + 1.0 / Math.pow(this.curTau, 2));
            double sigma = 1.0 / (1.0 / this.variances[i] + 1.0 / Math.pow(this.curTau, 2));
            nd = new NormalDistribution(mean, sigma);
            double lorMean = nd.sample();
            // 各ASE位点对数优势比及其期望值
            this.replicationLogOddRatio[i] = new NormalDistribution(lorMean, this.variances[i]).sample();
            this.replicationDataParams[i] = mean;
        }
    }

    /**
     * 返回replication data中的最小值
     * @return replication data中的最小值
     */
    private double getMinimumLogOddRatio(double[] logOddRatioList) {
        double minimum = logOddRatioList[0];
        for (double lor: logOddRatioList) {
            if (lor < minimum)
                minimum = lor;
        }

        return minimum;
    }
}
