package HierarchicalBayesianAnalysis;

import org.apache.commons.math3.distribution.NormalDistribution;
import java.util.HashMap;

public class HierarchicalBayesianModel {
    private TauSampling ts;
    private LogOddRatioSampling lors;
    private double curTau, curTauPosteriorDensity, miu;
    private int samplingTime, burnIn;
    private double[] logOddRatios, variances, samplingLORs;
    private int[] majorAlleleReads, minorAlleleReads;
    private double ratio = 0;

    /**
     * Constructor
     */
    public HierarchicalBayesianModel(double minTau, double maxTau, int samplingTime, int burnIn,
                                     int[] majorAlleleReads, int[] minorAlleleReads) {
        this.ts = new TauSampling(minTau, maxTau);
        this.lors = new LogOddRatioSampling();
        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
        this.majorAlleleReads = majorAlleleReads;
        this.minorAlleleReads = minorAlleleReads;
        this.samplingLORs = new double[samplingTime];
    }

    /**
     * 依据采样结果返回该peak的ASE显著性水平
     * @return p值
     */
    public double testSignificant() {
        this.initializer();
        this.sampling();
        for (int i=0; i<this.samplingLORs.length; i++) {
            if (i<this.burnIn)
                continue;
            double globalLOR = this.samplingLORs[i];
            if (globalLOR < 0)
                this.ratio++;
        }

        return 2 * this.ratio / this.samplingTime;
    }

    /**
     * 初始化参数tau等参数
     */
    private void initializer() {
        this.curTau = this.ts.randomTau();
        this.curTauPosteriorDensity = this.ts.priorTauDensity(this.curTau);
        HashMap<String, double[]> initLORAndVar = this.getInitLORAndVar();
        this.logOddRatios = initLORAndVar.get("LOR");
        this.variances = initLORAndVar.get("VAR");

        this.miu = 0;
        for (double variance: this.variances) {
            this.miu += 1.0 / (variance + Math.pow(this.curTau, 2));
        }
        this.miu = 1.0 / this.miu;
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
            double[] samplingRes = this.ts.sampling(prevTau, prevTauPosteriorDensity, this.logOddRatios,
                                                    this.variances, this.miu);
            this.curTau = samplingRes[0];
            this.curTauPosteriorDensity = samplingRes[1];

            // 对全局对数优势比进行采样
            double globalLORs = lors.globalLogOddRatioSampling(this.curTau, this.logOddRatios, this.variances);
            this.samplingLORs[i] = globalLORs;

            // 对各个ASE位点的对数优势比进行采样
            this.logOddRatios = lors.singleAseOddRatioSampling(this.curTau, globalLORs, this.logOddRatios, this.variances);
        }
    }
}
