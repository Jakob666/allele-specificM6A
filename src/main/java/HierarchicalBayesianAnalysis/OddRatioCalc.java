package HierarchicalBayesianAnalysis;

import java.util.HashMap;

/**
 * 计算某个peak覆盖范围内，各ASE位点上major allele reads的对数优势比及优势比的方差
 */
public class OddRatioCalc {
    private int[] majorSNPReads, minorSNPReads;

    /**
     * Constructor
     * @param majorSNPReads 某个m6A信号下各 ASE位点上 major allele reads的数目
     * @param minorSNPReads 某个m6A信号下各 ASE位点上 minor allele reads的数目
     */
    public OddRatioCalc(int[] majorSNPReads, int[] minorSNPReads) {
        this.majorSNPReads = majorSNPReads;
        this.minorSNPReads = minorSNPReads;
    }

    /**
     * 计算得到peak覆盖范围内每个位点上major allele reads的对数优势比及对数优势比的方差
     * @return {"LOR": [log orr ratios], "VAR": [variances]}
     */
    public HashMap<String, double[]> getLogOddRatio() {
        double[] logOR = new double[this.majorSNPReads.length];
        double[] variance = new double[this.majorSNPReads.length];
        int majorReads, minorReads;
        double logOddratio, var;
        for (int i=0; i<this.majorSNPReads.length; i++) {
            majorReads = this.majorSNPReads[i];
            minorReads = this.minorSNPReads[i];
            logOddratio = this.calculateLogOddRatio(majorReads, minorReads);
            var = this.calculateVariance(majorReads, minorReads);
            logOR[i] = logOddratio;
            variance[i] = var;
        }
        HashMap<String, double[]> calcRes = new HashMap<>();
        calcRes.put("LOR", logOR);
        calcRes.put("VAR", variance);

        return calcRes;
    }

    /**
     * 计算某个ASE位点上 major allele reads的对数优势比，计算公式为
     *      yi = ln[y_major/(total-y_major)] - ln[0.5*total/(total-0.5*total)]
     *         = ln[y_major/(total-y_major)]
     * @param majorAlleleReads ASE位点上 major allele reads的数目
     * @param minorAlleleReads ASE位点上 minor allele reads的数目
     * @return 该位点上 major allele reads的对数优势比
     */
    private double calculateLogOddRatio(int majorAlleleReads, int minorAlleleReads) {
        double oddRatio = (double) majorAlleleReads / (double) (minorAlleleReads);
        return Math.log(oddRatio);
    }

    /**
     * 计算某个位点上 major allele reads对数优势比的方差，计算公式为
     *      s = 1/y_major + 1/(total-y_major) + 1/0.5*total + 1/0.5*total
     * @param majorAlleleReads ASE位点上 major allele reads的数目
     * @param minorAlleleReads ASE位点上 minor allele reads的数目
     * @return 该位点上 major allele reads的对数优势比的方差
     */
    private double calculateVariance(int majorAlleleReads, int minorAlleleReads) {
        return 1/(double) majorAlleleReads + 1/(double) (minorAlleleReads) + 4/(double) (majorAlleleReads+minorAlleleReads);
    }
}
