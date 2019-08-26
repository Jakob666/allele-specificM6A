package HierarchicalBayesianAnalysis;

import java.util.HashMap;

/**
 * 计算某个peak覆盖范围内，各ASE位点上major allele reads的对数优势比及优势比的方差
 */
public class OddRatioCalc {
    private int[] majorSNPReads, minorSNPReads, majorAlleleBackground, minorAlleleBackground;

    /**
     * Constructor
     * @param majorSNPReads 某个m6A信号下各 ASE位点上 major allele reads的数目
     * @param minorSNPReads 某个m6A信号下各 ASE位点上 minor allele reads的数目
     * @param majorAlleleBackground WES数据得到的major allele背景域的reads数
     * @param minorAlleleBackground WES数据得到的minor allele背景域的reads数
     */
    public OddRatioCalc(int[] majorSNPReads, int[] minorSNPReads, int[] majorAlleleBackground, int[] minorAlleleBackground) {
        this.majorSNPReads = majorSNPReads;
        this.minorSNPReads = minorSNPReads;
        this.majorAlleleBackground = majorAlleleBackground;
        this.minorAlleleBackground = minorAlleleBackground;
    }

    /**
     * 计算得到peak覆盖范围内每个位点上major allele reads的对数优势比及对数优势比的方差
     * @return {"LOR": [log orr ratios], "VAR": [variances]}
     */
    public HashMap<String, double[]> getLogOddRatio() {
        double[] logOR = new double[this.majorSNPReads.length];
        double[] variance = new double[this.majorSNPReads.length];
        int majorReads, minorReads, majorBackground, minorBackground;
        double logOddratio, var;
        for (int i=0; i<this.majorSNPReads.length; i++) {
            majorReads = this.majorSNPReads[i];
            minorReads = this.minorSNPReads[i];
            if (this.majorAlleleBackground == null)
                majorBackground = (majorReads + minorReads) / 2;
            else
                majorBackground = this.majorAlleleBackground[i];
            if (this.minorAlleleBackground == null)
                minorBackground = (majorReads + minorReads) / 2;
            else
                minorBackground = this.minorAlleleBackground[i];
            logOddratio = this.calculateLogOddRatio(majorReads, minorReads, majorBackground, minorBackground);
            var = this.calculateVariance(majorReads, minorReads, majorBackground, minorBackground);
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
     * @param majorBackground WES得到的 major allele reads数目
     * @param minorBackground WES得到的 minor allele reads数目
     * @return 该位点上 major allele reads的对数优势比
     */
    private double calculateLogOddRatio(int majorAlleleReads, int minorAlleleReads, int majorBackground, int minorBackground) {
        double oddRatio = (double) majorAlleleReads / (double) (minorAlleleReads) / ((double) (majorBackground) / (double) minorBackground);
        return Math.log(oddRatio);
    }

    /**
     * 计算某个位点上 major allele reads对数优势比的方差，计算公式为
     *      s = 1/y_major + 1/(total-y_major) + 1/majorBackground + 1/minorBackground
     * @param majorAlleleReads ASE位点上 major allele reads的数目
     * @param minorAlleleReads ASE位点上 minor allele reads的数目
     * @param majorBackground WES得到的 major allele reads数目
     * @param minorBackground WES得到的 minor allele reads数目
     * @return 该位点上 major allele reads的对数优势比的方差
     */
    private double calculateVariance(int majorAlleleReads, int minorAlleleReads, int majorBackground, int minorBackground) {
        return 1/(double) majorAlleleReads + 1/(double) (minorAlleleReads) + 1/(double) (majorBackground) + 1/(double) (minorBackground);
    }
}
