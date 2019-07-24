package HierarchicalBayesianAnalysis;

import GTFComponent.VcfSnpMatchGene;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

/**
 * 通过SNP calling得到的VCF文件对ASE gene进行检验
 */
public class AseGeneDetection {
    // geneAlleleReads = {"geneId->geneName": {"major":[count1, count2,...], "minor": [count1, count2,...]}, ...}
    private String aseGeneFile;
    private HashMap<String, HashMap<String, ArrayList<Integer>>> geneAlleleReads;
    private HashMap<String, int[]> geneMajorStrand;
    private int samplingTime, burnIn;
    private HashMap<String, Double> geneAsePValue = new HashMap<>(), geneAseQValue = new HashMap<>();
    private DecimalFormat df = new DecimalFormat("0.0000");

    /**
     * Constructor
     * @param gtfFile GTF文件
     * @param vcfFile SNP calling得到的VCF文件
     * @param aseGeneFile 结果输出文件
     * @param samplingTime 采样次数
     * @param burnIn burn in次数
     */
    public AseGeneDetection(String gtfFile, String vcfFile, String aseGeneFile, int samplingTime, int burnIn) {
        VcfSnpMatchGene vsmg = new VcfSnpMatchGene(vcfFile, gtfFile);
        vsmg.parseVcfFile();
        this.geneAlleleReads = vsmg.getGeneAlleleReads();
        this.geneMajorStrand = vsmg.getGeneMajorAlleleStrand();
        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
        this.aseGeneFile = aseGeneFile;
    }

    public void getTestResult() {
        this.aseGeneTest();
        this.bhRecalibrationOfEachPeak();
        this.outputResult();
    }

    /**
     * 使用层次模型对基因的ASE进行检验得到显著性p值
     */
    private void aseGeneTest() {
        ArrayList<Integer> majorAlleleCount, minorAlleleCount;
        int[] majorCount, minorCount;
        HierarchicalBayesianModel hb;
        double p;
        for (String label: this.geneAlleleReads.keySet()) {
            majorAlleleCount = this.geneAlleleReads.get(label).get("major");
            minorAlleleCount = this.geneAlleleReads.get(label).get("minor");

            assert majorAlleleCount.size() == minorAlleleCount.size();
            majorCount = new int[majorAlleleCount.size()];
            minorCount = new int[minorAlleleCount.size()];
            for (int i=0; i<majorAlleleCount.size(); i++) {
                majorCount[i] = majorAlleleCount.get(i);
                minorCount[i] = minorAlleleCount.get(i);
            }
            hb = new HierarchicalBayesianModel(0, 1, this.samplingTime, this.burnIn, majorCount, minorCount);
            p = hb.testSignificant();
            this.geneAsePValue.put(label, p);
        }
    }

    /**
     * 对基因ASE的p值进行BH校正, 得到 Q value.
     */
    private void bhRecalibrationOfEachPeak() {
        List<Map.Entry<String, Double>> sortedByValue = new ArrayList<Map.Entry<String, Double>>(this.geneAsePValue.entrySet());
        // 依据ASE的p值从小到大排序
        Collections.sort(sortedByValue, new Comparator<Map.Entry<String, Double>>() {
            public int compare(Map.Entry<String, Double> o1,
                               Map.Entry<String, Double> o2) {
                return (o2.getValue()).compareTo(o1.getValue());
            }
        });

        int rankage = 0;
        int totalPeak = sortedByValue.size();
        double pValue, qValue, prevQVal = Collections.max(this.geneAsePValue.values());
        String label;
        for (Map.Entry<String, Double> record: sortedByValue) {
            rankage += 1;
            label = record.getKey();
            pValue = record.getValue();
            qValue = pValue * totalPeak / rankage;
            if (qValue > prevQVal) {
                qValue = prevQVal;
            } else {
                prevQVal = Math.min(1.0, qValue);
            }
            this.geneAseQValue.put(label, qValue);
        }
    }

    /**
     * 将检验结果写入文件
     */
    private void outputResult() {
        List<Map.Entry<String, Double>> sortedByValue = new ArrayList<Map.Entry<String, Double>>(this.geneAseQValue.entrySet());
        // 依据ASE的q值从小到大排序
        Collections.sort(sortedByValue, new Comparator<Map.Entry<String, Double>>() {
            public int compare(Map.Entry<String, Double> o1,
                               Map.Entry<String, Double> o2) {
                return (o1.getValue()).compareTo(o2.getValue());
            }
        });
        BufferedWriter bfw = null;
        try {
            bfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(this.aseGeneFile))));
            String line, label, geneId, geneName, majorAlleleStrand;
            String[] info;
            int[] alleleStrand;
            double qVal;
            ArrayList<Integer> majorAlleleCount, minorAlleleCount;
            bfw.write("#geneId\tgeneName\tq-value\tmajorAlleleReads\tminorAlleleReads\tmajorAlleleStrand\n");
            for (Map.Entry<String, Double> record: sortedByValue) {
                label = record.getKey();
                info = label.split("->");
                geneId = info[0];
                geneName = info[1];
                qVal = record.getValue();
                majorAlleleCount = this.geneAlleleReads.get(label).get("major");
                minorAlleleCount = this.geneAlleleReads.get(label).get("minor");
                alleleStrand = this.geneMajorStrand.get(geneId);
                if (alleleStrand[0] >= alleleStrand[1])
                    majorAlleleStrand = "+";
                else
                    majorAlleleStrand = "-";

                line = String.join("\t", new String[]{geneId, geneName, this.df.format(qVal),
                                    majorAlleleCount.toString(), minorAlleleCount.toString(), majorAlleleStrand});
                bfw.write(line);
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
}
