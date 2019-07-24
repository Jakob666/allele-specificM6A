package HierarchicalBayesianAnalysis;


import AseM6aPeakDetector.HeterozygoteReadsCount;
import org.apache.log4j.Logger;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

public class AsmPeakDetection {
    private String peakBedFile, vcfFile, asmPeakFile, peakCoveredSnpFile;
    private int samplingTime, burnIn;
    private Logger log;
    private HashMap<String, HashMap<String, HashMap<String, LinkedList<Integer>>>> peakSnpReadsCount;
    private HashMap<String, Double> asmPValue = new HashMap<>(), asmQValue = new HashMap<>();
    private HashMap<String, String> peakCoveredGene = new HashMap<>(), peakMajorAlleleStrand;
    private DecimalFormat df = new DecimalFormat("0.0000");

    /**
     * Constructor
     * @param peakBedFile peak calling得到的BED格式文件
     * @param vcfFile IP样本SNP Calling得到的VCF格式文件
     * @param peakCoveredSnpFile 记录peak信号覆盖的SNP位点的文件
     * @param asmPeakFile ASM peak检验结果输出文件
     * @param samplingTime 采样次数
     * @param burnIn burn in次数
     */
    public AsmPeakDetection(String peakBedFile, String vcfFile, String peakCoveredSnpFile, String asmPeakFile,
                            int samplingTime, int burnIn) {
        this.peakBedFile = peakBedFile;
        this.vcfFile = vcfFile;
        this.peakCoveredSnpFile = peakCoveredSnpFile;
        this.asmPeakFile = asmPeakFile;
        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
        this.log = initLog(new File(asmPeakFile).getParent());
    }

    public void getTestResult() {
        this.getPeakCoveredGene();;
        this.getPeakCoveredSnpResult();
        this.getPeakSNPReadsCount();
        this.asmPeakTest();
        this.bhRecalibrationOfEachPeak();
        this.outputResult();
    }

    /**
     * peak覆盖的基因 peakCoveredGene = {chr:peakStart:peakEnd -> geneId, ...}
     */
    private void getPeakCoveredGene() {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.peakBedFile))));
            String line = "", chrNum, peakStart, peakEnd, geneId, label;
            String[] info;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    chrNum = info[0];
                    peakStart = info[1];
                    peakEnd = info[2];
                    geneId = info[3];

                    label = String.join(":", new String[]{chrNum, peakStart, peakEnd});
                    this.peakCoveredGene.put(label, geneId);
                }
            }
            bfr.close();
        } catch (IOException ie) {
            ie.printStackTrace();
        } finally {
            if (bfr != null) {
                try {
                    bfr.close();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
    }

    /**
     * 获取被peak覆盖的SNP位点，并将记录写入文件
     */
    private void getPeakCoveredSnpResult() {
        PeakCoveredSNPRecord pcsr = new PeakCoveredSNPRecord(this.vcfFile, this.peakBedFile, this.peakCoveredSnpFile);
        pcsr.getPeakCoveredSNP();
    }

    /**
     * 获取每个m6A信号范围内ASE位点上的 major allele和 minor allele的reads count
     * @return [chr1: [peak1: [major: [count1, count2,...], minor: [count1, count2,...]]], chr2:....]
     */
    private HashMap<String, HashMap<String, HashMap<String, LinkedList<Integer>>>> getPeakSNPReadsCount() {
        HeterozygoteReadsCount hrc = new HeterozygoteReadsCount(this.peakCoveredSnpFile, this.log);
        this.peakMajorAlleleStrand = hrc.getPeakMajorAlleleStrand();
        return hrc.getMajorMinorHaplotype();
    }

    /**
     * 对每个覆盖SNP的peak进行检验
     */
    private void asmPeakTest() {
        // [chr1: [peak1: [major: [count1, count2,...], minor: [count1, count2,...]]], chr2:....]
        this.peakSnpReadsCount = this.getPeakSNPReadsCount();
        this.hierarchicalModelTest();
    }

    /**
     * 对每个覆盖了SNP的peak进行检验，得到ASM显著性p值
     */
    private void hierarchicalModelTest() {
        double pVal;
        Set<Map.Entry<String, HashMap<String, HashMap<String, LinkedList<Integer>>>>> chrPeaks = this.peakSnpReadsCount.entrySet();
        for (Map.Entry<String, HashMap<String, HashMap<String, LinkedList<Integer>>>> chrPeak: chrPeaks) {
            String chrNum = chrPeak.getKey();
            // 某条染色体上所有的peak
            HashMap<String, HashMap<String, LinkedList<Integer>>> peaksOnChr = chrPeak.getValue();

            Set<Map.Entry<String, HashMap<String, LinkedList<Integer>>>> peakRanges = peaksOnChr.entrySet();
            // 对每个peak求取p值
            for (Map.Entry<String, HashMap<String, LinkedList<Integer>>> peakRange: peakRanges) {
                // key记录peak的范围 => start:end
                String peakStartToEnd = peakRange.getKey();
                // 某个m6A信号下 major和 minor haplotype SNP reads数目
                HashMap<String, LinkedList<Integer>> haplotypeSnpReads = peakRange.getValue();

                LinkedList<Integer> majorSNPReadsCount = haplotypeSnpReads.get("major");
                int[] majorCount = new int[majorSNPReadsCount.size()];
                for (int i = 0 ; i < majorSNPReadsCount.size(); i++ ) {
                    majorCount[i] = majorSNPReadsCount.get(i);
                }
                majorSNPReadsCount = null;
                LinkedList<Integer> minorSNPReadsCount = haplotypeSnpReads.get("minor");
                int[] minorCount = new int[minorSNPReadsCount.size()];
                for (int i = 0 ; i < minorSNPReadsCount.size(); i++ ) {
                    minorCount[i] = minorSNPReadsCount.get(i);
                }
                minorSNPReadsCount = null;
                // 对peak下所有的SNP位点进行元分析, 得到该peak对应的 p value
                HierarchicalBayesianModel hb = new HierarchicalBayesianModel(0, 1, this.samplingTime,
                        this.burnIn, majorCount, minorCount);
                pVal = hb.testSignificant();
                this.asmPValue.put(chrNum+":"+peakStartToEnd, pVal);
                hb = null;
            }
        }
    }

    /**
     * 对Peak ASM的p值进行BH校正, 得到 Q value.
     */
    private void bhRecalibrationOfEachPeak() {
        List<Map.Entry<String, Double>> sortedByValue = new ArrayList<Map.Entry<String, Double>>(this.asmPValue.entrySet());
        // 依据ASE的p值从小到大排序
        Collections.sort(sortedByValue, new Comparator<Map.Entry<String, Double>>() {
            public int compare(Map.Entry<String, Double> o1,
                               Map.Entry<String, Double> o2) {
                return (o2.getValue()).compareTo(o1.getValue());
            }
        });

        int rankage = 0;
        int totalPeak = sortedByValue.size();
        double pValue, qValue, prevQVal = Collections.max(this.asmPValue.values());
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
            this.asmQValue.put(label, qValue);
        }
    }

    /**
     * 将ASM的检验结果输出到文件
     */
    private void outputResult() {
        List<Map.Entry<String, Double>> sortedByValue = new ArrayList<Map.Entry<String, Double>>(this.asmQValue.entrySet());
        // 依据ASE的q值从小到大排序
        Collections.sort(sortedByValue, new Comparator<Map.Entry<String, Double>>() {
            public int compare(Map.Entry<String, Double> o1,
                               Map.Entry<String, Double> o2) {
                return (o1.getValue()).compareTo(o2.getValue());
            }
        });
        BufferedWriter bfw = null;
        try {
            bfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(this.asmPeakFile))));
            String line, label, geneId, chrNum, peakStart, peakEnd, majorAlleleReads, minorAlleleReads, majorAlleleStrand;
            String[] info;
            double qVal;
            bfw.write("#chr\tpeakStart\tpeakEnd\tgeneId\tq-value\tmajorAlleleReads\tminorAlleleReads\tmajorAlleleStrand\n");
            for (Map.Entry<String, Double> record: sortedByValue) {
                label = record.getKey();
                // info = [chrNum, peakStart, peakEnd]
                info = label.split(":");
                chrNum = info[0];
                peakStart = info[1];
                peakEnd = info[2];
                geneId = this.peakCoveredGene.get(label);
                majorAlleleStrand = this.peakMajorAlleleStrand.get(label);
                qVal = record.getValue();
                majorAlleleReads = this.peakSnpReadsCount.get(chrNum).get(peakStart+":"+peakEnd).get("major").toString();
                minorAlleleReads = this.peakSnpReadsCount.get(chrNum).get(peakStart+":"+peakEnd).get("minor").toString();
                line = String.join("\t", new String[]{chrNum, peakStart, peakEnd, geneId, this.df.format(qVal),
                                                                 majorAlleleReads, minorAlleleReads, majorAlleleStrand});
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

    /**
     * 初始化log4j Logger 对象
     * @param logHome output directory of log file
     * @return Logger instance
     */
    private static Logger initLog(String logHome) {
        System.setProperty("log_home", logHome);
        return Logger.getLogger(AsmPeakDetection.class);
    }

}
