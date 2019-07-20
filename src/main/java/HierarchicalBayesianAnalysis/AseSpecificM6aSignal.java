package HierarchicalBayesianAnalysis;

import AseM6aPeakDetector.HeterozygoteReadsCount;
import org.apache.log4j.Logger;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

public class AseSpecificM6aSignal {
    private String ipPeakCoveredSNPFile, inputPeakCoveredSNPFile;
    private File ipAseM6aPeakFile, inputAseM6aPeakFile;
    private HashMap<String, HashMap<String, HashMap<String, LinkedList<Integer>>>> inputPeakSnpReadsCount, ipPeakSnpReadsCount;
    private Logger log;
    private DecimalFormat df = new DecimalFormat("0.0000");

    public AseSpecificM6aSignal(String inputPeakCoveredSNPFile, String ipPeakCoveredSNPFile, Logger log) {
        this.inputPeakCoveredSNPFile = inputPeakCoveredSNPFile;
        this.ipPeakCoveredSNPFile = ipPeakCoveredSNPFile;
        String inputOutputFileName = inputPeakCoveredSNPFile.substring(0, inputPeakCoveredSNPFile.lastIndexOf("_")) + "_asePeakBayesian.txt";
        String ipOutputFileName = ipPeakCoveredSNPFile.substring(0, ipPeakCoveredSNPFile.lastIndexOf("_")) + "_asmPeakBayesian.txt";
        this.inputAseM6aPeakFile = new File(inputOutputFileName);
        this.ipAseM6aPeakFile = new File(ipOutputFileName);
        this.log = log;
    }

    /**
     * 检验ASE特异的m6A信号并将其写入文件
     */
    public void detectAsePeak() {
        // 得到每个m6A信号的p值，进行BH校正。 [chrNum:peakStart:peakEnd -> pValue, ...]
        HashMap<String, HashMap<String, Double>> modelTestResult = this.getPValueOfEachPeak();
        HashMap<String, Double> ipPeakPValues = modelTestResult.get("ip");
        HashMap<String, Double> inputPeakPValues = modelTestResult.get("input");
        HashMap<String, HashMap<String, Double>> ipAseM6aPeaks = this.bhRecalibrationOfEachPeak(ipPeakPValues);
        HashMap<String, HashMap<String, Double>> inputAseM6aPeaks = this.bhRecalibrationOfEachPeak(inputPeakPValues);

        this.outputResult(this.inputAseM6aPeakFile, inputAseM6aPeaks, this.inputPeakSnpReadsCount);
        this.outputResult(this.ipAseM6aPeakFile, ipAseM6aPeaks, this.ipPeakSnpReadsCount);
    }

    /**
     * 获取每个m6A信号范围内ASE位点上的 major allele和 minor allele的reads count
     * @return [chr1: [peak1: [major: [count1, count2,...], minor: [count1, count2,...]]], chr2:....]
     */
    private HashMap<String, HashMap<String, HashMap<String, LinkedList<Integer>>>> getPeakSNPReadsCount(String peakCoveredSNPFile) {
        HeterozygoteReadsCount hrc = new HeterozygoteReadsCount(peakCoveredSNPFile, this.log);
        // chr: peakRange: major/minor: readsCounts

        return hrc.getMajorMinorHaplotype();
    }

    /**
     * 计算每个m6A信号等位基因特异性的p值
     * @return [chrNum:peakStart:peakEnd -> pValue, ...]
     */
    private HashMap<String, HashMap<String, Double>> getPValueOfEachPeak() {
        HashMap<String, Double> ipPValues, inputPValues;


        // [chr1: [peak1: [major: [count1, count2,...], minor: [count1, count2,...]]], chr2:....]
        this.inputPeakSnpReadsCount = this.getPeakSNPReadsCount(this.inputPeakCoveredSNPFile);
        this.ipPeakSnpReadsCount = this.getPeakSNPReadsCount(this.ipPeakCoveredSNPFile);
        inputPValues = this.hierarchicalModelTest(this.inputPeakSnpReadsCount);
        ipPValues = this.hierarchicalModelTest(this.ipPeakSnpReadsCount);

        HashMap<String, HashMap<String, Double>> result = new HashMap<>();
        result.put("ip", ipPValues);
        result.put("input", inputPValues);

        return result;
    }

    /**
     * 使用模型对IP和INPUT样本进行检验分别得到ASE和ASM的peak
     * @param peakSnpReadsCount 每个m6A信号下major allele和minor allele的统计
     * @return 各m6A信号的显著性p值
     */
    private HashMap<String, Double> hierarchicalModelTest(HashMap<String, HashMap<String, HashMap<String, LinkedList<Integer>>>> peakSnpReadsCount) {
        HashMap<String, Double> pValues = new HashMap<>();
        double pVal;
        Set<Map.Entry<String, HashMap<String, HashMap<String, LinkedList<Integer>>>>> chrPeaks = peakSnpReadsCount.entrySet();
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
                HierarchicalBayesianModel hb = new HierarchicalBayesianModel(0, 1, 5000,
                        200, majorCount, minorCount);
                pVal = hb.testSignificant();
                pValues.put(chrNum+":"+peakStartToEnd, pVal);
                hb = null;
            }
        }

        return pValues;
    }

    /**
     * 对m6A信号的p值进行BH校正, 得到 Q value.
     * @param peakPValues 每个m6A信号的p值
     * @return [chrNum:peakStart:peakEnd : [pVal: xxx, qVal: xxx], ...]
     */
    private HashMap<String, HashMap<String, Double>> bhRecalibrationOfEachPeak(HashMap<String, Double> peakPValues) {
        List<Map.Entry<String, Double>> sortedByValue = new ArrayList<Map.Entry<String, Double>>(peakPValues.entrySet());
        // 依据peak的p值从小到大排序
        Collections.sort(sortedByValue, new Comparator<Map.Entry<String, Double>>() {
            public int compare(Map.Entry<String, Double> o1,
                               Map.Entry<String, Double> o2) {
                return (o2.getValue()).compareTo(o1.getValue());
            }
        });

        HashMap<String, HashMap<String, Double>> aseM6aPeaks = new HashMap<>();
        int rankage = 0;
        int totalPeak = sortedByValue.size();
        double peakPValue, qValue, prevQVal = Collections.max(peakPValues.values());
        String peakLabel;
        for (Map.Entry<String, Double> peakRecord: sortedByValue) {
            rankage += 1;
            peakLabel = peakRecord.getKey();
            peakPValue = peakRecord.getValue();
            qValue = peakPValue * totalPeak / rankage;
            if (qValue > prevQVal) {
                qValue = prevQVal;
            } else {
                prevQVal = Math.min(1.0, qValue);
            }
            HashMap<String, Double> pAndQValue = new HashMap<>();
            pAndQValue.put("qVal", qValue);
            aseM6aPeaks.put(peakLabel, pAndQValue);
        }

        return aseM6aPeaks;
    }

    /**
     * 将检验结果写入到文件
     * @param outputResultFile 输出文件
     * @param aseM6aPeaks 每个peak经层次模型检验后的Q值
     * @param peakSnpReadsCount 每个peak覆盖范围内major allele和minor allele的reads记录
     */
    private void outputResult(File outputResultFile, HashMap<String, HashMap<String, Double>> aseM6aPeaks, HashMap<String,
                              HashMap<String, HashMap<String, LinkedList<Integer>>>> peakSnpReadsCount) {
        BufferedWriter bfw = null;
        try {
            bfw = new BufferedWriter(
                    new OutputStreamWriter(new FileOutputStream(outputResultFile))
            );

            // 依据peak的p值从小到大排序
            ArrayList<Map.Entry<String, HashMap<String, Double>>> peakResult = new ArrayList<>(aseM6aPeaks.entrySet());
            Collections.sort(peakResult, new Comparator<Map.Entry<String, HashMap<String, Double>>>() {
                public int compare(Map.Entry<String, HashMap<String, Double>> o1,
                                   Map.Entry<String, HashMap<String, Double>> o2) {
                    return (o1.getValue().get("qVal")).compareTo(o2.getValue().get("qVal"));
                }
            });
            double qVal;
            String peakLabel, writeOut;
            String[] info, outputLine;
            bfw.write("#chr\tpeakStart\tpeakEnd\tq-value\tmajorAlleleReads\tminorAlleleReads\n");
            for (Map.Entry<String, HashMap<String, Double>> record: peakResult) {
                peakLabel = record.getKey();
                info = peakLabel.split(":");
                String chrNum = info[0];
                String peakStart = info[1];
                String peakEnd = info[2];
                qVal = record.getValue().get("qVal");
                String majorAlleleReads = peakSnpReadsCount.get(chrNum).get(peakStart+":"+peakEnd).get("major").toString();
                String minorAlleleReads = peakSnpReadsCount.get(chrNum).get(peakStart+":"+peakEnd).get("minor").toString();
                outputLine = new String[]{chrNum, peakStart, peakEnd, this.df.format(qVal),
                        majorAlleleReads, minorAlleleReads};
                writeOut = String.join("\t", outputLine);
                bfw.write(writeOut);
                bfw.newLine();
            }
            bfw.close();
        } catch (IOException ie) {
            this.log.error("can not write ASE m6a peak result file");
            this.log.error(ie.getMessage());
        } finally {
            if (bfw != null) {
                try {
                    bfw.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }
}
