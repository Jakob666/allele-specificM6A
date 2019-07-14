package HierarchicalBayesianAnalysis;

import AseM6aPeakDetector.HeterozygoteReadsCount;
import org.apache.log4j.Logger;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

public class AseSpecificM6aSignal {
    private String peakCoveredSNPFile;
    private File aseM6aPeakFile;
    private Logger log;
    private DecimalFormat df = new DecimalFormat("0.0000");
    private HashMap<String, HashMap<String, HashMap<String, LinkedList<Integer>>>> peakSnpReadsCount;

    public AseSpecificM6aSignal(String peakCoveredSNPFile, Logger log) {
        this.peakCoveredSNPFile = peakCoveredSNPFile;
        String outputFileName = peakCoveredSNPFile.substring(0, peakCoveredSNPFile.lastIndexOf("_")) + "_asePeakBayesian.txt";
        this.aseM6aPeakFile = new File(outputFileName);
        this.log = log;
    }

    /**
     * 检验ASE特异的m6A信号并将其写入文件
     */
    public void detectAsePeak() {
        // 获取beta-binomial分布全局离散度rho
        BufferedWriter bfw = null;
        try {
            bfw = new BufferedWriter(
                    new OutputStreamWriter(new FileOutputStream(this.aseM6aPeakFile))
            );
            // 得到每个m6A信号的p值，进行BH校正。 [chrNum:peakStart:peakEnd -> pValue, ...]
            HashMap<String, Double> peakPValues = this.getPValueOfEachPeak();
            HashMap<String, HashMap<String, Double>> aseM6aPeaks = this.bhRecalibrationOfEachPeak(peakPValues);

            // 将结果写入文件
            Set<Map.Entry<String, HashMap<String, Double>>> records = aseM6aPeaks.entrySet();
            // 依据peak的p值从小到大排序
            ArrayList<Map.Entry<String, HashMap<String, Double>>> peakResult = new ArrayList<>(aseM6aPeaks.entrySet());
            Collections.sort(peakResult, new Comparator<Map.Entry<String, HashMap<String, Double>>>() {
                public int compare(Map.Entry<String, HashMap<String, Double>> o1,
                                   Map.Entry<String, HashMap<String, Double>> o2) {
                    return (o1.getValue().get("qVal")).compareTo(o2.getValue().get("qVal"));
                }
            });
            double pVal, qVal;
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
                String majorAlleleReads = this.peakSnpReadsCount.get(chrNum).get(peakStart+":"+peakEnd).get("major").toString();
                String minorAlleleReads = this.peakSnpReadsCount.get(chrNum).get(peakStart+":"+peakEnd).get("minor").toString();
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

    /**
     * 获取每个m6A信号范围内ASE位点上的 major allele和 minor allele的reads count
     * @return [chr1: [peak1: [major: [count1, count2,...], minor: [count1, count2,...]]], chr2:....]
     */
    private HashMap<String, HashMap<String, HashMap<String, LinkedList<Integer>>>> getPeakSNPReadsCount() {
        HeterozygoteReadsCount hrc = new HeterozygoteReadsCount(this.peakCoveredSNPFile, this.log);
        // chr: peakRange: major/minor: readsCounts
        HashMap<String, HashMap<String, HashMap<String, LinkedList<Integer>>>> majorMinorHaplotype = hrc.getMajorMinorHaplotype();
        hrc = null;

        return majorMinorHaplotype;
    }

    /**
     * 计算每个m6A信号等位基因特异性的p值
     * @return [chrNum:peakStart:peakEnd -> pValue, ...]
     */
    private HashMap<String, Double> getPValueOfEachPeak() {
        double pVal;
        HashMap<String, Double> pValues = new HashMap<>();

        // [chr1: [peak1: [major: [count1, count2,...], minor: [count1, count2,...]]], chr2:....]
        this.peakSnpReadsCount = this.getPeakSNPReadsCount();
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
}
