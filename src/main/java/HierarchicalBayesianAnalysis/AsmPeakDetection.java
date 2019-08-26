package HierarchicalBayesianAnalysis;


import AseM6aPeakDetector.HeterozygoteReadsCount;
import org.apache.log4j.Logger;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

public class AsmPeakDetection {
    private String peakBedFile, vcfFile, wesFile, asmPeakFile, peakCoveredSnpFile, peakCoveredWesSnpFile;
    private int samplingTime, burnIn;
    private Logger log;
    private HashMap<String, HashMap<String, HashMap<String, int[]>>> peakSnpReadsCount, peakSnpBackground;
    private HashMap<String, Double> asmPValue = new HashMap<>(), asmQValue = new HashMap<>();
    private HashMap<String, String> peakCoveredGene = new HashMap<>(), peakMajorAlleleStrand;
    private DecimalFormat df = new DecimalFormat("0.0000");

    /**
     * Constructor
     * @param peakBedFile peak calling得到的BED格式文件
     * @param vcfFile RNA-seq数据SNP Calling得到的VCF格式文件
     * @param wesFile WES数据SNP Calling得到的VCF格式文件
     * @param peakCoveredSnpFile 记录peak信号覆盖的RNA-seq SNP位点的文件
     * @param peakCoveredWesSnpFile 记录peak信号覆盖的WES SNP位点的文件
     * @param asmPeakFile ASM peak检验结果输出文件
     * @param samplingTime 采样次数
     * @param burnIn burn in次数
     */
    public AsmPeakDetection(String peakBedFile, String vcfFile, String wesFile, String peakCoveredSnpFile, String peakCoveredWesSnpFile,
                            String asmPeakFile, int samplingTime, int burnIn) {
        this.peakBedFile = peakBedFile;
        this.vcfFile = vcfFile;
        this.wesFile = wesFile;
        this.peakCoveredSnpFile = peakCoveredSnpFile;
        this.peakCoveredWesSnpFile = peakCoveredWesSnpFile;
        if (this.wesFile != null)
            assert this.peakCoveredWesSnpFile != null;
        this.asmPeakFile = asmPeakFile;
        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
        this.log = initLog(new File(asmPeakFile).getParent());
    }

    public void getTestResult() {
        this.getPeakCoveredGene();
        this.getPeakCoveredSnpResult();
        if (this.wesFile != null)
            this.getPeakCoveredSnpBackground();
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
     * 当存在WES数据SNP calling结果时，将结果写入文件
     */
    private void getPeakCoveredSnpBackground() {
        if (this.wesFile != null) {
            PeakCoveredSNPRecord pcsr = new PeakCoveredSNPRecord(this.wesFile, this.peakBedFile, this.peakCoveredWesSnpFile);
            pcsr.getPeakCoveredSNP();
        }
    }

    /**
     * 获取RNA-seq得到的每个m6A信号范围内SNV位点上的 major allele和 minor allele的reads count
     * @return [chr1: [peak1: position1: [majorCount, minorCount], position2:[majorCount, minorCount]]], chr2:....]
     */
    private HashMap<String, HashMap<String, HashMap<String, int[]>>> getPeakSNPReadsCount() {
        HeterozygoteReadsCount hrc = new HeterozygoteReadsCount(this.peakCoveredSnpFile, this.log);
        this.peakMajorAlleleStrand = hrc.getPeakMajorAlleleStrand();
        return hrc.getMajorMinorHaplotype();
    }

    /**
     * 获取WES得到的每个m6A信号范围内SNV位点上的 major allele和 minor allele的reads count
     * @return [chr1: [peak1: position1: [majorCount, minorCount], position2:[majorCount, minorCount]]], chr2:....]
     */
    private HashMap<String, HashMap<String, HashMap<String, int[]>>> getPeakSNPBackground() {
        HeterozygoteReadsCount hrc = new HeterozygoteReadsCount(this.peakCoveredWesSnpFile, this.log);
        return hrc.getMajorMinorHaplotype();
    }

    /**
     * 对每个覆盖SNP的peak进行检验
     */
    private void asmPeakTest() {
        // [chr1: [peak1: position1: [majorCount, minorCount], position2:[majorCount, minorCount]]], chr2:....]
        this.peakSnpReadsCount = this.getPeakSNPReadsCount();
        if (this.wesFile != null)
            this.peakSnpBackground = this.getPeakSNPBackground();
        else
            this.peakSnpBackground = null;
        this.hierarchicalModelTest();
    }

    /**
     * 对每个覆盖了SNP的peak进行检验，得到ASM显著性p值
     */
    private void hierarchicalModelTest() {
        double pVal;
        Set<Map.Entry<String, HashMap<String, HashMap<String, int[]>>>> rnaSeqChrPeaks = this.peakSnpReadsCount.entrySet();
        Set<Map.Entry<String, HashMap<String, HashMap<String, int[]>>>> wesSeqChrPeaks = (this.peakSnpBackground != null)? this.peakSnpBackground.entrySet():null;
        HashMap<String, HashMap<String, int[]>> rnaSeqPeakSnvAlleleReads, wesPeakSnvAlleleReads;
        HashMap<String, int[]> rnaSeqMutPositionAlleleReads, wesMutPositionAlleleReads;
        int[] rnaSeqReads, wesReads, majorCount, minorCount, majorBackground, minorBackground;
        int major, minor, backgroundMajor, backgroundMinor;
        ArrayList<Integer> rnaSeqMajor, rnaSeqMinor, wesMajor, wesMinor;
        for (String chrNum: this.peakSnpReadsCount.keySet()) {
            // 某条染色体上所有Peak记录
            rnaSeqPeakSnvAlleleReads = this.peakSnpReadsCount.get(chrNum);

            if (wesSeqChrPeaks != null)
                wesPeakSnvAlleleReads = this.peakSnpBackground.getOrDefault(chrNum, null);
            else
                wesPeakSnvAlleleReads = null;

            for (String peakRange: rnaSeqPeakSnvAlleleReads.keySet()) {
                // 染色体上某个Peak覆盖的SNV记录
                rnaSeqMutPositionAlleleReads = rnaSeqPeakSnvAlleleReads.get(peakRange);
                if (wesPeakSnvAlleleReads != null)
                    wesMutPositionAlleleReads = wesPeakSnvAlleleReads.getOrDefault(peakRange, null);
                else
                    wesMutPositionAlleleReads = null;

                rnaSeqMajor = new ArrayList<>();
                rnaSeqMinor = new ArrayList<>();
                wesMajor = new ArrayList<>();
                wesMinor = new ArrayList<>();
                for (String position: rnaSeqMutPositionAlleleReads.keySet()) {
                    rnaSeqReads = rnaSeqMutPositionAlleleReads.get(position);
                    if (wesMutPositionAlleleReads != null)
                        wesReads = wesMutPositionAlleleReads.getOrDefault(position, null);
                    else
                        wesReads = null;
                    major = rnaSeqReads[0];
                    minor = rnaSeqReads[1];
                    backgroundMajor = (wesReads != null)? wesReads[0]: (minor+major)/2;
                    backgroundMinor = (wesReads != null)? wesReads[1]: (minor+major)/2;
                    rnaSeqMajor.add(major);
                    rnaSeqMinor.add(minor);
                    wesMajor.add(backgroundMajor);
                    wesMinor.add(backgroundMinor);
                }
                majorCount = new int[rnaSeqMajor.size()];
                minorCount = new int[rnaSeqMinor.size()];
                majorBackground = new int[wesMajor.size()];
                minorBackground = new int[wesMinor.size()];
                assert majorCount.length == minorCount.length;
                assert majorCount.length == majorBackground.length;
                assert majorBackground.length == minorBackground.length;
                for (int i=0; i<majorCount.length; i++) {
                    majorCount[i] = rnaSeqMajor.get(i);
                    minorCount[i] = rnaSeqMinor.get(i);
                    majorBackground[i] = wesMajor.get(i);
                    minorBackground[i] = wesMinor.get(i);
                }
                rnaSeqMajor.clear();
                rnaSeqMinor.clear();
                wesMajor.clear();
                wesMinor.clear();
                // 对peak下所有的SNP位点进行元分析, 得到该peak对应的 p value
                HierarchicalBayesianModel hb = new HierarchicalBayesianModel(0, 1, this.samplingTime,
                        this.burnIn, majorCount, minorCount, majorBackground, minorBackground);
                pVal = hb.testSignificant();
                this.asmPValue.put(chrNum+":"+peakRange, pVal);
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
