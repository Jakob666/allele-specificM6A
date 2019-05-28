package AseM6aPeakDetector;

import betaBinomialMetaAnalysis.RhoEstimator;
import betaBinomialMetaAnalysis.SignificantTest;
import betaBinomialMetaAnalysis.betaBinomialDistribution;
import heterozygoteSiteAnalysis.PeakCoveredSNP;
import org.apache.commons.cli.*;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.*;

public class AseM6aPeakDetector {
    private String peakCoveredSNPFile;
    private File aseM6aPeakFile;
    private double initialRho, learningRate, unImproveThreshold;
    private Logger log;

    public static void main(String[] args) throws ParseException{
        Options options = new Options();
        CommandLine commandLine = setCommand(args, options);
        String vcfFile, peakFile;
        double rho = 0.01, lr = 0.0005, unImproveThreshold = 0.000001;

        vcfFile = commandLine.getOptionValue("v");
        peakFile = commandLine.getOptionValue("p");
        if (commandLine.hasOption("rho"))
            rho = Double.parseDouble(commandLine.getOptionValue("rho"));
        if (commandLine.hasOption("lr"))
            lr = Double.parseDouble(commandLine.getOptionValue("lr"));

        Logger log = initLog(new File(vcfFile).getParent());

        PeakCoveredSNP pcs = new PeakCoveredSNP(vcfFile, peakFile, log);
        pcs.filterSNPAndPeak();
        String peakCoverFile = vcfFile.substring(0, vcfFile.lastIndexOf("_")) + "_peakCoveredSNP.txt";

        AseM6aPeakDetector aseM6aPeakDetector = new AseM6aPeakDetector(peakCoverFile, rho, lr, unImproveThreshold, log);
        aseM6aPeakDetector.detectAsePeak();

    }

    public AseM6aPeakDetector(String peakCoveredSNPFile, double initialRho, double learningRate,
                              double unImproveThreshold, Logger log) {
        this.peakCoveredSNPFile = peakCoveredSNPFile;
        this.initialRho = initialRho;
        this.learningRate = learningRate;
        this.unImproveThreshold = unImproveThreshold;
        this.log = log;
        String outputFileName = peakCoveredSNPFile.substring(0, peakCoveredSNPFile.lastIndexOf("_")) + "_asePeak.txt";
        this.aseM6aPeakFile = new File(outputFileName);
    }

    /**
     * 检验ASE特异的m6A信号并将其写入文件
     */
    public void detectAsePeak() {
        // 获取beta-binomial分布全局离散度rho
        double rho = this.getGlobalRho();
        BufferedWriter bfw = null;
        try {
            bfw = new BufferedWriter(
                    new OutputStreamWriter(new FileOutputStream(this.aseM6aPeakFile))
            );
            // 得到每个m6A信号的p值，进行BH校正。 [chrNum:peakStart:peakEnd -> pValue, ...]
            HashMap<String, Double> peakPValues = this.getPValueOfEachPeak(rho);
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
            bfw.write("#chr\tpeakStart\tpeakEnd\tp-value\tq-value\n");
            for (Map.Entry<String, HashMap<String, Double>> record: peakResult) {
                peakLabel = record.getKey();
                info = peakLabel.split(":");
                String chrNum = info[0];
                String peakStart = info[1];
                String peakEnd = info[2];
                pVal = record.getValue().get("pVal");
                qVal = record.getValue().get("qVal");
                outputLine = new String[]{chrNum, peakStart, peakEnd, Double.toString(pVal), Double.toString(qVal)};
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
     * 获取每个m6A信号在 major haplotype和 minor haplotype上ASE位点的reads count
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
     * 得到 major haplotype和 minor haplotype上每个ASE位点的reads count
     * @return [major: [count1, count2..., countN], minor: [count1, count2..., countN]]
     */
    private HashMap<String, LinkedList<Integer>> getHaplotypeSNPReadsCount() {
        HaplotypeSNPReadsCount hsrc = new HaplotypeSNPReadsCount(this.peakCoveredSNPFile, log);
        HashMap<String, LinkedList<Integer>> haplotypeSNPReadsCount = hsrc.haplotypeSnpReadsCount();
        hsrc = null;

        return haplotypeSNPReadsCount;
    }

    /**
     * get the overdispersion value rho of the beta binomial distribution
     * @return rho
     */
    private double getGlobalRho() {
        HashMap<String, LinkedList<Integer>> haplotypeSNPReadsCount = this.getHaplotypeSNPReadsCount();
        // major haplotype上所有major allele
        LinkedList<Integer> majorHaplotype = haplotypeSNPReadsCount.get("major");
        int[] majorSNPReadsCount = new int[majorHaplotype.size()];
        for (int i = 0 ; i < majorHaplotype.size(); i++) {
            majorSNPReadsCount[i] = majorHaplotype.get(i);
        }
        majorHaplotype = null;
        // minor haplotype上所有minor allele
        LinkedList<Integer> minorHaplotype = haplotypeSNPReadsCount.get("minor");
        int[] minorSNPReadsCount = new int[minorHaplotype.size()];
        for (int i = 0 ; i < minorHaplotype.size(); i++ ) {
            minorSNPReadsCount[i] = minorHaplotype.get(i);
        }
        minorHaplotype = null;
        // 初步计算beta-binomial分布的离散度rho
        RhoEstimator re = new RhoEstimator(majorSNPReadsCount, minorSNPReadsCount, this.initialRho,
                                           this.learningRate, this.unImproveThreshold);
        double rho = re.gradientAscend();
//        System.out.println("first rho: " + rho);

        // 使用目前的rho计算每个位点major SNP概率
        betaBinomialDistribution bd = new betaBinomialDistribution();
        int major, minor, total;
        double fMajor, pValue;
        HashMap<String, Double> pValueList = new HashMap<>();
        fMajor = (double) getSum(majorSNPReadsCount) / (double) (getSum(majorSNPReadsCount) + getSum(minorSNPReadsCount));
        for (int i = 0; i < majorSNPReadsCount.length; i++) {
            major = majorSNPReadsCount[i];
            minor = minorSNPReadsCount[i];
            total = major + minor;
            pValue = 1.0 - bd.betaBinomialCdf(0, major-1, total, fMajor, rho);
//            System.out.println("major haplotype: "+major+"\tminor haplotype: "+minor+"\tp:"+pValue);
            pValueList.put(Integer.toString(i), pValue);
        }
        // BH校正后除去p值小于0.05的SNP位点重新求取beta-binomial分布的离散度rho
        this.bhRecalibrationOfSNP(pValueList);
        Set<String> idxs = pValueList.keySet();
        int[] newMajorSNPReadsCount = new int[idxs.size()];
        int[] newMinorSNPReadsCount = new int[idxs.size()];
        int i = 0;
        for (String idx: idxs) {
            int index = Integer.parseInt(idx);
            newMajorSNPReadsCount[i] = majorSNPReadsCount[index];
            newMinorSNPReadsCount[i] = minorSNPReadsCount[index];
            i++;
        }
        majorSNPReadsCount = null;
        minorSNPReadsCount = null;
        re = new RhoEstimator(newMajorSNPReadsCount, newMinorSNPReadsCount, this.initialRho,
                              this.learningRate, this.unImproveThreshold);
        rho = re.gradientAscend();
//        System.out.println(rho);

        newMajorSNPReadsCount = null;
        newMinorSNPReadsCount = null;

        return rho;
    }

    /**
     * 计算每个m6A信号等位基因特异性的p值
     * @param rho beta binomial分布的全局离散度
     * @return [chrNum:peakStart:peakEnd -> pValue, ...]
     */
    private HashMap<String, Double> getPValueOfEachPeak(double rho) {
        double pVal;
        SignificantTest st;
        HashMap<String, Double> pValues = new HashMap<>();

        // [chr1: [peak1: [major: [count1, count2,...], minor: [count1, count2,...]]], chr2:....]
        HashMap<String, HashMap<String, HashMap<String, LinkedList<Integer>>>> peakSnpReadsCount = this.getPeakSNPReadsCount();
        Set<Map.Entry<String, HashMap<String, HashMap<String, LinkedList<Integer>>>>> chrPeaks = peakSnpReadsCount.entrySet();
        for (Map.Entry<String, HashMap<String, HashMap<String, LinkedList<Integer>>>> chrPeak: chrPeaks) {
            String chrNum = chrPeak.getKey();
            // 某条染色体的peak的记录
            HashMap<String, HashMap<String, LinkedList<Integer>>> peaksOnChr = chrPeak.getValue();

            Set<Map.Entry<String, HashMap<String, LinkedList<Integer>>>> peakRanges = peaksOnChr.entrySet();
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
                st = new SignificantTest(majorCount, minorCount, rho);
                pVal = st.testSignificant();
                pValues.put(chrNum+":"+peakStartToEnd, pVal);
                st = null;
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
                return (o1.getValue()).compareTo(o2.getValue());
            }
        });

        HashMap<String, HashMap<String, Double>> aseM6aPeaks = new HashMap<>();
        int rankage = 0;
        int totalPeak = sortedByValue.size();
        double peakPValue;
        String peakLabel;
        for (Map.Entry<String, Double> peakRecord: sortedByValue) {
            rankage += 1;
            peakLabel = peakRecord.getKey();
            peakPValue = peakRecord.getValue();
            double qValue = SignificantTest.BHRecalibration(peakPValue, rankage, totalPeak);
            HashMap<String, Double> pAndQValue = new HashMap<>();
            pAndQValue.put("pVal", peakPValue);
            pAndQValue.put("qVal", qValue);
            aseM6aPeaks.put(peakLabel, pAndQValue);
        }

        return aseM6aPeaks;
    }

    /**
     * 对SNP位点的p value进行BH校正
     * @param snpPValues [SNP位点序号: p value, ...]
     */
    private void bhRecalibrationOfSNP(HashMap<String, Double> snpPValues) {
        List<Map.Entry<String, Double>> sortedByValue = new ArrayList<Map.Entry<String, Double>>(snpPValues.entrySet());
        // 通过 p values 大小排序，升序排列
        Collections.sort(sortedByValue, new Comparator<Map.Entry<String, Double>>() {
            public int compare(Map.Entry<String, Double> o1,
                               Map.Entry<String, Double> o2) {
                return (o1.getValue()).compareTo(o2.getValue());
            }
        });
        int rankage = 0, totalSnp = sortedByValue.size();
        double snpPValue;
        String snpLabel;
        for (Map.Entry<String, Double> entry: sortedByValue) {
            rankage++;
            snpLabel = entry.getKey();
            snpPValue = entry.getValue();
            double snpQValue = SignificantTest.BHRecalibration(snpPValue, rankage, totalSnp);
            if (snpPValue < 0.05)
                snpPValues.remove(snpLabel);
            else
                snpPValues.put(snpLabel, snpQValue);
        }
    }

    private int getSum(int[] SNPReads) {
        int sum = 0;
        for (int i: SNPReads) {
            sum = sum + i;
        }
        return sum;
    }

    private static CommandLine setCommand(String[] args, Options options) throws ParseException {
        Option option = new Option("v", "vcf_file", true, "SNP calling VCF file path");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("p", "peak_file", true, "peak calling BED file path");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("rho", "initial_rho", true, "initial rho of the beta-binomial distribution, default 0.001");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("lr", "learning_rate", true, "learning rate for calculating global rho");
        option.setRequired(false);
        options.addOption(option);

        CommandLineParser parser = new DefaultParser();

        return parser.parse(options, args);
    }

    /**
     * 初始化log4j Logger 对象
     * @param logHome output directory of log file
     * @return Logger instance
     */
    private static Logger initLog(String logHome) {
        System.setProperty("log_home", logHome);
        return Logger.getLogger(AseM6aPeakDetector.class);
    }
}
