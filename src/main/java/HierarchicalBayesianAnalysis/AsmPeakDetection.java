package HierarchicalBayesianAnalysis;


import AseM6aPeakDetector.HeterozygoteReadsCount;
import org.apache.commons.cli.*;
import org.apache.log4j.Logger;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

public class AsmPeakDetection {
    private String peakBedFile, vcfFile, wesFile, asmPeakFile, peakCoveredSnpFile, peakCoveredWesSnpFile;
    private int ipSNPReadInfimum, wesSNPReadInfimum, samplingTime, burnIn;
    double tauInfimum, tauSupremum;
    private Logger log;
    private HashMap<String, HashMap<String, HashMap<String, HashMap<String, Integer>>>> peakSnpReadsCount, peakSnpBackground;
    private HashMap<Double, ArrayList<String>> asmPValue = new HashMap<>();
    private HashMap<String, HashMap<String, Integer>> peakMajorMinorAlleleCount = new HashMap<>(), peakMajorMinorBackground = new HashMap<>();
    private HashMap<String, Double> peakMajorAlleleFrequency = new HashMap<>();
    private HashMap<String, String> peakCoveredGene = new HashMap<>(), peakMajorAlleleStrand;
    private HashMap<String, Integer> peakSNVNum = new HashMap<>();
    private ArrayList<String> asmQValue = new ArrayList<>();
    private DecimalFormat df = new DecimalFormat("0.0000");

    /**
     * Constructor
     * @param peakBedFile peak calling得到的BED格式文件
     * @param vcfFile RNA-seq数据SNP Calling得到的VCF格式文件
     * @param wesFile WES数据SNP Calling得到的VCF格式文件
     * @param peakCoveredSnpFile 记录peak信号覆盖的RNA-seq SNP位点的文件
     * @param peakCoveredWesSnpFile 记录peak信号覆盖的WES SNP位点的文件
     * @param asmPeakFile ASM peak检验结果输出文件
     * @param tauInfimum tau采样时的下确界
     * @param tauSupremum tau采样时的上确界
     * @param ipSNPReadInfimum 记录IP样本Peak覆盖的SNP位点时，对SNP位点筛选所用的阈值
     * @param wesSNPReadInfimum 记录WES样本Peak覆盖的SNP位点时，对SNP位点筛选所用的阈值
     * @param samplingTime 采样次数
     * @param burnIn burn in次数
     */
    public AsmPeakDetection(String peakBedFile, String vcfFile, String wesFile, String peakCoveredSnpFile,
                            String peakCoveredWesSnpFile, String asmPeakFile, double tauInfimum, double tauSupremum,
                            int ipSNPReadInfimum, int wesSNPReadInfimum, int samplingTime, int burnIn) {
        this.peakBedFile = peakBedFile;
        this.vcfFile = vcfFile;
        this.wesFile = wesFile;
        this.peakCoveredSnpFile = peakCoveredSnpFile;
        this.peakCoveredWesSnpFile = peakCoveredWesSnpFile;
        if (this.wesFile != null)
            assert this.peakCoveredWesSnpFile != null;
        this.asmPeakFile = asmPeakFile;
        this.tauInfimum = tauInfimum;
        this.tauSupremum = tauSupremum;
        this.ipSNPReadInfimum = ipSNPReadInfimum;
        this.wesSNPReadInfimum = wesSNPReadInfimum;
        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
        this.log = initLog(new File(asmPeakFile).getParent());
    }

    public static void main(String[] args) throws ParseException {
        Options options = new Options();
        CommandLine commandLine = setCommandLine(args, options);

        String bedFile = null, aseVcfFile = null, wesVcfFile = null, outputFile, outputDir,
               peakCoveredSnpFile, peakCoveredSnpBackgroundFile;
        int ipSNPCoverageInfimum = 0, wesSNPCoverageInfimum = 0, samplingTime = 5000, burn_in = 200;
        double infimum = 0, supremum = 2;

        if (!commandLine.hasOption("o"))
            outputFile = new File(System.getProperty("user.dir"), "asmPeak.txt").getAbsolutePath();
        else
            outputFile = commandLine.getOptionValue("o");

        outputDir = new File(outputFile).getParent();
        Logger logger = initLog(outputDir);

        if (!commandLine.hasOption("bed")) {
            logger.error("Peak calling BED format file can not be empty");
            System.exit(2);
        } else {
            File bed = new File(commandLine.getOptionValue("bed"));
            if (!bed.exists() || !bed.isFile()) {
                logger.error("invalid file path: " + bed.getAbsolutePath());
                System.exit(2);
            }
            bedFile = bed.getAbsolutePath();
        }

        if (!commandLine.hasOption("vcf")) {
            logger.error("ASE SNP calling VCF file can not be empty");
            System.exit(2);
        } else {
            File vcf = new File(commandLine.getOptionValue("vcf"));
            if (!vcf.exists() || !vcf.isFile()) {
                logger.error("invalid file path: " + vcf.getAbsolutePath());
                System.exit(2);
            }
            aseVcfFile = vcf.getAbsolutePath();
        }

        if (commandLine.hasOption("wes")) {
            File vcf = new File(commandLine.getOptionValue("wes"));
            if (!vcf.exists() || !vcf.isFile()) {
                logger.error("invalid file path: " + vcf.getAbsolutePath());
                System.exit(2);
            }
            wesVcfFile = vcf.getAbsolutePath();
        }

        if (commandLine.hasOption("peak_cover_snp"))
            peakCoveredSnpFile = commandLine.getOptionValue("peak_cover_snp");
        else
            peakCoveredSnpFile = new File(outputDir, "PeakCoveredSNP.txt").getAbsolutePath();

        if (commandLine.hasOption("peak_cover_snp_bkg"))
            peakCoveredSnpBackgroundFile = commandLine.getOptionValue("peak_cover_snp_bkg");
        else
            peakCoveredSnpBackgroundFile = new File(outputDir, "PeakCoveredSNPBackground.txt").getAbsolutePath();

        if (commandLine.hasOption("s"))
            samplingTime = Integer.parseInt(commandLine.getOptionValue("s"));
        if (commandLine.hasOption("b"))
            burn_in = Integer.parseInt(commandLine.getOptionValue("b"));
        if (samplingTime <= 500 || burn_in <= 100) {
            logger.error("sampling times larger than 500 and burn in times at least 100");
            System.exit(2);
        }
        if (commandLine.hasOption("tl"))
            infimum = Double.parseDouble(commandLine.getOptionValue("tl"));
        if (commandLine.hasOption("th"))
            supremum = Double.parseDouble(commandLine.getOptionValue("th"));
        if (infimum >= supremum) {
            System.out.println("invalid uniform distribution parameter for tau sampling.");
            System.exit(2);
        }
        if (commandLine.hasOption("ip_cov_infimum"))
            ipSNPCoverageInfimum = Integer.parseInt(commandLine.getOptionValue("ip_cov_infimum"));
        if (commandLine.hasOption("wes_cov_infimum"))
            wesSNPCoverageInfimum = Integer.parseInt(commandLine.getOptionValue("wes_cov_infimum"));

        AsmPeakDetection apd = new AsmPeakDetection(bedFile, aseVcfFile, wesVcfFile, outputFile, peakCoveredSnpFile,
                                                    peakCoveredSnpBackgroundFile, infimum, supremum, ipSNPCoverageInfimum,
                                                    wesSNPCoverageInfimum, samplingTime, burn_in);
        apd.getTestResult();
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
     * 通过INPUT样本的VCF结果和IP样本的BED结果，得到被peak覆盖的SNP位点，并将记录写入文件。文件记录内容
     * #chr	SNP strand	position	peakStart	peakEnd	majorAlleleStrand	majorCount	minorCount
     */
    private void getPeakCoveredSnpResult() {
        PeakCoveredSNPRecord pcsr = new PeakCoveredSNPRecord(this.vcfFile, this.peakBedFile,
                                                             this.peakCoveredSnpFile, this.ipSNPReadInfimum);
        pcsr.getPeakCoveredSNP();
    }

    /**
     * 当存在WES数据SNP calling结果时，得到被peak覆盖的SNP位点，将结果写入文件作为背景。文件记录内容
     * #chr	SNP strand	position	peakStart	peakEnd	majorAlleleStrand	majorCount	minorCount
     */
    private void getPeakCoveredSnpBackground() {
        if (this.wesFile != null) {
            PeakCoveredSNPRecord pcsr = new PeakCoveredSNPRecord(this.wesFile, this.peakBedFile,
                                                                 this.peakCoveredWesSnpFile, this.wesSNPReadInfimum);
            pcsr.getPeakCoveredSNP();
        }
    }

    /**
     * 获取RNA-seq得到的每个m6A信号范围内SNV位点上的 major allele和 minor allele的reads count
     * @return [chr1: [peak1: position1: [major: count, minor: count], position2:[major: count, minor:count]], chr2:....]
     */
    private HashMap<String, HashMap<String, HashMap<String, HashMap<String, Integer>>>> getPeakSNPReadsCount() {
        HeterozygoteReadsCount hrc = new HeterozygoteReadsCount(this.peakCoveredSnpFile, this.log);
        this.peakMajorAlleleStrand = hrc.getPeakMajorAlleleStrand();
        return hrc.getMajorMinorHaplotype();
    }

    /**
     * 获取WES得到的每个m6A信号范围内SNV位点上的 major allele和 minor allele的reads count
     * @return [chr1: [peak1: position1: [major: count, minor: count], position2:[major: count, minor:count]], chr2:....]
     */
    private HashMap<String, HashMap<String, HashMap<String, HashMap<String, Integer>>>> getPeakSNPBackground() {
        HeterozygoteReadsCount hrc = new HeterozygoteReadsCount(this.peakCoveredWesSnpFile, this.log);
        return hrc.getMajorMinorHaplotype();
    }

    /**
     * 对每个覆盖SNP的peak进行检验
     */
    private void asmPeakTest() {
        // [chr1: [peak1: position1: [major: count, minor: count], position2:[major: count, minor:count]], chr2:....]
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
        double pVal, maf;
        HashMap<String, HashMap<String, HashMap<String, Integer>>> rnaSeqPeakSnvAlleleReads, wesPeakSnvAlleleReads;
        HashMap<String, HashMap<String, Integer>> rnaSeqMutPositionAlleleReads, wesMutPositionAlleleReads;
        HashMap<String, Integer> rnaSeqReads, wesReads;
        int[] majorCount, minorCount, majorBackground, minorBackground;
        int major, minor, backgroundMajor, backgroundMinor;
        String majorAllele, minorAllele, label;
        ArrayList<Integer> rnaSeqMajor, rnaSeqMinor, wesMajor, wesMinor;

        for (String chrNum: this.peakSnpReadsCount.keySet()) {
            // 某条染色体上所有Peak记录
            rnaSeqPeakSnvAlleleReads = this.peakSnpReadsCount.get(chrNum);

            if (this.peakSnpBackground != null)
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

                    // 依据reads count从大到小排序
                    List<Map.Entry<String, Integer>> nucleotides = new ArrayList<>(rnaSeqReads.entrySet());
                    Collections.sort(nucleotides, new Comparator<Map.Entry<String, Integer>>() {
                        @Override
                        public int compare(Map.Entry<String, Integer> o1, Map.Entry<String, Integer> o2) {
                            return o1.getValue().compareTo(o2.getValue());
                        }
                    });

                    ArrayList<String> ncs = new ArrayList<>();
                    ArrayList<Integer> counts = new ArrayList<>();
                    for (Map.Entry<String, Integer> entry: nucleotides) {
                        String nc = entry.getKey();
                        Integer reads = entry.getValue();
                        ncs.add(nc);
                        counts.add(reads);
                    }

                    major = counts.get(0);
                    minor = counts.get(1);
                    majorAllele = ncs.get(0);
                    minorAllele = ncs.get(1);
                    rnaSeqMajor.add(major);
                    rnaSeqMinor.add(minor);
                    if (wesReads != null && wesReads.keySet().contains(majorAllele) && wesReads.keySet().contains(minorAllele)) {
                        backgroundMajor = wesReads.get(majorAllele);
                        backgroundMinor = wesReads.get(minorAllele);
                    } else {
                        backgroundMajor = (minor+major)/2;
                        backgroundMinor = (minor+major)/2;
                    }
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

                label = String.join(":", new String[]{chrNum, peakRange});
                maf = this.calcMajorAlleleFrequency(majorCount, minorCount);
                this.peakMajorAlleleFrequency.put(label, maf);
                this.peakSNVNum.put(label, majorCount.length);

                Integer totalMajor = this.getSum(majorCount), totalMinor = this.getSum(minorCount);
                HashMap<String, Integer> record = new HashMap<>();
                record.put("major", totalMajor);
                record.put("minor", totalMinor);
                this.peakMajorMinorAlleleCount.put(label, record);

                Integer totalMajorBkg, totalMinorBkg;
                if (wesMutPositionAlleleReads == null) {
                    totalMajorBkg = 0;
                    totalMinorBkg = 0;
                } else {
                    totalMajorBkg = this.getSum(majorBackground);
                    totalMinorBkg = this.getSum(minorBackground);
                }
                HashMap<String, Integer> background = new HashMap<>();
                background.put("major", totalMajorBkg);
                background.put("minor", totalMinorBkg);
                this.peakMajorMinorBackground.put(label, background);

                rnaSeqMajor.clear();
                rnaSeqMinor.clear();
                wesMajor.clear();
                wesMinor.clear();

                // 对peak下所有的SNP位点进行元分析, 得到该peak对应的 p value
                HierarchicalBayesianModel hb = new HierarchicalBayesianModel(this.tauInfimum, this.tauSupremum, this.samplingTime,
                        this.burnIn, majorCount, minorCount, majorBackground, minorBackground);
                pVal = hb.testSignificant();
                ArrayList<String> samePValPeaks = this.asmPValue.getOrDefault(pVal, new ArrayList<>());
                samePValPeaks.add(label);
                this.asmPValue.put(pVal, samePValPeaks);
                hb = null;
            }
        }
    }

    private double calcMajorAlleleFrequency(int[] majorCounts, int[] minorCounts) {
        int major = 0, minor = 0;
        for (int i=0; i<majorCounts.length; i++) {
            major += majorCounts[i];
            minor += minorCounts[i];
        }

        return (double) major / (double) (major + minor);
    }

    /**
     * 对Peak ASM的p值进行BH校正, 得到 Q value.
     */
    private void bhRecalibrationOfEachPeak() {
        ArrayList<Map.Entry<Double, ArrayList<String>>> sortedPVals = new ArrayList<>(this.asmPValue.entrySet());
        // p值从小到大排序
        Collections.sort(sortedPVals, new Comparator<Map.Entry<Double, ArrayList<String>>>() {
            @Override
            public int compare(Map.Entry<Double, ArrayList<String>> o1, Map.Entry<Double, ArrayList<String>> o2) {
                return o2.getKey().compareTo(o1.getKey());
            }
        });


        int totalPeak = sortedPVals.size(), rankage = totalPeak;
        double qValue;
        String pValString, qValString;
        for (Map.Entry<Double, ArrayList<String>> entry: sortedPVals) {
            Double pVal = entry.getKey();
            ArrayList<String> samePValPeaks = entry.getValue();

            // 相同p值的Peak上SNV的数目
            HashMap<String, Integer> samePValPeaksSNVs = new HashMap<>();
            for (String peak: samePValPeaks)
                samePValPeaksSNVs.put(peak, this.peakSNVNum.get(peak));
            // 相同p值的Peak的major allele frequency
            HashMap<String, Double> samePValPeakMajorAlleleFrequency = new HashMap<>();
            for (String peak: samePValPeaks)
                samePValPeakMajorAlleleFrequency.put(peak, this.peakMajorAlleleFrequency.get(peak));

            List<Map.Entry<String, Integer>> samePValPeakEntry = new ArrayList<>(samePValPeaksSNVs.entrySet());
            Collections.sort(samePValPeakEntry, new Comparator<Map.Entry<String, Integer>>() {
                @Override
                public int compare(Map.Entry<String, Integer> o1, Map.Entry<String, Integer> o2) {

                    // 首先按照MAF进行排序，若MAF相同，则按照SNV的数目进行排序
                    String peak1 = o1.getKey(), peak2 = o2.getKey();
                    Double peak1MAF = samePValPeakMajorAlleleFrequency.get(peak1), peak2MAF = samePValPeakMajorAlleleFrequency.get(peak2);
                    if (Math.abs(peak1MAF - peak2MAF) < 0.00001) {
                        Integer peak1SNVs = o1.getValue(), peak2SNVs = o2.getValue();
                        return peak2SNVs.compareTo(peak1SNVs);
                    } else
                        return peak2MAF.compareTo(peak1MAF);
                }
            });

            for (Map.Entry<String, Integer> geneEntry: samePValPeakEntry) {
                String peak = geneEntry.getKey();
                qValue = Math.min(1.0, pVal * totalPeak / rankage);
//                System.out.println(geneName + "\t" + geneEntry.getValue() + "\t" + pVal + "\t" + qValue + samePValGeneMajorAlleleFrequency.get(geneEntry.getKey()) + "\t" + rankage);
                rankage--;

                pValString = Double.toString(pVal);
                qValString = this.df.format(qValue);
                this.asmQValue.add(String.join("->", new String[]{peak, pValString, qValString}));
            }
        }
    }

    /**
     * 将ASM的检验结果输出到文件
     */
    private void outputResult() {
        ArrayList<String[]> outputRecord = new ArrayList<>();
        BufferedWriter bfw = null;
        try {
            bfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(this.asmPeakFile))));
            String line, label, geneId, chrNum, peakStart, peakEnd, majorAlleleReads, minorAlleleReads,
                   majorAlleleBackground, minorAlleleBackground, majorAlleleStrand, pValue, qValue;
            String[] info, rec, finalInfo;
            int snvNum;
            double majorAlleleFrequency;
            bfw.write("#chr\tpeakStart\tpeakEnd\tgeneId\tp-value\tq-value\tsnvNum\tmajor/minorAlleleReads\tmajor/minorBackground\tMajorAlleleFrequency\tmajorAlleleStrand\n");
            for (String record: this.asmQValue) {
                rec = record.split("->");
                label = rec[0];
                // info = [chrNum, peakStart, peakEnd]
                info = label.split(":");
                chrNum = info[0];
                peakStart = info[1];
                peakEnd = info[2];
                pValue = rec[1];
                qValue = rec[2];
                geneId = this.peakCoveredGene.get(label);
                majorAlleleFrequency = this.peakMajorAlleleFrequency.get(label);
                majorAlleleStrand = this.peakMajorAlleleStrand.get(label);
                snvNum = this.peakSNVNum.get(label);
                majorAlleleReads = this.peakMajorMinorAlleleCount.get(label).get("major").toString();
                minorAlleleReads = this.peakMajorMinorAlleleCount.get(label).get("minor").toString();
                majorAlleleBackground = this.peakMajorMinorBackground.get(label).get("major").toString();
                minorAlleleBackground = this.peakMajorMinorBackground.get(label).get("minor").toString();
                finalInfo = new String[]{chrNum, peakStart, peakEnd, geneId, pValue, qValue,
                                   Integer.toString(snvNum), majorAlleleReads + "," + minorAlleleReads,
                                   majorAlleleBackground+","+minorAlleleBackground, Double.toString(majorAlleleFrequency),
                                   majorAlleleStrand};
                outputRecord.add(finalInfo);
            }

            Collections.sort(outputRecord, new Comparator<String[]>() {
                @Override
                public int compare(String[] o1, String[] o2) {
                    Double q1 = Double.parseDouble(o1[5]), q2 = Double.parseDouble(o2[5]);
                    if (!q1.equals(q2))
                        return q1.compareTo(q2);
                    // BH校正后q-value有相同值，此时依据SNV的数目进行降序排列
                    Integer snvCount1 = Integer.parseInt(o1[6]), snvCount2 = Integer.parseInt(o2[6]);
                    if (snvCount1 - snvCount2 != 0)
                        return snvCount2.compareTo(snvCount1);
                    // 若SNV数目相同，则依据major reads和minor reads的差异大小进行排序，差异大的靠前

                    Double maf1 = Double.parseDouble(o1[9]), maf2 = Double.parseDouble(o2[9]);
                    return maf2.compareTo(maf1);
                }
            });
            for (String[] record: outputRecord) {
                line = String.join("\t", record);
                bfw.write(line);
                bfw.newLine();
            }
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

    private Integer getSum(int[] readsCount) {
        int total = 0;
        for (int i: readsCount) {
            total += i;
        }

        return total;
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

    private static CommandLine setCommandLine(String[] args, Options options) throws ParseException {
        Option option = new Option("vcf", "vcf_file", true, "IP sample SNP calling result VCF file");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("wes", "wes_vcf_file", true, "WES SNP calling VCF format file");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("bed", "peak_bed_file", true, "Peak calling output result in BED format");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("o", "output", true, "ASE gene test output file, default ./asmPeak.txt");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("peak_cover_snp", "peak_cover_snp_record_file", true, "Peak covered INPUT sample SNP record");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("peak_cover_snp_bkg", "peak_cover_snp_background_record_file", true, "Peak covered WES sample SNP record");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("tl", "tau_low", true, "infimum of uniform distribution for sampling model hyper-parameter tau, default 0");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("th", "tau_high", true, "supremum of uniform distribution for sampling model hyper-parameter tau, default 2");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("ip_cov_infimum", "ip_snp_coverage_infimum", true, "IP sample SNP site coverage infimum, default 0");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("wes_cov_infimum", "wes_snp_coverage_infimum", true, "WES sample SNP site coverage infimum, default 0");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("s", "sampling", true, "sampling times, larger than 500, default 5000");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("b", "burn", true, "burn-in times, more than 100 and less than sampling times. Default 200");
        option.setRequired(false);
        options.addOption(option);

        CommandLineParser parser = new DefaultParser();

        return parser.parse(options, args);
    }
}
