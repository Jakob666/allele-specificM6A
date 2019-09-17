package HierarchicalBayesianAnalysis;

import GTFComponent.VcfSnpMatchGene;
import org.apache.commons.cli.*;
import org.apache.log4j.Logger;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

/**
 * 通过SNP calling得到的VCF文件对ASE gene进行检验
 */
public class AseGeneDetection {
    private String aseGeneFile;
    // geneAlleleReads = {"geneId->geneName": {pos1: [refAllele:count, altAllele: count1]}, ...}
    private HashMap<String, HashMap<Integer, String[]>> geneAlleleReads, geneBackgroundReads;
    private HashMap<String, HashMap<Integer, String>> geneMajorStrand;
    private HashMap<String, HashMap<String, Integer>> geneReadsCount = new HashMap<>(), geneBackgroundCount = new HashMap<>();
    private int samplingTime, burnIn;
    private HashMap<Double, ArrayList<String>> geneAsePValue = new HashMap<>();
    private HashMap<String, Integer> geneSNVs = new HashMap<>();
    private HashMap<String, Double> geneMajorAlleleFrequency = new HashMap<>();
    private HashMap<String, String> geneMajorNucleotide = new HashMap<>();
    private ArrayList<String> geneAseQValue;
    private DecimalFormat df = new DecimalFormat("0.0000");
    private double infimum, supremum;

    /**
     * Constructor
     * @param gtfFile GTF文件
     * @param vcfFile INPUT样本 SNP calling得到的VCF文件
     * @param wesFile WES样本 SNP calling得到的VCF文件
     * @param aseGeneFile 结果输出文件
     * @param infimum tau采样时的下确界
     * @param supremum tau采样时的上确界
     * @param readsCoverageThreshold SNV筛选时设置的reads coverage阈值
     * @param samplingTime 采样次数
     * @param burnIn burn in次数
     */
    public AseGeneDetection(String gtfFile, String vcfFile, String wesFile, String aseGeneFile,
                            double infimum, double supremum, int readsCoverageThreshold, int samplingTime, int burnIn) {
        VcfSnpMatchGene vsmg = new VcfSnpMatchGene(vcfFile, gtfFile, readsCoverageThreshold);
        vsmg.parseVcfFile();
        this.geneAlleleReads = vsmg.getGeneAlleleReads();
        this.geneMajorStrand = vsmg.getGeneMajorAlleleNucleotide();
        if (wesFile != null) {
            vsmg = new VcfSnpMatchGene(wesFile, gtfFile, readsCoverageThreshold);
            vsmg.parseVcfFile();
            this.geneBackgroundReads = vsmg.getGeneAlleleReads();
        } else {
            this.geneBackgroundReads = null;
        }

        this.infimum = infimum;
        this.supremum = supremum;
        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
        this.aseGeneFile = aseGeneFile;
    }

    public static void main(String[] args) throws ParseException {
        Options options = new Options();
        CommandLine commandLine = setCommandLine(args, options);
        String gtfFile = null, aseVcfFile = null, wesVcfFile = null, outputFile, outputDir;
        int samplingTime = 5000, burn_in = 200, readsCoverageThreshold = 10;
        double infimum = 0, supremum = 2;
        if (!commandLine.hasOption("o"))
            outputFile = new File(System.getProperty("user.dir"), "aseGene.txt").getAbsolutePath();
        else
            outputFile = commandLine.getOptionValue("o");
        outputDir = new File(outputFile).getParent();
        Logger logger = initLog(outputDir);

        if (!commandLine.hasOption("g")) {
            logger.error("GTF annotation file can not be empty");
            System.exit(2);
        } else {
            File gtf = new File(commandLine.getOptionValue("g"));
            if (!gtf.exists() || !gtf.isFile()) {
                logger.error("invalid file path: " + gtf.getAbsolutePath());
                System.exit(2);
            }
            gtfFile = gtf.getAbsolutePath();
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
        if (commandLine.hasOption("rc"))
            readsCoverageThreshold = Integer.valueOf(commandLine.getOptionValue("rc"));
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

        AseGeneDetection agd = new AseGeneDetection(gtfFile, aseVcfFile, wesVcfFile, outputFile,
                                                    infimum, supremum, readsCoverageThreshold, samplingTime, burn_in);
        agd.getTestResult();
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
        HashMap<Integer, String[]> geneSNVs, wesSNVs;
        ArrayList<Integer> majorAlleleCount, minorAlleleCount, majorBackgroundCount, minorBackgroundCount;
        int[] majorCount, minorCount, majorBackground, minorBackground;
        HashMap<Integer, String> refOrAlt;
        HierarchicalBayesianModel hb;
        double p;
        String geneId;
        for (String label: this.geneAlleleReads.keySet()) {
            geneId = label.split("->")[0];
            refOrAlt = this.geneMajorStrand.get(geneId);
            int majorReadsCount = 0, minorReadsCount = 0;
            majorAlleleCount = new ArrayList<>();
            minorAlleleCount = new ArrayList<>();
            majorBackgroundCount = new ArrayList<>();
            minorBackgroundCount = new ArrayList<>();

            // 获取RNA-seq和WES经过SNP calling流程得到的某个基因的reference allele和alternative allele 的类型和reads数目
            geneSNVs = this.geneAlleleReads.get(label);
            if (this.geneBackgroundReads != null)
                wesSNVs = this.geneBackgroundReads.getOrDefault(label, null);
            else
                wesSNVs = null;

            // 记录major allele的位置和对应的碱基
            LinkedList<String> majorNcRecords = new LinkedList<>();
            for (Integer mutPosition: geneSNVs.keySet()) {
                String[] nucleotideReadsCount = geneSNVs.get(mutPosition);
                String majorAlleleRecord, minorAlleleRecord, majorNC, minorNC;
                // 获取major allele和minor allele及reads
                if ((refOrAlt.get(mutPosition)).equals("ref")) {
                    majorAlleleRecord = nucleotideReadsCount[0];
                    minorAlleleRecord = nucleotideReadsCount[1];
                } else {
                    minorAlleleRecord = nucleotideReadsCount[0];
                    majorAlleleRecord = nucleotideReadsCount[1];
                }
                majorNC = majorAlleleRecord.split(":")[0];
                minorNC = minorAlleleRecord.split(":")[0];
                int major = Integer.parseInt(majorAlleleRecord.split(":")[1]);
                int minor = Integer.parseInt(minorAlleleRecord.split(":")[1]);
                majorAlleleCount.add(major);
                minorAlleleCount.add(minor);
                majorNcRecords.add(String.join(":", new String[] {Integer.toString(mutPosition), majorNC}));

                // WES数据是否存在相应碱基reads count的background
                if (wesSNVs == null) {
                    int alleleCount = (major + minor) / 2;
                    majorBackgroundCount.add(alleleCount);
                    minorBackgroundCount.add(alleleCount);
                } else {
                    String[] snvSite = wesSNVs.getOrDefault(mutPosition, null);
                    int alleleCount = (major + minor) / 2;
                    if (snvSite == null) {
                        majorBackgroundCount.add(alleleCount);
                        minorBackgroundCount.add(alleleCount);
                    } else {
                        int[] wesReads = this.getWesReads(snvSite, majorNC, minorNC);
                        if (wesReads != null) {
                            majorBackgroundCount.add(wesReads[0]);
                            minorBackgroundCount.add(wesReads[1]);
                        } else {
                            majorBackgroundCount.add(alleleCount);
                            minorBackgroundCount.add(alleleCount);
                        }

                    }
                }
            }
            this.geneMajorNucleotide.put(geneId, this.getString(majorNcRecords));

            HashMap<String, Integer> readsCount = new HashMap<>();
            readsCount.put("major", this.getSum(majorAlleleCount));
            readsCount.put("minor", this.getSum(minorAlleleCount));
            this.geneReadsCount.put(label, readsCount);
            HashMap<String, Integer> backgroundCount = new HashMap<>();
            backgroundCount.put("major", this.getSum(majorBackgroundCount));
            backgroundCount.put("minor", this.getSum(minorBackgroundCount));
            this.geneBackgroundCount.put(label, backgroundCount);

            assert majorAlleleCount.size() == minorAlleleCount.size();
            assert majorAlleleCount.size() == majorBackgroundCount.size();
            majorCount = new int[majorAlleleCount.size()];
            minorCount = new int[minorAlleleCount.size()];
            majorBackground = new int[majorBackgroundCount.size()];
            minorBackground = new int[minorBackgroundCount.size()];
            for (int i=0; i<majorAlleleCount.size(); i++) {
                majorCount[i] = majorAlleleCount.get(i);
                majorReadsCount += majorAlleleCount.get(i);
                minorCount[i] = minorAlleleCount.get(i);
                minorReadsCount += minorAlleleCount.get(i);

                majorBackground[i] = majorBackgroundCount.get(i);
                minorBackground[i] = minorBackgroundCount.get(i);
            }

            this.geneMajorAlleleFrequency.put(label, (double) majorReadsCount / (double) (majorReadsCount + minorReadsCount));
            hb = new HierarchicalBayesianModel(this.infimum, supremum, this.samplingTime, this.burnIn,
                                                majorCount, minorCount, majorBackground, minorBackground);
            p = hb.testSignificant();

            ArrayList<String> samePValGenes = this.geneAsePValue.getOrDefault(p, new ArrayList<>());
            samePValGenes.add(label);
            this.geneAsePValue.put(p, samePValGenes);
            this.geneSNVs.put(label, majorCount.length);
        }
    }

    private int[] getWesReads(String[] wesRecord, String majorNC, String minorNC) {
        int[] res = null;
        HashMap<String, Integer> wesNC = new HashMap<>();
        for (String nc: wesRecord) {
            wesNC.put(nc.split(":")[0], Integer.parseInt(nc.split(":")[1]));
        }
        if (wesNC.keySet().contains(majorNC) && wesNC.keySet().contains(minorNC))
            res = new int[] {wesNC.get(majorNC), wesNC.get(minorNC)};
        wesNC.clear();

        return res;
    }

    /**
     * 对基因ASE的p值进行BH校正, 得到 Q value.
     */
    private void bhRecalibrationOfEachPeak() {
        ArrayList<Map.Entry<Double, ArrayList<String>>> sortedPVals = new ArrayList<>(this.geneAsePValue.entrySet());
        // p值从小到大排序
        Collections.sort(sortedPVals, new Comparator<Map.Entry<Double, ArrayList<String>>>() {
            @Override
            public int compare(Map.Entry<Double, ArrayList<String>> o1, Map.Entry<Double, ArrayList<String>> o2) {
                return o2.getKey().compareTo(o1.getKey());
            }
        });

        int totalGenes = this.geneSNVs.size(), rankage = totalGenes;
        double qValue;
        String pValString, qValString;
        this.geneAseQValue = new ArrayList<>(totalGenes);

        // 对相同p值的基因进行排序
        for (Map.Entry<Double, ArrayList<String>> entry: sortedPVals) {
            Double pVal = entry.getKey();
            ArrayList<String> samePValGenes = entry.getValue();

            // 相同p值的基因上SNV的数目
            HashMap<String, Integer> samePValGeneSNVs = new HashMap<>();
            for (String gene: samePValGenes)
                samePValGeneSNVs.put(gene, this.geneSNVs.get(gene));
            // 相同p值的基因的major allele frequency
            HashMap<String, Double> samePValGeneMajorAlleleFrequency = new HashMap<>();
            for (String gene: samePValGenes)
                samePValGeneMajorAlleleFrequency.put(gene, this.geneMajorAlleleFrequency.get(gene));

            List<Map.Entry<String, Integer>> samePValGeneEntry = new ArrayList<>(samePValGeneSNVs.entrySet());
            Collections.sort(samePValGeneEntry, new Comparator<Map.Entry<String, Integer>>() {
                @Override
                public int compare(Map.Entry<String, Integer> o1, Map.Entry<String, Integer> o2) {

                    // 首先按照MAF进行排序，若MAF相同，则按照SNV的数目进行排序
                    String gene1 = o1.getKey(), gene2 = o2.getKey();
                    Double gene1MAF = samePValGeneMajorAlleleFrequency.get(gene1), gene2MAF = samePValGeneMajorAlleleFrequency.get(gene2);
                    if (Math.abs(gene1MAF - gene2MAF) < 0.00001) {
                        Integer gene1SNVs = o1.getValue(), gene2SNVs = o2.getValue();
                        return gene2SNVs.compareTo(gene1SNVs);
                    } else
                        return gene2MAF.compareTo(gene1MAF);
                }
            });

            for (Map.Entry<String, Integer> geneEntry: samePValGeneEntry) {
                String geneName = geneEntry.getKey();
                qValue = Math.min(1.0, pVal * totalGenes / rankage);
//                System.out.println(geneName + "\t" + geneEntry.getValue() + "\t" + pVal + "\t" + qValue + samePValGeneMajorAlleleFrequency.get(geneEntry.getKey()) + "\t" + rankage);
                rankage--;

                pValString = Double.toString(pVal);
                qValString = this.df.format(qValue);
                this.geneAseQValue.add(String.join("->", new String[]{geneName, pValString, qValString}));
            }
        }
    }

    /**
     * 将检验结果写入文件
     */
    private void outputResult() {
        HashMap<String, String[]> finalRecords = new HashMap<>();
        BufferedWriter bfw = null;
        try {
            bfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(this.aseGeneFile))));
            String line, label, geneId, geneName, majorAlleleStrand, pVal, qVal;
            String[] info;
            int snvNum, majorCount, minorCount, majorBackground, minorBackground;
            double majorAlleleFrequency;
            bfw.write("#geneId\tgeneName\tp-value\tq-value\tsnvNum\tmajor/minorAlleleReads\tmajor/minorBackground\tMajorAlleleFrequency\tmajorAlleleNC\n");
            for (String record: this.geneAseQValue) {
                info = record.split("->");
                geneId = info[0];
                geneName = info[1];
                pVal = info[2];
                qVal = info[3];
                label = String.join("->", new String[]{geneId, geneName});
                majorCount = this.geneReadsCount.get(label).get("major");
                minorCount = this.geneReadsCount.get(label).get("minor");
                majorAlleleFrequency = (double) majorCount / (double) (majorCount + minorCount);
                majorAlleleStrand = this.geneMajorNucleotide.get(geneId);
                snvNum = this.geneSNVs.get(label);
                majorBackground = (this.geneBackgroundReads == null)? 0: this.geneBackgroundCount.get(label).get("major");
                minorBackground = (this.geneBackgroundReads == null)? 0: this.geneBackgroundCount.get(label).get("minor");

                finalRecords.put(geneName, new String[]{geneId, geneName, pVal, qVal, Integer.toString(snvNum),
                        Integer.toString(majorCount), Integer.toString(minorCount), Integer.toString(majorBackground),
                        Integer.toString(minorBackground), this.df.format(majorAlleleFrequency), majorAlleleStrand});
            }

            List<Map.Entry<String, String[]>> records = new ArrayList<>(finalRecords.entrySet());
            // 按照q值进行排序
            Collections.sort(records, new Comparator<Map.Entry<String, String[]>>() {
                @Override
                public int compare(Map.Entry<String, String[]> o1, Map.Entry<String, String[]> o2) {
                    String[] data1 = o1.getValue(), data2 = o2.getValue();
                    Double q1 = Double.parseDouble(data1[3]), q2 = Double.parseDouble(data2[3]);
                    if (!q1.equals(q2))
                        return q1.compareTo(q2);
                    // BH校正后q-value有相同值，此时依据SNV的数目进行降序排列
                    Integer snvCount1 = Integer.parseInt(data1[4]), snvCount2 = Integer.parseInt(data2[4]);
                    if (snvCount1 - snvCount2 != 0)
                        return snvCount2.compareTo(snvCount1);
                    // 若SNV数目相同，则依据major reads和minor reads的差异大小进行排序，差异大的靠前
                    Integer majorCount1 = Integer.parseInt(data1[5]), minorCount1 = Integer.parseInt(data1[6]),
                            majorCount2 = Integer.parseInt(data2[5]), minorCount2 = Integer.parseInt(data2[6]);
                    Double maf1 = (double) majorCount1 / (double)(majorCount1 + minorCount1),
                            maf2 = (double) majorCount2 / (double) (majorCount2 + minorCount2);
                    return maf2.compareTo(maf1);
                }
            });

            for (Map.Entry<String, String[]> rec: records) {
                String[] data = rec.getValue();
                String[] lineInfo = new String[]{data[0], data[1], data[2], data[3], data[4], data[5] + "," + data[6], data[7]+ "," + data[8], data[9], data[10]};
                line = String.join("\t", lineInfo);
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

    private Integer getSum(ArrayList<Integer> list) {
        Integer total = 0;
        for (Integer i: list)
            total += i;

        return total;
    }

    private String getString(LinkedList<String> list) {
        Collections.sort(list, new Comparator<String>() {
            @Override
            public int compare(String o1, String o2) {
                Integer pos1 = Integer.parseInt(o1.split(":")[0]);
                Integer pos2 = Integer.parseInt(o2.split(":")[0]);
                return pos2.compareTo(pos1);
            }
        });

        String[] str = new String[list.size()];
        for (int i=0; i<list.size(); i++) {
            str[i] = list.get(i);
        }
        String res = String.join(";", str);
        str = null;

        return res;
    }

    private static Logger initLog(String logHome) {
        System.setProperty("log_home", logHome);
        return Logger.getLogger(AseGeneDetection.class);
    }

    private static CommandLine setCommandLine(String[] args, Options options) throws ParseException {
        Option option = new Option("vcf", "vcf_file", true, "INPUT sample SNP calling result VCF file");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("wes", "wes_vcf_file", true, "WES SNP calling VCF format file");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("g", "gtf", true, "GTF annotation file");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("o", "output", true, "ASE gene test output file, default ./aseGene.txt");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("tl", "tau_low", true, "infimum of uniform distribution for sampling model hyper-parameter tau, default 0");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("th", "tau_high", true, "supremum of uniform distribution for sampling model hyper-parameter tau, default 2");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("rc", "readsCoverage", true, "reads coverage threshold using for filter SNV records, default 10");
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
