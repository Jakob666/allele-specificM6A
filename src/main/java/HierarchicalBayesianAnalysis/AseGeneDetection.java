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
    // geneAlleleReads = {"geneId->geneName": {"major":[count1, count2,...], "minor": [count1, count2,...]}, ...}
    private String aseGeneFile;
    private HashMap<String, HashMap<String, ArrayList<Integer>>> geneAlleleReads;
    private HashMap<String, int[]> geneMajorStrand;
    private int samplingTime, burnIn;
    private HashMap<Double, ArrayList<String>> geneAsePValue = new HashMap<>();
    private HashMap<String, Integer> geneSNVs = new HashMap<>();
    private HashMap<String, Double> geneMajorAlleleFrequency = new HashMap<>();
    private ArrayList<String> geneAseQValue;
    private DecimalFormat df = new DecimalFormat("0.000000");

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

    public static void main(String[] args) throws ParseException {
        Options options = new Options();
        CommandLine commandLine = setCommandLine(args, options);

        String gtfFile = null, aseVcfFile = null, outputFile, outputDir;
        int samplingTime = 5000, burn_in = 200;
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

        if (!commandLine.hasOption("ase")) {
            logger.error("ASE SNP calling VCF file can not be empty");
            System.exit(2);
        } else {
            File vcf = new File(commandLine.getOptionValue("ase"));
            if (!vcf.exists() || !vcf.isFile()) {
                logger.error("invalid file path: " + vcf.getAbsolutePath());
                System.exit(2);
            }
            aseVcfFile = vcf.getAbsolutePath();
        }

        if (commandLine.hasOption("s"))
            samplingTime = Integer.parseInt(commandLine.getOptionValue("s"));
        if (commandLine.hasOption("b"))
            burn_in = Integer.parseInt(commandLine.getOptionValue("b"));
        if (samplingTime <= 500 || burn_in <= 100) {
            logger.error("sampling times larger than 500 and burn in times at least 100");
            System.exit(2);
        }

        AseGeneDetection agd = new AseGeneDetection(gtfFile, aseVcfFile, outputFile, samplingTime, burn_in);
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
        ArrayList<Integer> majorAlleleCount, minorAlleleCount;
        int[] majorCount, minorCount;
        HierarchicalBayesianModel hb;
        double p;
        for (String label: this.geneAlleleReads.keySet()) {
            int majorReadsCount = 0, minorReadsCount = 0;
            majorAlleleCount = this.geneAlleleReads.get(label).get("major");
            minorAlleleCount = this.geneAlleleReads.get(label).get("minor");

            assert majorAlleleCount.size() == minorAlleleCount.size();
            majorCount = new int[majorAlleleCount.size()];
            minorCount = new int[minorAlleleCount.size()];
            for (int i=0; i<majorAlleleCount.size(); i++) {
                majorCount[i] = majorAlleleCount.get(i);
                majorReadsCount += majorAlleleCount.get(i);
                minorCount[i] = minorAlleleCount.get(i);
                minorReadsCount += minorAlleleCount.get(i);
            }
            if (majorCount.length == 1 && (minorCount[0] ==0 && majorCount[0] < 10))  // majorCount[0] - minorCount[0] <= 4
                continue;
            this.geneMajorAlleleFrequency.put(label, (double) majorReadsCount / (double) (majorReadsCount + minorReadsCount));
            hb = new HierarchicalBayesianModel(0, 1, this.samplingTime, this.burnIn, majorCount, minorCount);
            p = hb.testSignificant();

            ArrayList<String> samePValGenes = this.geneAsePValue.getOrDefault(p, new ArrayList<>());
            samePValGenes.add(label);
            this.geneAsePValue.put(p, samePValGenes);
            this.geneSNVs.put(label, majorCount.length);
        }
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
        double qValue, prevQVal = sortedPVals.get(0).getKey();
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
                System.out.println(geneEntry.getKey() + "\t" + geneEntry.getValue() + "\t" + pVal + "\t" + qValue + samePValGeneMajorAlleleFrequency.get(geneEntry.getKey()) + "\t" + rankage);
                rankage--;

                pValString = Double.toString(pVal);
                qValString = Double.toString(qValue);
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
            int snvNum, majorCount, minorCount;
            int[] alleleStrand;
            double majorAlleleFrequency;
            ArrayList<Integer> majorAlleleCount, minorAlleleCount;
            bfw.write("#geneId\tgeneName\tp-value\tq-value\tsnvNum\tmajor/minorAlleleReads\tMajorAlleleFrequency\tmajorAlleleStrand\n");
            for (String record: this.geneAseQValue) {
                info = record.split("->");
                geneId = info[0];
                geneName = info[1];
                pVal = info[2];
                qVal = info[3];
                label = String.join("->", new String[]{geneId, geneName});
                majorAlleleCount = this.geneAlleleReads.get(label).get("major");
                minorAlleleCount = this.geneAlleleReads.get(label).get("minor");
                snvNum = majorAlleleCount.size();
                majorCount = this.getSum(majorAlleleCount);
                minorCount = this.getSum(minorAlleleCount);
                majorAlleleFrequency = (double) majorCount / (double) (majorCount + minorCount);
                alleleStrand = this.geneMajorStrand.get(geneId);
                if (alleleStrand[0] >= alleleStrand[1])
                    majorAlleleStrand = "+";
                else
                    majorAlleleStrand = "-";

                finalRecords.put(geneName, new String[]{geneId, geneName, pVal, qVal, Integer.toString(snvNum),
                        Integer.toString(majorCount), Integer.toString(minorCount),
                        Double.toString(majorAlleleFrequency), majorAlleleStrand});
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
                    Integer diff1 = majorCount1 - minorCount1, diff2 = majorCount2 - minorCount2;
                    return diff2.compareTo(diff1);
                }
            });

            for (Map.Entry<String, String[]> rec: records) {
                String[] data = rec.getValue();
                String[] lineInfo = new String[]{data[0], data[1], data[2], data[3], data[4], data[5] + "," + data[6], data[7], data[8]};
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

    private static Logger initLog(String logHome) {
        System.setProperty("log_home", logHome);
        return Logger.getLogger(AseGeneDetection.class);
    }

    private static CommandLine setCommandLine(String[] args, Options options) throws ParseException {
        Option option = new Option("ase", "ase_vcf_file", true, "INPUT sample SNP calling result VCF file");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("g", "gtf", true, "GTF annotation file");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("o", "output", true, "ASE gene test output file, default ./aseGene.txt");
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
