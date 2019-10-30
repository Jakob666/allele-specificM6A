package HierarchicalBayesianAnalysis;

import GTFComponent.VcfSnpMatchGene;
import heterozygoteSiteAnalysis.DbsnpAnnotation;
import org.apache.commons.cli.*;
import org.apache.log4j.Logger;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

/**
 * Test ASE genes with VCF format file
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
    private HashMap<String, ArrayList<int[]>> statisticForTest = new HashMap<>();
    private ArrayList<String> geneAseQValue;
    private DecimalFormat df = new DecimalFormat("0.0000");
    private double degreeOfFreedom;
    private Logger logger;

    /**
     * Constructor
     * @param gtfFile GTF annotation file
     * @param vcfFile VCF format file via MeRIP-seq INPUT data
     * @param wesFile VCF format file via WES data, optional
     * @param dbsnpFile dbsnp annotation file for SNP filtering
     * @param aseGeneFile test result output file
     * @param degreeOfFreedom the degree of freedom of inverse-Chi-square distribution, default 10
     * @param readsCoverageThreshold reads coverage threshold when filter INPUT sample SNV sites, default 10
     * @param wesCoverageThreshold reads coverage threshold when filter WES SNV sites, default 30
     * @param samplingTime sampling time, default 5000
     * @param burnIn burn in time, default 200
     */
    public AseGeneDetection(String gtfFile, String vcfFile, String wesFile, String dbsnpFile, String aseGeneFile,
                            double degreeOfFreedom, int readsCoverageThreshold, int wesCoverageThreshold,
                            int samplingTime, int burnIn, Logger logger) {
        HashMap<String, LinkedList<DbsnpAnnotation.DIYNode>> dbsnpRecord = null;
        if (dbsnpFile != null) {
            DbsnpAnnotation da = new DbsnpAnnotation(dbsnpFile, logger);
            da.parseDbsnpFile();
            dbsnpRecord = da.getDbsnpRecord();
        }

        VcfSnpMatchGene vsmg = new VcfSnpMatchGene(vcfFile, gtfFile, readsCoverageThreshold);
        vsmg.parseVcfFile(dbsnpRecord);
        this.geneAlleleReads = vsmg.getGeneAlleleReads();
        // {geneId: {pos: ref/alt, pos:ref/alt...}, ...}
        this.geneMajorStrand = vsmg.getGeneMajorAlleleNucleotide();
        if (wesFile != null) {
            vsmg = new VcfSnpMatchGene(wesFile, gtfFile, wesCoverageThreshold);
            vsmg.parseVcfFile(dbsnpRecord);
            this.geneBackgroundReads = vsmg.getGeneAlleleReads();
        } else {
            this.geneBackgroundReads = null;
        }
        if (dbsnpRecord != null)
            dbsnpRecord = null;

        this.degreeOfFreedom = degreeOfFreedom;
        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
        this.aseGeneFile = aseGeneFile;
        this.logger = logger;
    }

    public static void main(String[] args) throws ParseException {
        Options options = new Options();
        CommandLine commandLine = setCommandLine(args, options);
        String gtfFile = null, aseVcfFile = null, wesVcfFile = null, dbsnpFile = null, outputFile, outputDir;
        int samplingTime = 5000, burn_in = 200, readsCoverageThreshold = 10, wesCoverageThreshold = 30;
        double degreeOfFreedom = 10;
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
        if (commandLine.hasOption("db")) {
            File dbsnp = new File(commandLine.getOptionValue("db"));
            if (!dbsnp.exists() || !dbsnp.isFile()) {
                logger.error("invalid file path: " + dbsnp.getAbsolutePath());
                System.exit(2);
            }
            dbsnpFile = dbsnp.getAbsolutePath();
        }

        if (commandLine.hasOption("rc"))
            readsCoverageThreshold = Integer.valueOf(commandLine.getOptionValue("rc"));
        if (commandLine.hasOption("wc"))
            wesCoverageThreshold = Integer.valueOf(commandLine.getOptionValue("wc"));
        if (commandLine.hasOption("s"))
            samplingTime = Integer.parseInt(commandLine.getOptionValue("s"));
        if (commandLine.hasOption("b"))
            burn_in = Integer.parseInt(commandLine.getOptionValue("b"));
        if (samplingTime <= 500 || burn_in <= 100) {
            logger.error("sampling times larger than 500 and burn in times at least 100");
            System.exit(2);
        }
        if (commandLine.hasOption("df"))
            degreeOfFreedom = Double.parseDouble(commandLine.getOptionValue("df"));
        if (degreeOfFreedom <= 1) {
            System.out.println("invalid inverse-Chi-square distribution parameter for tau sampling. Must larger than 1.0");
            System.exit(2);
        }

        AseGeneDetection agd = new AseGeneDetection(gtfFile, aseVcfFile, wesVcfFile, dbsnpFile, outputFile,
                degreeOfFreedom, readsCoverageThreshold, wesCoverageThreshold, samplingTime, burn_in, logger);
        agd.getTestResult();
    }

    public void getTestResult() {
        this.aseGeneTest();
        this.bhRecalibrationOfEachPeak();
        this.outputResult();
    }

    /**
     * get ASE significant p value with hierarchical model
     */
    private void aseGeneTest() {
        HashMap<Integer, String[]> geneSNVs, wesSNVs;
        ArrayList<Integer> majorAlleleCount, minorAlleleCount, majorBackgroundCount, minorBackgroundCount;
        int[] majorCount, minorCount, majorBackground, minorBackground;
        HashMap<Integer, String> refOrAlt;
        HierarchicalBayesianModel hb;
        double p;
        String geneId, type;
        for (String label: this.geneAlleleReads.keySet()) {
            // WES data exists, but there is no WES SNV sites locate in the gene
            if (this.geneBackgroundReads != null) {
                if (this.geneBackgroundReads.keySet().contains(label))
                    wesSNVs = this.geneBackgroundReads.getOrDefault(label, null);
                else
                    continue;
            } else
                wesSNVs = null;

            geneId = label.split("->")[0];
            refOrAlt = this.geneMajorStrand.get(geneId);
            int majorReadsCount = 0, minorReadsCount = 0;
            majorAlleleCount = new ArrayList<>();
            minorAlleleCount = new ArrayList<>();
            majorBackgroundCount = new ArrayList<>();
            minorBackgroundCount = new ArrayList<>();

            // get reference allele and alternative allele nucleotide and reads count
            geneSNVs = this.geneAlleleReads.get(label);

            // record the genome positions of major alleles and the corresponding nucleotides
            LinkedList<String> majorNcRecords = new LinkedList<>();
            for (Integer mutPosition : geneSNVs.keySet()) {

                if (wesSNVs != null && !wesSNVs.keySet().contains(mutPosition))
                    continue;

                String[] nucleotideReadsCount = geneSNVs.get(mutPosition);
                String majorAlleleRecord, minorAlleleRecord, majorNC, minorNC;
                // get major allele and minor allele with their reads
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
                type = this.geneMajorStrand.get(geneId).get(mutPosition);
                majorNcRecords.add(String.join(":", new String[]{Integer.toString(mutPosition), type, majorNC}));

                int alleleCount = (major + minor) / 2;
                majorBackgroundCount.add(alleleCount);
                minorBackgroundCount.add(alleleCount);
            }
            if (majorAlleleCount.size() == 0) {
                majorAlleleCount = null;
                minorAlleleCount = null;
                majorBackgroundCount = null;
                minorBackgroundCount = null;
                continue;
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
            for (int i = 0; i < majorAlleleCount.size(); i++) {
                majorCount[i] = majorAlleleCount.get(i);
                majorReadsCount += majorAlleleCount.get(i);
                minorCount[i] = minorAlleleCount.get(i);
                minorReadsCount += minorAlleleCount.get(i);

                majorBackground[i] = majorBackgroundCount.get(i);
                minorBackground[i] = minorBackgroundCount.get(i);
            }

            double majorAlleleFrequency = (double) majorReadsCount / (double) (majorReadsCount + minorReadsCount);
            // MAF threshold for reducing false positive ASE genes. Genes with MAF < 0.55 are seen as normal mutated gene
            // cause by alignment error
            if (majorAlleleFrequency - 0.55 < 0.00001)
                continue;
            // MAF threshold for reducing false positive ASE genes. Genes with MAF > 0.95 are seen as homo-zygote gene
            // cause by alignment error
            if (majorAlleleFrequency - 0.95 > 0.00001)
                continue;
            this.geneMajorAlleleFrequency.put(label, majorAlleleFrequency);
            ArrayList<int[]> statistic = new ArrayList<>(4);
            statistic.add(majorCount);
            statistic.add(minorCount);
            statistic.add(majorBackground);
            statistic.add(minorBackground);
            this.statisticForTest.put(label, statistic);
        }

        double lorStd = this.calcLorStd();

        // test ASE gene with Hierarchical model
        for (String label: this.statisticForTest.keySet()) {
            ArrayList<int[]> statistic = this.statisticForTest.get(label);
            majorCount = statistic.get(0);
            minorCount = statistic.get(1);
            majorBackground = statistic.get(2);
            minorBackground = statistic.get(3);
            hb = new HierarchicalBayesianModel(lorStd, this.degreeOfFreedom, this.samplingTime, this.burnIn,
                                                majorCount, minorCount, majorBackground, minorBackground);
            p = hb.testSignificant();

            ArrayList<String> samePValGenes = this.geneAsePValue.getOrDefault(p, new ArrayList<>());
            samePValGenes.add(label);
            this.geneAsePValue.put(p, samePValGenes);
            this.geneSNVs.put(label, majorCount.length);
        }
        this.statisticForTest.clear();
    }

    /**
     * calculate the standard deviation of LOR of all SNV sites on genome
     * @return LOR Std
     */
    private double calcLorStd() {
        ArrayList<Double> lorList = new ArrayList<>();
        int[] majorCount, minorCount, majorBackground, minorBackground;
        double lor, cum = 0;
        for (String label: this.statisticForTest.keySet()) {
            ArrayList<int[]> statistic = this.statisticForTest.get(label);
            majorCount = statistic.get(0);
            minorCount = statistic.get(1);
            majorBackground = statistic.get(2);
            minorBackground = statistic.get(3);

            for (int i=0; i<majorCount.length; i++) {
                double major = majorCount[i], minor = minorCount[i],
                       majorBack = majorBackground[i], minorBack = minorBackground[i];
                if ((minor - 0) < 0.00001)
                    minor = 0.1;
                if ((minorBack - 0) < 0.00001) {
                    majorBack = (major + minor) / 2;
                    minorBack = (major + minor) / 2;
                }

                lor = (major / minor) / (majorBack / minorBack);
                lor = Math.log(lor);
                lorList.add(lor);
                cum += lor;
            }
        }
        double lorMean = cum / lorList.size();
        double variance = 0.0;
        for (Double val: lorList) {
            variance += Math.pow((val - lorMean), 2);
        }

        double lorStd = Math.sqrt(variance / lorList.size());
        lorList.clear();

        // avoid org.apache.commons.math3.exception.NotStrictlyPositiveException: standard deviation (0)
        return lorStd + 0.0001;
    }

    /**
     * Deprecated!
     * @param wesRecord Deprecated!
     * @param majorNC Deprecated!
     * @param minorNC Deprecated!
     * @return Deprecated!
     */
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
     * recalibrate p value with BH method, and get corresponding q value
     */
    private void bhRecalibrationOfEachPeak() {
        ArrayList<Map.Entry<Double, ArrayList<String>>> sortedPVals = new ArrayList<>(this.geneAsePValue.entrySet());
        // sort p value from small to large
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

        // sort genes with same p value
        for (Map.Entry<Double, ArrayList<String>> entry: sortedPVals) {
            Double pVal = entry.getKey();
            ArrayList<String> samePValGenes = entry.getValue();

            // sort by SNV numbers
            HashMap<String, Integer> samePValGeneSNVs = new HashMap<>();
            for (String gene: samePValGenes)
                samePValGeneSNVs.put(gene, this.geneSNVs.get(gene));
            // sort by major allele frequency(MAF)
            HashMap<String, Double> samePValGeneMajorAlleleFrequency = new HashMap<>();
            for (String gene: samePValGenes)
                samePValGeneMajorAlleleFrequency.put(gene, this.geneMajorAlleleFrequency.get(gene));

            List<Map.Entry<String, Integer>> samePValGeneEntry = new ArrayList<>(samePValGeneSNVs.entrySet());
            Collections.sort(samePValGeneEntry, new Comparator<Map.Entry<String, Integer>>() {
                @Override
                public int compare(Map.Entry<String, Integer> o1, Map.Entry<String, Integer> o2) {

                    // first sort genes with same p value with their MAF, then sort by SNV numbers if get same MAF
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
                rankage--;

                pValString = Double.toString(pVal);
                qValString = this.df.format(qValue);
                this.geneAseQValue.add(String.join("->", new String[]{geneName, pValString, qValString}));
            }
        }
    }

    /**
     * write the test result into file
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
            // sort items with its q value
            Collections.sort(records, new Comparator<Map.Entry<String, String[]>>() {
                @Override
                public int compare(Map.Entry<String, String[]> o1, Map.Entry<String, String[]> o2) {
                    String[] data1 = o1.getValue(), data2 = o2.getValue();
                    Double q1 = Double.parseDouble(data1[3]), q2 = Double.parseDouble(data2[3]);
                    if (!q1.equals(q2))
                        return q1.compareTo(q2);
                    // sort gene records with SNV numbers if has same q value
                    Integer snvCount1 = Integer.parseInt(data1[4]), snvCount2 = Integer.parseInt(data2[4]);
                    if (snvCount1 - snvCount2 != 0)
                        return snvCount2.compareTo(snvCount1);
                    // sort gene records with MAF if has same q value
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

        option = new Option("db", "dbsnp", true, "dbsnp file for filtering SNP site");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("o", "output", true, "ASE gene test output file, default ./aseGene.txt");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("df", "degree_of_freedom", true, "degree of freedom of inverse-Chi-square distribution, default 10");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("rc", "reads_coverage", true, "RNA-seq coverage threshold using for filter SNV records, default 10");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("wc", "wes_coverage", true, "WES coverage threshold using for filter SNV records, default 30");
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
