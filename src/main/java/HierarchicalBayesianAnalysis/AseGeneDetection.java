package HierarchicalBayesianAnalysis;

import GTFComponent.GeneSNVRecord;
import GTFComponent.VcfSnpMatchGene;
import heterozygoteSiteAnalysis.DbsnpAnnotation;
import org.apache.commons.cli.*;
import org.apache.log4j.Logger;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.locks.ReentrantLock;

/**
 * Test ASE genes with VCF format file
 */
public class AseGeneDetection {
    private String aseGeneFile;
    // geneAlleleReads = {"geneId->geneName": {pos1: [majorAllele:count, minorAllele: count1]}, ...}
    private HashMap<String, HashMap<Integer, String[]>> geneAlleleReads, geneBackgroundReads;
    private HashMap<String, HashMap<String, Integer>> geneReadsCount = new HashMap<>(), geneBackgroundCount = new HashMap<>();
    private int samplingTime, burnIn, threadNumber;
    private HashMap<Double, ArrayList<String>> geneAsePValue = new HashMap<>();
    private HashMap<String, Integer> genesSNVs = new HashMap<>();
    private HashMap<String, Double> geneMajorAlleleFrequency = new HashMap<>();
    private HashMap<String, String> geneMajorNucleotide = new HashMap<>();
    private HashMap<String, ArrayList<int[]>> statisticForTest = new HashMap<>();
    private ArrayList<String> geneAseQValue;
    private DecimalFormat df = new DecimalFormat("0.0000");
    private ReentrantLock lock;
    private Logger logger;

    /**
     * Constructor
     * @param gtfFile GTF annotation file
     * @param vcfFile VCF format file via MeRIP-seq INPUT data
     * @param wesFile VCF format file via WES data, optional
     * @param dbsnpFile dbsnp annotation file for SNP filtering
     * @param aseGeneFile test result output file
     * @param readsCoverageThreshold reads coverage threshold when filter INPUT sample SNV sites, default 10
     * @param wesCoverageThreshold reads coverage threshold when filter WES SNV sites, default 30
     * @param samplingTime sampling time, default 10000
     * @param burnIn burn in time, default 2000
     * @param threadNumber thread number, default 2
     * @param logger log4j instance
     */
    public AseGeneDetection(String gtfFile, String vcfFile, String wesFile, String dbsnpFile, String aseGeneFile,
                            int readsCoverageThreshold, int wesCoverageThreshold,
                            int samplingTime, int burnIn, int threadNumber, Logger logger) {
        HashMap<String, LinkedList<DbsnpAnnotation.DIYNode>> dbsnpRecord = null;
        if (dbsnpFile != null) {
            DbsnpAnnotation da = new DbsnpAnnotation(dbsnpFile, logger);
            da.parseDbsnpFile();
            dbsnpRecord = da.getDbsnpRecord();
        }
        VcfSnpMatchGene vsmg = new VcfSnpMatchGene(vcfFile, gtfFile, readsCoverageThreshold);
        vsmg.parseVcfFile(dbsnpRecord);
        this.geneAlleleReads = vsmg.getGeneAlleleReads();
        if (wesFile != null) {
            vsmg = new VcfSnpMatchGene(wesFile, gtfFile, wesCoverageThreshold);
            vsmg.parseVcfFile(dbsnpRecord);
            this.geneBackgroundReads = vsmg.getGeneAlleleReads();
        } else {
            this.geneBackgroundReads = null;
        }
        if (dbsnpRecord != null)
            dbsnpRecord = null;

        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
        this.aseGeneFile = aseGeneFile;
        this.logger = logger;
        if (this.aseGeneFile == null) {
            this.logger.error("output file path can not be null");
            System.exit(2);
        }
        this.logger.debug("locate SNP record in " + vcfFile + " to corresponding genes");
        String snvLocationFile = new File(new File(aseGeneFile).getParent(), "snv_location.txt").getAbsolutePath();
        GeneSNVRecord gsr = new GeneSNVRecord(gtfFile, vcfFile, snvLocationFile);
        gsr.locateSnv();
        gsr = null;
        this.threadNumber = threadNumber;
        this.logger.debug("SNP locations in " + snvLocationFile);
    }

    /**
     * Constructor
     * @param gtfFile GTF annotation file
     * @param vcfFile VCF format file via MeRIP-seq INPUT data
     * @param wesFile VCF format file via WES data, optional
     * @param dbsnpFile dbsnp annotation file for SNP filtering
     * @param snvLocationFile SNVs location file
     * @param aseGeneFile test result output file
     * @param readsCoverageThreshold reads coverage threshold when filter INPUT sample SNV sites, default 10
     * @param wesCoverageThreshold reads coverage threshold when filter WES SNV sites, default 30
     * @param samplingTime sampling time, default 10000
     * @param burnIn burn in time, default 2000
     * @param threadNumber thread number, default 2
     * @param onlyExonMutation only remain SNPs locate in exon region, default false
     * @param logger log4j instance
     */
    public AseGeneDetection(String gtfFile, String vcfFile, String wesFile, String dbsnpFile, String snvLocationFile,
                            String aseGeneFile, int readsCoverageThreshold, int wesCoverageThreshold,
                            int samplingTime, int burnIn, int threadNumber, Logger logger) {
        HashMap<String, LinkedList<DbsnpAnnotation.DIYNode>> dbsnpRecord = null;
        if (dbsnpFile != null) {
            DbsnpAnnotation da = new DbsnpAnnotation(dbsnpFile, logger);
            da.parseDbsnpFile();
            dbsnpRecord = da.getDbsnpRecord();
        }
        VcfSnpMatchGene vsmg = new VcfSnpMatchGene(vcfFile, gtfFile, readsCoverageThreshold);
        vsmg.parseVcfFile(dbsnpRecord);
        this.geneAlleleReads = vsmg.getGeneAlleleReads();
        if (wesFile != null) {
            vsmg = new VcfSnpMatchGene(wesFile, gtfFile, wesCoverageThreshold);
            vsmg.parseVcfFile(dbsnpRecord);
            this.geneBackgroundReads = vsmg.getGeneAlleleReads();
        } else {
            this.geneBackgroundReads = null;
        }
        if (dbsnpRecord != null)
            dbsnpRecord = null;

        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
        this.aseGeneFile = aseGeneFile;
        this.logger = logger;
        if (this.aseGeneFile == null) {
            this.logger.error("output file path can not be null");
            System.exit(2);
        }
        this.logger.debug("locate SNP record in " + vcfFile + " to corresponding genes");
        GeneSNVRecord gsr = new GeneSNVRecord(gtfFile, vcfFile, snvLocationFile);
        gsr.locateSnv();
        gsr = null;
        this.threadNumber = threadNumber;
        this.logger.debug("SNP locations in " + snvLocationFile);
    }

    public static void main(String[] args) {
        Options options = new Options();
        CommandLine commandLine = null;
        HelpFormatter help = new HelpFormatter();
        String header = "AseGeneDetection contains following parameters: ";
        String footer = "";

        try {
            commandLine = setCommandLine(args, options);
        } catch (ParseException pe) {
            System.err.println(pe.getMessage());
            help.printHelp("java -jar renlabm6a_allele.jar AseGeneDetection", header, options, footer, true);
            System.exit(2);
        }

        if (commandLine.hasOption("h")) {
            help.printHelp("java -jar renlabm6a_allele.jar AseGeneDetection", header, options, footer, true);
            System.exit(0);
        }

        // default parameters
        String gtfFile = null, aseVcfFile = null, wesVcfFile = null, dbsnpFile = null, outputFile, outputDir;
        int samplingTime = 50000, burn_in = 10000, readsCoverageThreshold = 5, wesCoverageThreshold = 30, threadNumber = 2;

        if (!commandLine.hasOption("o"))
            outputFile = new File(System.getProperty("user.dir"), "aseGene.txt").getAbsolutePath();
        else
            outputFile = commandLine.getOptionValue("o");
        outputDir = new File(outputFile).getParent();
        Logger logger = initLog(outputDir);

        if (!commandLine.hasOption("g")) {
            logger.error("GTF annotation file can not be empty");
            help.printHelp("java -jar renlabm6a_allele.jar AseGeneDetection", header, options, footer, true);
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
            help.printHelp("java -jar renlabm6a_allele.jar AseGeneDetection", header, options, footer, true);
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
        if (commandLine.hasOption("bc"))
            wesCoverageThreshold = Integer.valueOf(commandLine.getOptionValue("bc"));
        if (commandLine.hasOption("s"))
            samplingTime = Integer.parseInt(commandLine.getOptionValue("s"));
        if (commandLine.hasOption("b"))
            burn_in = Integer.parseInt(commandLine.getOptionValue("b"));
        if (samplingTime <= 500 || burn_in <= 100) {
            logger.error("sampling times larger than 500 and burn in times at least 100");
            System.exit(2);
        }
        if (commandLine.hasOption("t")) {
            if (Integer.valueOf(commandLine.getOptionValue("t")) <= 0) {
                System.err.println("invalid thread number, should be a positive integer");
                System.exit(2);
            }
            threadNumber = Integer.valueOf(commandLine.getOptionValue("t"));
        }

        AseGeneDetection agd = new AseGeneDetection(gtfFile, aseVcfFile, wesVcfFile, dbsnpFile, outputFile,
                readsCoverageThreshold, wesCoverageThreshold, samplingTime, burn_in, threadNumber, logger);
        agd.getTestResult();
    }

    public HashMap<String, String> getGeneMajorNucleotide() {
        return geneMajorNucleotide;
    }

    public HashMap<String, HashMap<Integer, String[]>> getGeneAlleleReads() {
        return this.geneAlleleReads;
    }

    public void getTestResult() {
        this.dataPreparation();
        this.aseGeneTest();
        this.bhRecalibrationOfEachGene();
        this.outputResult();
    }

    /**
     * data preparation for hierarchical test
     */
    public void dataPreparation() {
        HashMap<Integer, String[]> geneSNVs, wesSNVs;
        ArrayList<Integer> majorAlleleCount, minorAlleleCount, majorBackgroundCount, minorBackgroundCount;
        int[] majorCount, minorCount, majorBackground, minorBackground;
        String geneId;
        // geneAlleleReads = {"geneId->geneName": {pos1: [majorAllele:count, minorAllele: count]}, ...}
        this.logger.debug("calculate LOR for genes to be tested");
        for (String label : this.geneAlleleReads.keySet()) {
            // WES data exists, but there is no WES SNV sites locate in the gene
            if (this.geneBackgroundReads != null) {
                if (this.geneBackgroundReads.keySet().contains(label))
                    wesSNVs = this.geneBackgroundReads.getOrDefault(label, null);
                else
                    continue;
            } else
                wesSNVs = null;

            geneId = label.split("->")[0];
            majorAlleleCount = new ArrayList<>();
            minorAlleleCount = new ArrayList<>();
            majorBackgroundCount = new ArrayList<>();
            minorBackgroundCount = new ArrayList<>();

            // {pos1: [majorAllele:count, minorAllele: count], pos2: [majorAllele: count, minorAllele: count]...}
            geneSNVs = this.geneAlleleReads.get(label);

            // record the genome positions of major alleles and the corresponding nucleotides
            LinkedList<String> majorNcRecords = new LinkedList<>();
            for (Integer mutPosition : geneSNVs.keySet()) {

                if (wesSNVs != null && !wesSNVs.keySet().contains(mutPosition))
                    continue;
                // [majorAllele:count, minorAllele: count1]
                String[] nucleotideReadsCount = geneSNVs.get(mutPosition);
                String majorAlleleRecord, minorAlleleRecord, majorNC, minorNC;
                majorAlleleRecord = nucleotideReadsCount[0];
                minorAlleleRecord = nucleotideReadsCount[1];
                majorNC = majorAlleleRecord.split(":")[0];
                minorNC = minorAlleleRecord.split(":")[0];
                int major = Integer.valueOf(majorAlleleRecord.split(":")[1]);
                int minor = Integer.valueOf(minorAlleleRecord.split(":")[1]);
                majorAlleleCount.add(major);
                minorAlleleCount.add(minor);
                majorNcRecords.add(String.join(":", new String[]{Integer.toString(mutPosition), majorNC}));

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
                minorCount[i] = minorAlleleCount.get(i);

                majorBackground[i] = majorBackgroundCount.get(i);
                minorBackground[i] = minorBackgroundCount.get(i);
            }

            ArrayList<int[]> statistic = new ArrayList<>(4);
            statistic.add(majorCount);
            statistic.add(minorCount);
            statistic.add(majorBackground);
            statistic.add(minorBackground);
            this.statisticForTest.put(label, statistic);
        }
        if (this.statisticForTest.isEmpty()) {
            this.logger.error("contains no genes with SNV sites for hierarchical test, please check the input data");
            System.exit(2);
        }
    }

    /**
     * get ASE significant p value with hierarchical model
     */
    private void aseGeneTest() {
        // test ASE gene with Hierarchical model
        this.logger.debug("hierarchical Bayesian model test start");
        // init thread pool and lock
        ExecutorService threadPoolExecutor = Executors.newFixedThreadPool(this.threadNumber);
        this.lock = new ReentrantLock();
        CountDownLatch countDown = new CountDownLatch(this.statisticForTest.size());
        long tenPercent = Math.round(this.statisticForTest.size() * 0.1);

        RunTest task = (String name) -> {
            return new Runnable() {
                @Override
                public void run() {
                    ArrayList<int[]> statistic;
                    double df, scaleParam, aveDepth, epsilon = 0.00000000001, mafLowerBound = 0.49, asLowBound = 0.59;
                    int[] majorCount, minorCount, majorBackground, minorBackground;
                    HierarchicalBayesianModel hb;
                    try {
                        statistic = statisticForTest.get(name);
                        majorCount = statistic.get(0);
                        minorCount = statistic.get(1);
                        majorBackground = statistic.get(2);
                        minorBackground = statistic.get(3);
                        if (majorCount.length == 1 && minorCount[0] == 0 && majorCount[0] <= 35)
                            df = 2;
                        else
                            df = Math.max(3, majorCount.length);
                        aveDepth = Arrays.stream(majorCount).average().getAsDouble();
                        if (majorCount.length == 1)
                            scaleParam = 50;
                        else
                            scaleParam = (aveDepth - 15 < epsilon)? 50: 100;
                        hb = new HierarchicalBayesianModel(df, scaleParam, samplingTime, burnIn,
                                majorCount, minorCount, majorBackground, minorBackground);
                        boolean wellDone = false;
                        int maxTrail = 0;
                        double p = 0, geneOddRatio, geneMAF = 0;
                        while (!wellDone && maxTrail < 20) {    // avoid MCMC failed
                            p = hb.testSignificant();
                            geneOddRatio = Math.exp(hb.quantifyGeneLOR());
                            geneMAF = Math.min(1.0, geneOddRatio / (geneOddRatio + 1));
                            Double maf = new Double(geneMAF);
                            if (!maf.equals(Double.NaN) && !(!hb.isAllZero() && Math.abs(geneMAF - 1.0) < epsilon) &&
                                !(geneMAF - mafLowerBound < epsilon) && !(geneMAF-asLowBound>epsilon && p - 0.05 > epsilon))
                                wellDone = true;
                            maxTrail++;
                        }

                        lock.lock();

                        ArrayList<String> samePValGenes = geneAsePValue.getOrDefault(p, new ArrayList<>());
                        samePValGenes.add(name);
                        geneAsePValue.put(p, samePValGenes);
                        genesSNVs.put(name, majorCount.length);
                        geneMajorAlleleFrequency.put(name, geneMAF);
                    } catch (Exception e) {
                        logger.error("error occurs on record " + name);
                        e.printStackTrace();
                        logger.error(e.getMessage());
                    } finally {
                        countDown.countDown();
                        if (countDown.getCount() % tenPercent == 0) {
                            double proportion = 100 - 10.0 * countDown.getCount() / tenPercent;
                            if (proportion >= 0)
                                logger.debug(proportion + "% completed");
                        }
                        lock.unlock();
                        hb = null;
                        majorCount = null;
                        minorCount = null;
                        majorBackground = null;
                        minorBackground = null;
                    }
                }
            };
        };

        this.logger.debug(this.statisticForTest.size() + " genes to be tested");
        for (String label: this.statisticForTest.keySet()) {
            Runnable runnable = task.runTask(label);
            threadPoolExecutor.submit(runnable);
        }
        try {
            countDown.await();
        } catch (InterruptedException ie) {
            this.logger.error("analysis interrupted");
            this.logger.error(ie.getMessage());
        } finally {
            this.statisticForTest = null;
            try {
                threadPoolExecutor.shutdown();
                if (!threadPoolExecutor.awaitTermination(1000, TimeUnit.MILLISECONDS))
                    threadPoolExecutor.shutdownNow();
            } catch (InterruptedException ie) {
                threadPoolExecutor.shutdownNow();
            }
            this.lock = null;
        }
        this.logger.debug("model test complete");
    }

    /**
     * recalibrate p value with BH method, and get corresponding q value
     */
    private void bhRecalibrationOfEachGene() {
        this.logger.debug("start recalibrating p values of hierarchical model");
        this.logger.debug("sorting test result in order");
        ArrayList<Map.Entry<Double, ArrayList<String>>> sortedPVals = new ArrayList<>(this.geneAsePValue.entrySet());
        // sort p value from large to small
        Collections.sort(sortedPVals, new Comparator<Map.Entry<Double, ArrayList<String>>>() {
            @Override
            public int compare(Map.Entry<Double, ArrayList<String>> o1, Map.Entry<Double, ArrayList<String>> o2) {
                return o2.getKey().compareTo(o1.getKey());
            }
        });

        int totalGenes = this.genesSNVs.size(), rankage = totalGenes;
        double prevQValue = 1.0, qValue;
        String pValString, qValString;
        this.geneAseQValue = new ArrayList<>(totalGenes);

        for (Map.Entry<Double, ArrayList<String>> entry: sortedPVals) {
            Double pVal = entry.getKey();
            ArrayList<String> samePValGenes = entry.getValue();
            // sort by SNV numbers
            HashMap<String, Integer> samePValGeneSNVs = new HashMap<>();
            for (String gene: samePValGenes)
                samePValGeneSNVs.put(gene, this.genesSNVs.get(gene));
            // sort by major allele frequency(MAF)
            HashMap<String, Double> samePValGeneMajorAlleleFrequency = new HashMap<>();
            for (String gene: samePValGenes)
                samePValGeneMajorAlleleFrequency.put(gene, this.geneMajorAlleleFrequency.get(gene));

            List<Map.Entry<String, Integer>> samePValGeneEntry = new ArrayList<>(samePValGeneSNVs.entrySet());
            Collections.sort(samePValGeneEntry, new Comparator<Map.Entry<String, Integer>>() {
                @Override
                public int compare(Map.Entry<String, Integer> o1, Map.Entry<String, Integer> o2) {

                    // first sort genes with same p value with their SNV number, then sort by MAF if get same SNV number
                    // both SNV number and MAF are sorted from large to small
                    String gene1 = o1.getKey(), gene2 = o2.getKey();
                    Double gene1MAF = samePValGeneMajorAlleleFrequency.get(gene1), gene2MAF = samePValGeneMajorAlleleFrequency.get(gene2);
                    Integer gene1SNVs = o1.getValue(), gene2SNVs = o2.getValue();
                    if (gene1SNVs.equals(gene2SNVs)) {
                        return gene2MAF.compareTo(gene1MAF);
                    } else
                        return gene2SNVs.compareTo(gene1SNVs);
                }
            });

            for (Map.Entry<String, Integer> geneEntry: samePValGeneEntry) {
                String geneName = geneEntry.getKey();
                qValue = Math.min(prevQValue, pVal * totalGenes / rankage);
                if (qValue - prevQValue < 0.00001)
                    prevQValue = qValue;
                rankage--;

                pValString = this.df.format(pVal);
                qValString = this.df.format(qValue);
                this.geneAseQValue.add(String.join("->", new String[]{geneName, pValString, qValString}));
            }
        }
        this.logger.debug("recalibration complete.");
    }

    /**
     * write the test result into file
     */
    private void outputResult() {
        HashMap<String, String[]> finalRecords = new HashMap<>();
        BufferedWriter bfw = null;
        try {
            bfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(this.aseGeneFile))));
            String line, label, geneId, geneName, majorAlleleRecord, pVal, qVal;
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
                majorAlleleFrequency = this.geneMajorAlleleFrequency.get(label);
                majorAlleleRecord = this.geneMajorNucleotide.get(geneId);

                snvNum = this.genesSNVs.get(label);
                majorBackground = (this.geneBackgroundReads == null)? 0: this.geneBackgroundCount.get(label).get("major");
                minorBackground = (this.geneBackgroundReads == null)? 0: this.geneBackgroundCount.get(label).get("minor");

                finalRecords.put(label, new String[]{geneId, geneName, pVal, qVal, Integer.toString(snvNum),
                        Integer.toString(majorCount), Integer.toString(minorCount), Integer.toString(majorBackground),
                        Integer.toString(minorBackground), this.df.format(majorAlleleFrequency), majorAlleleRecord});
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
            this.logger.debug("result file " + this.aseGeneFile);
        } catch (IOException ie) {
            this.logger.error(ie.getMessage());
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
        Option option = new Option("vcf", "vcf_file", true, "VCF format file generate by RNA-seq or MeRIP-seq data SNP calling process, required");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("wes", "wes_vcf_file", true, "VCF format file generate by WES data SNP calling process, optional");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("g", "gtf", true, "GTF annotation file, required");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("db", "dbsnp", true, "big scale SNV annotation data set, like dbsnp, 1000Genome etc. Optional, the file format see https://github.com/Jakob666/allele-specificM6A");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("o", "output", true, "ASE gene test output file. Optional, default ./aseGene.txt");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("rc", "reads_coverage", true, "reads coverage threshold using for filter RNA-seq or MeRIP-seq data SNVs in VCF file (aim for reducing FDR). Optional, default 5");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("bc", "bkg_coverage", true, "reads coverage threshold using for filter WES data SNVs in VCF file (aim for reducing FDR). Optional, default 30");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("s", "sampling", true, "sampling times, larger than 500. Optional, default 50000");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("b", "burn", true, "burn-in times, more than 100 and less than sampling times. Optional, default 10000");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("t", "thread", true, "thread number for running test. Optional, default 2");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("h", "help", false, "help message of AseGeneDetection");
        option.setRequired(false);
        options.addOption(option);

        CommandLineParser parser = new DefaultParser();

        return parser.parse(options, args);
    }
}
