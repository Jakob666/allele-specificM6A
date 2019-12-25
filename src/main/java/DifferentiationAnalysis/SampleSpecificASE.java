package DifferentiationAnalysis;

import HierarchicalBayesianAnalysis.AseGeneDetection;
import HierarchicalBayesianAnalysis.HierarchicalBayesianModel;
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
import java.util.function.Function;
import java.util.stream.Collectors;

public class SampleSpecificASE {
    private String outputFile;
    private int samplingTime, burnIn, threadNumber;
    private double degreeOfFreedom, lorStd;
    private Logger logger;
    private HashMap<String, HashMap<Integer, String[]>> sample1GeneAlleleReads, sample2GeneAlleleReads;
    private AseGeneDetection sample1GeneDetector, sample2GeneDetector;
    private HashMap<String, int[][]> statisticForTest;
    private HashMap<String, int[]> geneTotalReads;
    private HashMap<String, Integer> genesSNVs;
    private HashMap<String, Double> geneMajorAlleleFrequency;
    private HashMap<Double, ArrayList<String>> geneSampleSpecificPValue;
    private HashMap<String, String> geneMajorNucleotide;
    private ArrayList<String> geneSampleSpecificQValue;
    private DecimalFormat df = new DecimalFormat("0.0000");
    private ReentrantLock lock;

    /**
     * Constructor
     * @param gtfFile GTF annotation file
     * @param sample1VcfFile sample1 VCF format file via MeRIP-seq INPUT data
     * @param sample1WesFile sample1 VCF format file via WES data, optional
     * @param sample2VcfFile sample2 VCF format file via MeRIP-seq INPUT data
     * @param sample2WesFile sample2 VCF format file via WES data, optional
     * @param dbsnpFile dbsnp annotation file for SNP filtering
     * @param outputFile result output file
     * @param degreeOfFreedom the degree of freedom of inverse-Chi-square distribution, default 10
     * @param readsCoverageThreshold reads coverage threshold when filter INPUT sample SNV sites, default 10
     * @param wesCoverageThreshold reads coverage threshold when filter WES SNV sites, default 30
     * @param samplingTime sampling time, default 10000
     * @param burnIn burn in time, default 2000
     * @param threadNumber thread number, default 2
     * @param logger log4j instance
     */
    public SampleSpecificASE(String gtfFile, String sample1VcfFile, String sample1WesFile, String sample2VcfFile, String sample2WesFile,
                             String dbsnpFile, String outputFile, double degreeOfFreedom, int readsCoverageThreshold, int wesCoverageThreshold,
                             int samplingTime, int burnIn, int threadNumber, Logger logger) {
        this.logger = logger;
        this.outputFile = outputFile;
        this.degreeOfFreedom = degreeOfFreedom;
        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
        this.threadNumber = threadNumber;
        String snvLocationFile = new File(new File(outputFile).getParent(), "sample1_snv_location.txt").getAbsolutePath();
        this.sample1GeneDetector = new AseGeneDetection(gtfFile, sample1VcfFile, sample1WesFile, dbsnpFile, snvLocationFile, outputFile, degreeOfFreedom, readsCoverageThreshold, wesCoverageThreshold, samplingTime, burnIn, threadNumber, logger);
        snvLocationFile = new File(new File(outputFile).getParent(), "sample2_snv_location.txt").getAbsolutePath();
        this.sample2GeneDetector = new AseGeneDetection(gtfFile, sample2VcfFile, sample2WesFile, dbsnpFile, snvLocationFile, outputFile, degreeOfFreedom, readsCoverageThreshold, wesCoverageThreshold, samplingTime, burnIn, threadNumber, logger);
    }

    public static void main(String[] args) {
        Options options = new Options();
        CommandLine commandLine = null;
        HelpFormatter help = new HelpFormatter();
        String header = "SampleSpecificASE contains following parameters: ";
        String footer = "";

        try {
            commandLine = setCommandLine(args, options);
        } catch (ParseException pe) {
            System.err.println(pe.getMessage());
            help.printHelp("java -jar renlabm6a_allele.jar SampleSpecificASE", header, options, footer, true);
            System.exit(2);
        }

        if (commandLine.hasOption("h")) {
            help.printHelp("java -jar renlabm6a_allele.jar SampleSpecificASE", header, options, footer, true);
            System.exit(0);
        }

        // default parameters
        String gtfFile = null, sample1VcfFile = null, sample1WesVcfFile = null, sample2VcfFile = null, sample2WesVcfFile = null,
               dbsnpFile = null, outputFile, outputDir;
        int samplingTime = 50000, burn_in = 5000, readsCoverageThreshold = 10, wesCoverageThreshold = 30, threadNumber = 2;
        double degreeOfFreedom = 10;

        if (!commandLine.hasOption("o"))
            outputFile = new File(System.getProperty("user.dir"), "sampleSpecificASE.txt").getAbsolutePath();
        else
            outputFile = commandLine.getOptionValue("o");
        outputDir = new File(outputFile).getParent();
        Logger logger = initLog(outputDir);

        if (!commandLine.hasOption("g")) {
            logger.error("GTF annotation file can not be empty");
            help.printHelp("java -jar renlabm6a_allele.jar SampleSpecificASE", header, options, footer, true);
            System.exit(2);
        } else {
            File gtf = new File(commandLine.getOptionValue("g"));
            if (!gtf.exists() || !gtf.isFile()) {
                logger.error("invalid file path: " + gtf.getAbsolutePath());
                System.exit(2);
            }
            gtfFile = gtf.getAbsolutePath();
        }

        if (!commandLine.hasOption("s1_vcf")) {
            logger.error("sample1 SNP calling VCF file can not be empty");
            help.printHelp("java -jar renlabm6a_allele.jar SampleSpecificASE", header, options, footer, true);
            System.exit(2);
        } else {
            File vcf = new File(commandLine.getOptionValue("s1_vcf"));
            if (!vcf.exists() || !vcf.isFile()) {
                logger.error("invalid file path: " + vcf.getAbsolutePath());
                System.exit(2);
            }
            sample1VcfFile = vcf.getAbsolutePath();
        }

        if (!commandLine.hasOption("s2_vcf")) {
            logger.error("sample2 SNP calling VCF file can not be empty");
            help.printHelp("java -jar renlabm6a_allele.jar SampleSpecificASE", header, options, footer, true);
            System.exit(2);
        } else {
            File vcf = new File(commandLine.getOptionValue("s2_vcf"));
            if (!vcf.exists() || !vcf.isFile()) {
                logger.error("invalid file path: " + vcf.getAbsolutePath());
                System.exit(2);
            }
            sample2VcfFile = vcf.getAbsolutePath();
        }

        if (commandLine.hasOption("s1_wes")) {
            File vcf = new File(commandLine.getOptionValue("s1_wes"));
            if (!vcf.exists() || !vcf.isFile()) {
                logger.error("invalid file path: " + vcf.getAbsolutePath());
                System.exit(2);
            }
            sample1WesVcfFile = vcf.getAbsolutePath();
        }

        if (commandLine.hasOption("s2_wes")) {
            File vcf = new File(commandLine.getOptionValue("s2_wes"));
            if (!vcf.exists() || !vcf.isFile()) {
                logger.error("invalid file path: " + vcf.getAbsolutePath());
                System.exit(2);
            }
            sample2WesVcfFile = vcf.getAbsolutePath();
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
        if (commandLine.hasOption("df"))
            degreeOfFreedom = Double.parseDouble(commandLine.getOptionValue("df"));
        if (degreeOfFreedom <= 1) {
            System.out.println("invalid inverse-Chi-square distribution parameter for tau sampling. Must larger than 1.0");
            System.exit(2);
        }
        if (commandLine.hasOption("t")) {
            if (Integer.valueOf(commandLine.getOptionValue("t")) <= 0) {
                System.err.println("invalid thread number, should be a positive integer");
                System.exit(2);
            }
            threadNumber = Integer.valueOf(commandLine.getOptionValue("t"));
        }

        SampleSpecificASE ssase = new SampleSpecificASE(gtfFile, sample1VcfFile, sample1WesVcfFile, sample2VcfFile, sample2WesVcfFile,
                dbsnpFile, outputFile, degreeOfFreedom, readsCoverageThreshold, wesCoverageThreshold, samplingTime, burn_in, threadNumber, logger);
        ssase.testSampleSpecificAse();
    }

    public void testSampleSpecificAse() {
        this.dataPreparation();
        this.calcLORStd();
        this.aseDifferentiationAnalysis();
        this.bhRecalibrationOfEachGene();
        this.output();
    }

    /**
     * sample specific ASE detection
     */
    private void dataPreparation() {
        this.sample1GeneDetector.dataPreparation();
        this.sample2GeneDetector.dataPreparation();
        this.logger.debug("find common genes between two samples");
        // geneAlleleReads = {"geneId->geneName": {pos1: [majorAllele:count, minorAllele: count1]}, ...}
        this.sample1GeneAlleleReads = this.sample1GeneDetector.getGeneAlleleReads();
        this.sample2GeneAlleleReads = this.sample2GeneDetector.getGeneAlleleReads();
        this.geneMajorNucleotide = new HashMap<>();
        HashSet<String> sample1Genes = new HashSet<> (sample1GeneAlleleReads.keySet()), sample2Genes = new HashSet<>(sample2GeneAlleleReads.keySet());
        sample1Genes.retainAll(sample2Genes);
        if (sample1Genes.isEmpty()) {
            this.logger.debug("detect no common genes for allele-specific expression differentiation analysis between two samples");
            System.exit(0);
        } else
            this.logger.debug("start data preparation");

        this.statisticForTest = new HashMap<>();
        HashMap<String, String> sample1GeneMajorNucleotide = this.sample1GeneDetector.getGeneMajorNucleotide();
        HashMap<String, String> sample2GeneMajorNucleotide = this.sample2GeneDetector.getGeneMajorNucleotide();

        Function<String, HashSet<Integer>> getGeneSNVPositions = (String record) -> {
            HashSet<Integer> positions = new HashSet<>();
            // record format: 41346154:A;41341577:G;41341541:G
            Arrays.stream(record.split(";")).forEach(rec -> positions.add(Integer.valueOf(rec.split(":")[0])));
            return positions;
        };
        int[][] sample1Reads, sample2Reads;
        int[] sample1MajorCounts, sample1MinorCounts, sample2MajorCounts, sample2MinorCounts;
        String[] info;
        String geneId;
        for (String gene: sample1Genes) {
            info = gene.split("->");
            geneId = info[0];
            // prepare data for differentiation analysis
            String sample1GeneRecord = sample1GeneMajorNucleotide.get(geneId);
            HashSet<Integer> sample1GeneSNVPositions = getGeneSNVPositions.apply(sample1GeneRecord);
            String sample2GeneRecord = sample2GeneMajorNucleotide.get(geneId);
            HashSet<Integer> sample2GeneSNVPositions = getGeneSNVPositions.apply(sample2GeneRecord);
            // remain common SNV sites
            sample1GeneSNVPositions.retainAll(sample2GeneSNVPositions);
            if (sample1GeneSNVPositions.isEmpty())
                continue;
            sample1Reads = this.getGeneSNVReadsCount(this.sample1GeneAlleleReads.get(gene), sample1GeneSNVPositions);
            sample1MajorCounts = sample1Reads[0];
            sample1MinorCounts = sample1Reads[1];
            sample2Reads = this.getGeneSNVReadsCount(this.sample2GeneAlleleReads.get(gene), sample1GeneSNVPositions);
            sample2MajorCounts = sample2Reads[0];
            sample2MinorCounts = sample2Reads[1];
            // common SNV sites major allele nucleotide
            this.formMajorAlleleRecord(gene, sample1GeneSNVPositions);

            this.statisticForTest.put(gene, new int[][] {sample1MajorCounts, sample1MinorCounts, sample2MajorCounts, sample2MinorCounts});
        }
        if (this.statisticForTest.isEmpty()) {
            this.logger.error("contains no common genes with SNV sites between two samples for differentiation test, please check the input data");
            System.exit(2);
        } else
            this.logger.debug("data preparation completed");
    }

    /**
     * calculate the standard deviation of LOR of all SNV sites on genome
     */
    private void calcLORStd() {
        this.logger.debug("calculate LOR standard as parameter of Inverse chi-square distribution");
        ArrayList<Double> lorList = new ArrayList<>();
        int[] sample1MajorCount, sample1MinorCount, sample2MajorCount, sample2MinorCount;
        double lor, cum = 0;
        for (String label: this.statisticForTest.keySet()) {
            int[][] statistic = this.statisticForTest.get(label);
            sample1MajorCount = statistic[0];
            sample1MinorCount = statistic[1];
            sample2MajorCount = statistic[2];
            sample2MinorCount = statistic[3];

            for (int i=0; i<sample1MajorCount.length; i++) {
                double sample1Major = sample1MajorCount[i], sample1Minor = sample1MinorCount[i],
                        sample2Major = sample2MajorCount[i], sample2Minor = sample2MinorCount[i];

                // Haldane's correction, adding 0.5 to all of the cells of a contingency table
                // if any of the cell expectations would cause a division by zero error.
                if (sample1Minor == 0 | sample2Minor == 0)
                    lor = ((sample1Major + 0.5) / (sample1Minor + 0.5)) / ((sample2Major + 0.5) / (sample2Minor + 0.5));
                else
                    lor = sample1Major / sample1Minor / (sample2Major / sample2Minor);
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

        this.lorStd = Math.sqrt(variance / lorList.size());
        lorList = null;

        // avoid org.apache.commons.math3.exception.NotStrictlyPositiveException: standard deviation (0)
        this.lorStd = (Math.abs(this.lorStd-0)<0.00001)? 0.0001: this.lorStd;
    }

    /**
     * run differentiation analysis
     */
    private void aseDifferentiationAnalysis() {
        this.logger.debug(this.statisticForTest.size() + " genes can be used for sample specific ASE test");
        this.lock = new ReentrantLock();
        this.geneSampleSpecificPValue = new HashMap<>();
        this.genesSNVs = new HashMap<>();
        this.geneMajorAlleleFrequency = new HashMap<>();
        this.geneTotalReads = new HashMap<>();
        ExecutorService threadPoolExecutor = Executors.newFixedThreadPool(this.threadNumber);
        CountDownLatch countDown = new CountDownLatch(this.statisticForTest.size());

        for (String label: this.statisticForTest.keySet()) {
            Runnable runnable = () -> {
                HierarchicalBayesianModel hb;
                int[][] geneReads;
                int[] sample1MajorReads, sample1MinorReads, sample2MajorReads, sample2MinorReads;
                try {
                    geneReads = this.statisticForTest.get(label);
                    sample1MajorReads = geneReads[0];
                    sample1MinorReads = geneReads[1];
                    sample2MajorReads = geneReads[2];
                    sample2MinorReads = geneReads[3];
                    int sample1Major = Arrays.stream(sample1MajorReads).reduce((x, y) -> x + y).getAsInt();
                    int sample1Minor = Arrays.stream(sample1MinorReads).reduce((x, y) -> x + y).getAsInt();
                    int sample2Major = Arrays.stream(sample2MajorReads).reduce((x, y) -> x + y).getAsInt();
                    int sample2Minor = Arrays.stream(sample2MinorReads).reduce((x, y) -> x + y).getAsInt();
                    this.geneTotalReads.put(label, new int[]{sample1Major, sample1Minor, sample2Major, sample2Minor});

                    if (sample1Minor == 0 && sample2Minor == 0) // gene shows ASE in two samples simultaneously
                        hb = new HierarchicalBayesianModel(this.lorStd, this.degreeOfFreedom,
                                this.samplingTime, this.burnIn, sample1MajorReads, sample1MajorReads, sample2MajorReads, sample2MajorReads);
                    else
                        hb = new HierarchicalBayesianModel(this.lorStd, this.degreeOfFreedom,
                                this.samplingTime, this.burnIn, sample1MajorReads, sample1MinorReads, sample2MajorReads, sample2MinorReads);

                    double p = hb.testSignificant();
                    double oddRatio = Math.exp(hb.quantifyGeneLOR());
                    double geneMAF = Math.min(1.0, oddRatio / (oddRatio + 1));
                    lock.lock();

                    ArrayList<String> samePValGenes = this.geneSampleSpecificPValue.getOrDefault(p, new ArrayList<>());
                    samePValGenes.add(label);
                    this.geneSampleSpecificPValue.put(p, samePValGenes);
                    this.genesSNVs.put(label, sample1MajorReads.length);
                    this.geneMajorAlleleFrequency.put(label, geneMAF);
                } catch (Exception e) {
                    this.logger.error("error occurs on record " + label);
                    this.logger.error(e.getMessage());
                } finally {
                    lock.unlock();
                    countDown.countDown();
                    hb = null;
                    geneReads = null;
                    sample1MajorReads = null;
                    sample1MinorReads = null;
                    sample2MajorReads = null;
                    sample2MinorReads = null;
                }
            };
            threadPoolExecutor.submit(runnable);
        }

        try {
            countDown.await();
        } catch (InterruptedException ie) {
            this.logger.error("analysis interrupted");
            this.logger.error(ie.getMessage());
        } finally {
            try {
                threadPoolExecutor.shutdown();
                if (!threadPoolExecutor.awaitTermination(1000, TimeUnit.MILLISECONDS))
                    threadPoolExecutor.shutdownNow();
            } catch (InterruptedException ie) {
                threadPoolExecutor.shutdownNow();
            }
            this.lock = null;
        }
        this.logger.debug("differentiation analysis complete");
    }

    /**
     * recalibrate p value with BH method, and get corresponding q value
     */
    private void bhRecalibrationOfEachGene() {
        this.logger.debug("start recalibrating sample specific ASE gene p values");
        this.logger.debug("sorting test result in order");
        ArrayList<Map.Entry<Double, ArrayList<String>>> sortedPVals = new ArrayList<>(this.geneSampleSpecificPValue.entrySet());
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
        this.geneSampleSpecificQValue = new ArrayList<>(totalGenes);
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
                this.geneSampleSpecificQValue.add(String.join("->", new String[]{geneName, pValString, qValString}));
            }
        }
        this.logger.debug("recalibration complete.");
    }

    /**
     * output sample specific test result
     */
    private void output() {
        HashMap<String, String[]> finalRecords = new HashMap<>();
        BufferedWriter bfw = null;
        try {
            bfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(this.outputFile))));
            String line, label, geneId, geneName, majorAlleleRecord, pVal, qVal;
            String[] info;
            int[] readsCount;
            int snvNum, sample1MajorCount, sample1MinorCount, sample2MajorCount, sample2MinorCount;
            bfw.write("#geneId\tgeneName\tp-value\tq-value\tsnvNum\tsample1Major/minorAlleleReads\tsample2Major/minorAlleleReads\tmajorAlleleNC\n");

            for (String record: this.geneSampleSpecificQValue) {
                info = record.split("->");
                geneId = info[0];
                geneName = info[1];
                pVal = info[2];
                qVal = info[3];
                label = String.join("->", new String[]{geneId, geneName});

                readsCount = this.geneTotalReads.get(label);
                sample1MajorCount = readsCount[0];
                sample1MinorCount = readsCount[1];
                sample2MajorCount = readsCount[2];
                sample2MinorCount = readsCount[3];
                majorAlleleRecord = this.geneMajorNucleotide.get(label);
                snvNum = this.genesSNVs.get(label);

                finalRecords.put(geneName, new String[]{geneId, geneName, pVal, qVal, Integer.toString(snvNum),
                        Integer.toString(sample1MajorCount), Integer.toString(sample1MinorCount),
                        Integer.toString(sample2MajorCount), Integer.toString(sample2MinorCount), majorAlleleRecord});
            }

            // sort items with its q value
            List<Map.Entry<String, String[]>> records = new ArrayList<>(finalRecords.entrySet());
            Collections.sort(records, new Comparator<Map.Entry<String, String[]>>() {
                @Override
                public int compare(Map.Entry<String, String[]> o1, Map.Entry<String, String[]> o2) {
                    String[] data1 = o1.getValue(), data2 = o2.getValue();
                    Double q1 = Double.parseDouble(data1[3]), q2 = Double.parseDouble(data2[3]);
                    if (!q1.equals(q2))
                        return q1.compareTo(q2);
                    Double p1 = Double.parseDouble(data1[2]), p2 = Double.parseDouble(data2[2]);
                    if (!p1.equals(p2))
                        return p1.compareTo(p2);
                    // sort gene records with SNV numbers if has same q value
                    Integer snvCount1 = Integer.parseInt(data1[4]), snvCount2 = Integer.parseInt(data2[4]);
                    return snvCount2.compareTo(snvCount1);
                }
            });
            for (Map.Entry<String, String[]> rec: records) {
                String[] data = rec.getValue();
                String[] lineInfo = new String[]{data[0], data[1], data[2], data[3], data[4], data[5] + "," + data[6], data[7]+ "," + data[8], data[9]};
                line = String.join("\t", lineInfo);
                bfw.write(line);
                bfw.newLine();
            }
            this.logger.debug("result file " + this.outputFile);
        } catch (IOException ie) {
            this.logger.error(ie.getMessage());
            System.exit(2);
        } finally {
            if (bfw != null) {
                try {
                    bfw.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
            this.statisticForTest = null;
            finalRecords = null;
        }
    }

    /**
     * get SNV reads count
     * @param records {pos1: [majorAllele:count, minorAllele: count1], ...}
     * @param geneSNVPosition positions list
     * @return major reads count and minor reads count
     */
    private int[][] getGeneSNVReadsCount(HashMap<Integer, String[]> records, HashSet<Integer> geneSNVPosition) {
        int count = geneSNVPosition.size();
        Integer[] positions = geneSNVPosition.toArray(new Integer[count]);
        int[] sampleMajorReads = new int[count], sampleMinorReads = new int[count];
        int pos;
        String majorRec, minorRec;
        for (int i=0; i<count; i++) {
            pos = positions[i];
            majorRec = records.get(pos)[0];
            minorRec = records.get(pos)[1];
            sampleMajorReads[i] = Integer.valueOf(majorRec.split(":")[1]);
            sampleMinorReads[i] = Integer.valueOf(minorRec.split(":")[1]);
        }
        positions = null;

        return new int[][] {sampleMajorReads, sampleMinorReads};
    }

    private void formMajorAlleleRecord(String gene, HashSet<Integer> positions) {
        List<Integer> sortedPositions = positions.stream().sorted().collect(Collectors.toList());
        ArrayList<String> rec = new ArrayList<>();
        for (Integer pos: sortedPositions) {
            String majorNC = this.sample1GeneAlleleReads.get(gene).get(pos)[0].split(":")[0];
            rec.add(pos+":"+majorNC);
        }
        String record = rec.stream().collect(Collectors.joining(";"));
        this.geneMajorNucleotide.put(gene, record);
    }

    private static Logger initLog(String logHome) {
        System.setProperty("log_home", logHome);
        return Logger.getLogger(SampleSpecificASE.class);
    }

    private static CommandLine setCommandLine(String[] args, Options options) throws ParseException {
        Option option = new Option("s1_vcf", "sample1_vcf_file", true, "VCF format file generate by sample1 RNA-seq or MeRIP-seq data SNP calling process, required");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("s1_wes", "sample1_wes_vcf_file", true, "VCF format file generate by sample1 WES data SNP calling process, optional");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("s2_vcf", "sample2_vcf_file", true, "VCF format file generate by sample2 RNA-seq or MeRIP-seq data SNP calling process, required");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("s2_wes", "sample2_wes_vcf_file", true, "VCF format file generate by sample2 WES data SNP calling process, optional");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("g", "gtf", true, "GTF annotation file, required");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("db", "dbsnp", true, "big scale SNV annotation data set, like dbsnp, 1000Genome etc. Optional, the file format see https://github.com/Jakob666/allele-specificM6A)");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("o", "output", true, "ASE gene test output file. Optional, default ./sampleSpecificASE.txt");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("df", "degree_of_freedom", true, "degree of freedom of inverse-Chi-square distribution. Optional, default 10");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("rc", "reads_coverage", true, "reads coverage threshold using for filter RNA-seq or MeRIP-seq data SNVs in VCF file (aim for reducing FDR). Optional, default 10");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("bc", "bkg_coverage", true, "reads coverage threshold using for filter WES data SNVs in VCF file (aim for reducing FDR). Optional, default 30");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("s", "sampling", true, "sampling times, larger than 500. Optional, default 50000");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("b", "burn", true, "burn-in times, more than 100 and less than sampling times. Optional, default 5000");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("t", "thread", true, "thread number for running test. Optional, default 2");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("h", "help", false, "help message of SampleSpecificASE");
        option.setRequired(false);
        options.addOption(option);

        CommandLineParser parser = new DefaultParser();

        return parser.parse(options, args);
    }
}
