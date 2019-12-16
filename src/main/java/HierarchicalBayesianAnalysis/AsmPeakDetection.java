package HierarchicalBayesianAnalysis;


import AseM6aPeakDetector.HeterozygoteReadsCount;
import heterozygoteSiteAnalysis.DbsnpAnnotation;
import heterozygoteSiteAnalysis.PeakWithSNV;
import org.apache.commons.cli.*;
import org.apache.log4j.Logger;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.locks.ReentrantLock;

/**
 * Test ASM peaks with VCF format file and BED format file
 */
public class AsmPeakDetection {
    private String gtfFile, peakBedFile, wesFile, asmPeakFile, peakCoveredSnpFile, peakCoveredWesSnpFile;
    private int ipSNPReadInfimum, wesSNPReadInfimum, samplingTime, burnIn, totalPeakCount = 0, threadNumber;
    private double degreeOfFreedom, lorStd;
    private Logger log;
    private HashMap<String, HashMap<String, HashMap<String, HashMap<String, Integer>>>> peakSnpReadsCount, peakSnpBackground;
    private HashMap<Double, ArrayList<String>> asmPValue = new HashMap<>();
    private HashMap<String, HashMap<String, Integer>> peakMajorMinorAlleleCount = new HashMap<>(), peakMajorMinorBackground = new HashMap<>();
    private HashMap<String, Double> peakMajorAlleleFrequency = new HashMap<>();
    private HashMap<String, String> peakCoveredGene = new HashMap<>();
    private HashMap<String, ArrayList<int[]>> statisticForTest = new HashMap<>();
    private HashMap<String, LinkedList<String>> peakMajorAlleleNucleotide;
    private HashMap<String, Integer> peakSNVNum = new HashMap<>();
    private HashMap<String, LinkedList<DbsnpAnnotation.DIYNode>> dbsnpRecord;
    private HashMap<String, HashMap<String, String>> geneNames;
    private HashMap<String, HashSet<Integer>> snvForTest;
    private ArrayList<String> asmQValue = new ArrayList<>();
    private DecimalFormat df = new DecimalFormat("0.0000");
    private ReentrantLock lock;

    /**
     * Constructor
     * @param gtfFile GTF annotation file
     * @param peakBedFile BED format file via MeRIP-seq IP data
     * @param vcfFile VCF format file via MeRIP-seq INPUT data
     * @param wesFile VCF format file via WES data, optional
     * @param dbsnpFile dbsnp file for SNP filtering
     * @param asmPeakFile test result output file
     * @param degreeOfFreedom the degree of freedom of inverse-Chi-square distribution, default 10
     * @param ipSNPReadInfimum reads coverage threshold when filter INPUT sample SNV sites, default 10
     * @param wesSNPReadInfimum reads coverage threshold when filter WES SNV sites, default 30
     * @param samplingTime sampling time, default 10000
     * @param burnIn burn in time, default 2000
     * @param threadNumber hread number, default 2
     * @param log log4j instance
     */
    public AsmPeakDetection(String gtfFile, String peakBedFile, String vcfFile, String wesFile, String dbsnpFile,
                            String asmPeakFile, double degreeOfFreedom, int ipSNPReadInfimum, int wesSNPReadInfimum,
                            int samplingTime, int burnIn, int threadNumber, Logger log) {
        this.gtfFile = gtfFile;
        this.peakBedFile = peakBedFile;
        this.wesFile = wesFile;
        String outputDir = new File(asmPeakFile).getParent();
        this.log = log;
        if (dbsnpFile != null) {
            DbsnpAnnotation da = new DbsnpAnnotation(dbsnpFile, this.log);
            da.parseDbsnpFile();
            this.dbsnpRecord = da.getDbsnpRecord();
        } else
            this.dbsnpRecord = null;
        this.asmPeakFile = asmPeakFile;
        this.degreeOfFreedom = degreeOfFreedom;
        this.ipSNPReadInfimum = ipSNPReadInfimum;
        this.wesSNPReadInfimum = wesSNPReadInfimum;
        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
        this.threadNumber = threadNumber;

        this.log.debug("locate SNP record in " + vcfFile + " to corresponding m6A signal peaks");
        this.peakCoveredSnpFile = new File(outputDir, "peak_with_snv.txt").getAbsolutePath();
        PeakWithSNV pws = new PeakWithSNV(gtfFile, peakBedFile, vcfFile, this.peakCoveredSnpFile, true);
        pws.locateSnvInPeak();
        this.log.debug("SNP locations in " + this.peakCoveredSnpFile);
        this.peakCoveredWesSnpFile = new File(outputDir, "peak_with_snv_bkg.txt").getAbsolutePath();
        if (wesFile != null) {
            this.log.debug("locate WES data SNP record in " + vcfFile + " to corresponding m6A signal peaks");
            pws = new PeakWithSNV(gtfFile, peakBedFile, vcfFile, this.peakCoveredWesSnpFile, true);
            pws.locateSnvInPeak();
            this.log.debug("SNP locations in " + this.peakCoveredWesSnpFile);
        }
        pws = null;
        this.snvForTest = new HashMap<>();
    }

    /**
     * Constructor
     * @param gtfFile GTF annotation file
     * @param peakBedFile BED format file via MeRIP-seq IP data
     * @param vcfFile VCF format file via MeRIP-seq INPUT data
     * @param wesFile VCF format file via WES data, optional
     * @param dbsnpFile dbsnp file for SNP filtering
     * @param peakWithSnvFile peak contains SNVs
     * @param asmPeakFile test result output file
     * @param degreeOfFreedom the degree of freedom of inverse-Chi-square distribution, default 10
     * @param ipSNPReadInfimum reads coverage threshold when filter INPUT sample SNV sites, default 10
     * @param wesSNPReadInfimum reads coverage threshold when filter WES SNV sites, default 30
     * @param samplingTime sampling time, default 10000
     * @param burnIn burn in time, default 2000
     * @param threadNumber hread number, default 2
     * @param log log4j instance
     */
    public AsmPeakDetection(String gtfFile, String peakBedFile, String vcfFile, String wesFile, String dbsnpFile,
                            String peakWithSnvFile, String peakWithSnvBkgFile, String asmPeakFile, double degreeOfFreedom, int ipSNPReadInfimum, int wesSNPReadInfimum,
                            int samplingTime, int burnIn, int threadNumber, Logger log) {
        this.gtfFile = gtfFile;
        this.peakBedFile = peakBedFile;
        this.wesFile = wesFile;
        this.log = log;
        if (dbsnpFile != null) {
            DbsnpAnnotation da = new DbsnpAnnotation(dbsnpFile, this.log);
            da.parseDbsnpFile();
            this.dbsnpRecord = da.getDbsnpRecord();
        } else
            this.dbsnpRecord = null;
        this.asmPeakFile = asmPeakFile;
        this.degreeOfFreedom = degreeOfFreedom;
        this.ipSNPReadInfimum = ipSNPReadInfimum;
        this.wesSNPReadInfimum = wesSNPReadInfimum;
        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
        this.threadNumber = threadNumber;

        this.log.debug("locate SNP record in " + vcfFile + " to corresponding m6A signal peaks");
        this.peakCoveredSnpFile = peakWithSnvFile;
        PeakWithSNV pws = new PeakWithSNV(gtfFile, peakBedFile, vcfFile, this.peakCoveredSnpFile, true);
        pws.locateSnvInPeak();
        this.log.debug("SNP locations in " + this.peakCoveredSnpFile);
        this.peakCoveredWesSnpFile = peakWithSnvBkgFile;
        if (wesFile != null) {
            this.log.debug("locate WES data SNP record in " + vcfFile + " to corresponding m6A signal peaks");
            pws = new PeakWithSNV(gtfFile, peakBedFile, vcfFile, this.peakCoveredWesSnpFile, true);
            pws.locateSnvInPeak();
            this.log.debug("SNP locations in " + this.peakCoveredWesSnpFile);
        }
        pws = null;
        this.snvForTest = new HashMap<>();
    }

    public static void main(String[] args) {
        Options options = new Options();
        CommandLine commandLine = null;
        HelpFormatter help = new HelpFormatter();
        String header = "AsmPeakDetection contains following parameters: ";
        String footer = "";

        try {
            commandLine = setCommandLine(args, options);
        } catch (ParseException pe) {
            System.err.println(pe.getMessage());
            help.printHelp("java -jar renlabm6a_allele.jar AsmPeakDetection", header, options, footer, true);
            System.exit(2);
        }

        if (commandLine.hasOption("h")) {
            help.printHelp("java -jar renlabm6a_allele.jar AsmPeakDetection", header, options, footer, true);
            System.exit(0);
        }

        // default parameters
        String gtfFile = null, bedFile = null, aseVcfFile = null, wesVcfFile = null, dbsnpFile = null, outputFile, outputDir;
        int ipSNPCoverageInfimum = 10, wesSNPCoverageInfimum = 30, samplingTime = 10000, burn_in = 2000, threadNumber = 2;
        double degreeOfFreedom = 10;

        if (!commandLine.hasOption("o"))
            outputFile = new File(System.getProperty("user.dir"), "asmPeak.txt").getAbsolutePath();
        else
            outputFile = commandLine.getOptionValue("o");

        outputDir = new File(outputFile).getParent();
        Logger logger = initLog(outputDir);

        if (!commandLine.hasOption("bed")) {
            logger.error("Peak calling BED format file can not be empty");
            help.printHelp("java -jar renlabm6a_allele.jar AsmPeakDetection", header, options, footer, true);
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
            help.printHelp("java -jar renlabm6a_allele.jar AsmPeakDetection", header, options, footer, true);
            System.exit(2);
        } else {
            File vcf = new File(commandLine.getOptionValue("vcf"));
            if (!vcf.exists() || !vcf.isFile()) {
                logger.error("invalid file path: " + vcf.getAbsolutePath());
                System.exit(2);
            }
            aseVcfFile = vcf.getAbsolutePath();
        }

        if (!commandLine.hasOption("g")) {
            logger.error("GTF format file can not be empty");
            help.printHelp("java -jar renlabm6a_allele.jar AsmPeakDetection", header, options, footer, true);
            System.exit(2);
        } else {
            File gtf = new File(commandLine.getOptionValue("g"));
            if (!gtf.exists() || !gtf.isFile()) {
                logger.error("invalid file path: " + gtf.getAbsolutePath());
                System.exit(2);
            }
            gtfFile = gtf.getAbsolutePath();
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
        if (commandLine.hasOption("rc"))
            ipSNPCoverageInfimum = Integer.parseInt(commandLine.getOptionValue("rc"));
        if (commandLine.hasOption("bc"))
            wesSNPCoverageInfimum = Integer.parseInt(commandLine.getOptionValue("bc"));
        if (commandLine.hasOption("t")) {
            if (Integer.valueOf(commandLine.getOptionValue("t")) < 0) {
                System.err.println("invalid thread number, should be a positive integer");
                System.exit(2);
            }
            threadNumber = Integer.valueOf(commandLine.getOptionValue("t"));
        }

        AsmPeakDetection apd = new AsmPeakDetection(gtfFile, bedFile, aseVcfFile, wesVcfFile, dbsnpFile, outputFile,
                                                    degreeOfFreedom, ipSNPCoverageInfimum, wesSNPCoverageInfimum,
                                                    samplingTime, burn_in, threadNumber, logger);
        apd.getTestResult();
    }

    public HashMap<String, String> getPeakCoveredGenes() {
        return this.peakCoveredGene;
    }

    public HashMap<String, HashMap<String, HashMap<String, HashMap<String, Integer>>>> getPeakSnpReadsCount() {
        return this.peakSnpReadsCount;
    }

    public HashMap<String, HashSet<Integer>> getSnvForTest() {
        return this.snvForTest;
    }

    public void getTestResult() {
        this.getPeakCoveredGene();
        this.parseGTFFile();
        this.asmPeakTest();
        this.bhRecalibrationOfEachPeak();
        this.outputResult();
    }

    /**
     * locate m6A signal in BED format file to corresponding gene, peakCoveredGene = {chr:peakStart:peakEnd -> geneId, ...}
     */
    public void getPeakCoveredGene() {
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
     * get major allele and minor allele nucleotide and reads counts of each MeRIP-seq INPUT sample SNV sites covered by m6A peak
     * [chr1:peakStart:peakEnd -> [position1:majorNC, position2:majorNC],....]
     */
    public void getPeakSNPReadsCount() {
        HeterozygoteReadsCount hrc = new HeterozygoteReadsCount(this.peakCoveredSnpFile, this.ipSNPReadInfimum, this.log);
        this.peakSnpReadsCount = hrc.getMajorMinorHaplotype();
        this.peakMajorAlleleNucleotide = hrc.getPeakMajorAlleleNucleotides();
    }

    /**
     * get major allele and minor allele nucleotide and reads counts of each WES SNV sites covered by m6A peak
     *  [chr1: [peak1: position1: [major: count, minor: count], position2:[major: count, minor:count]], chr2:....]
     */
    private void getPeakSNPBackground() {
        if (this.wesFile != null) {
            HeterozygoteReadsCount hrc = new HeterozygoteReadsCount(this.peakCoveredWesSnpFile, this.wesSNPReadInfimum, this.log);
            this.peakSnpBackground = hrc.getMajorMinorHaplotype();
        } else
            this.peakSnpBackground = null;
    }

    /**
     * test the ASM significant p value of a m6A peak
     */
    public void asmPeakTest() {
        // [chr1: [peak1: position1: [major: count, minor: count], position2:[major: count, minor:count]], chr2:....]
        this.getPeakSNPReadsCount();
        this.getPeakSNPBackground();
        this.dataPreparation();
        this.calcLorStd();
        this.hierarchicalModelTest();
    }

    /**
     * prepare data for hierarchical test
     */
    public void dataPreparation() {
        HashMap<String, HashMap<String, HashMap<String, Integer>>> rnaSeqPeakSnvAlleleReads, wesPeakSnvAlleleReads;
        HashMap<String, HashMap<String, Integer>> rnaSeqMutPositionAlleleReads, wesMutPositionAlleleReads;
        HashMap<String, Integer> rnaSeqReads, wesReads;
        int[] majorCount, minorCount, majorBackground, minorBackground;
        int major, minor, backgroundMajor, backgroundMinor;
        String majorAllele, minorAllele, label;
        ArrayList<Integer> rnaSeqMajor, rnaSeqMinor, wesMajor, wesMinor;


        for (String chrNum : this.peakSnpReadsCount.keySet()) {
            // m6A signal on a particular chromosome
            rnaSeqPeakSnvAlleleReads = this.peakSnpReadsCount.get(chrNum);

            if (this.peakSnpBackground != null)
                wesPeakSnvAlleleReads = this.peakSnpBackground.getOrDefault(chrNum, null);
            else
                wesPeakSnvAlleleReads = null;

            for (String peakRange : rnaSeqPeakSnvAlleleReads.keySet()) {
                // SNV sites covered by the m6A signal
                rnaSeqMutPositionAlleleReads = rnaSeqPeakSnvAlleleReads.get(peakRange);

                if (wesPeakSnvAlleleReads != null)
                    wesMutPositionAlleleReads = wesPeakSnvAlleleReads.getOrDefault(peakRange, null);
                else
                    wesMutPositionAlleleReads = null;

                rnaSeqMajor = new ArrayList<>();
                rnaSeqMinor = new ArrayList<>();
                wesMajor = new ArrayList<>();
                wesMinor = new ArrayList<>();
                for (String position : rnaSeqMutPositionAlleleReads.keySet()) {
                    // filter SNV sites which
                    if (this.dbsnpRecord != null && !DbsnpAnnotation.getSearchRes(this.dbsnpRecord, chrNum, position))
                        continue;

                    rnaSeqReads = rnaSeqMutPositionAlleleReads.get(position);
                    if (wesMutPositionAlleleReads != null)
                        wesReads = wesMutPositionAlleleReads.getOrDefault(position, null);
                    else
                        wesReads = null;

                    // sorted by reads counts
                    List<Map.Entry<String, Integer>> nucleotides = new ArrayList<>(rnaSeqReads.entrySet());
                    Collections.sort(nucleotides, new Comparator<Map.Entry<String, Integer>>() {
                        @Override
                        public int compare(Map.Entry<String, Integer> o1, Map.Entry<String, Integer> o2) {
                            return o2.getValue().compareTo(o1.getValue());
                        }
                    });

                    ArrayList<String> ncs = new ArrayList<>();
                    ArrayList<Integer> counts = new ArrayList<>();
                    for (Map.Entry<String, Integer> entry : nucleotides) {
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
                        backgroundMajor = (minor + major) / 2;
                        backgroundMinor = (minor + major) / 2;
                    }
                    wesMajor.add(backgroundMajor);
                    wesMinor.add(backgroundMinor);
                    HashSet<Integer> snvs = this.snvForTest.getOrDefault(chrNum, new HashSet<>());
                    snvs.add(Integer.valueOf(position));
                    this.snvForTest.put(chrNum, snvs);
                    ncs = null;
                    counts = null;
                }
                // with no mutation site record in dbsnp
                if (rnaSeqMajor.size() == 0)
                    continue;
                majorCount = new int[rnaSeqMajor.size()];
                minorCount = new int[rnaSeqMinor.size()];
                majorBackground = new int[wesMajor.size()];
                minorBackground = new int[wesMinor.size()];
                assert majorCount.length == minorCount.length;
                assert majorCount.length == majorBackground.length;
                assert majorBackground.length == minorBackground.length;
                for (int i = 0; i < majorCount.length; i++) {
                    majorCount[i] = rnaSeqMajor.get(i);
                    minorCount[i] = rnaSeqMinor.get(i);
                    majorBackground[i] = wesMajor.get(i);
                    minorBackground[i] = wesMinor.get(i);
                }

                label = String.join(":", new String[]{chrNum, peakRange});

                this.peakSNVNum.put(label, majorCount.length);

                ArrayList<int[]> statistic = new ArrayList<>(4);
                statistic.add(majorCount);
                statistic.add(minorCount);
                statistic.add(majorBackground);
                statistic.add(minorBackground);
                this.statisticForTest.put(label, statistic);

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

                rnaSeqMajor = null;
                rnaSeqMinor = null;
                wesMajor = null;
                wesMinor = null;
            }
        }
        this.dbsnpRecord = null;
        if (this.statisticForTest.isEmpty()) {
            this.log.error("contains no peaks with SNV sites for hierarchical test, please check the input data");
            System.exit(2);
        }
    }

    /**
     * calculate the standard deviation of LOR of all SNV sites covered by m6A signals on genome
     */
    private void calcLorStd() {
        this.log.debug("calculate LOR standard as parameter of Inverse chi-square distribution");
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
                    minor = 0.5;

                // WES SNVs only used for annotating RNA-seq SNVs. Because of the alignment error, background reads count do not take part in LOR Std calculation
                // we assume the odd ratio of background major and minor reads equals 1
                lor = (major / minor);  // (majorBack / minorBack);
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
        lorList.clear();

        // avoid org.apache.commons.math3.exception.NotStrictlyPositiveException: standard deviation (0)
        this.lorStd = (Math.abs(this.lorStd-0)<0.00001)? 0.0001: this.lorStd;
    }

    /**
     * calculate m6A signal p value via hierarchical model
     */
    private void hierarchicalModelTest() {
        this.log.debug("hierarchical Bayesian model test start");
        ExecutorService threadPoolExecutor = Executors.newFixedThreadPool(this.threadNumber);
        this.lock = new ReentrantLock();
        CountDownLatch countDown = new CountDownLatch(this.statisticForTest.size());

        RunTest task = (String name) -> {
          return new Runnable() {
              @Override
              public void run() {
                  ArrayList<int[]> statistic = statisticForTest.get(name);
                  int[] majorCount = statistic.get(0);
                  int[] minorCount = statistic.get(1);
                  int[] majorBackground = statistic.get(2);
                  int[] minorBackground = statistic.get(3);
                  // get p value via hierarchical model
                  HierarchicalBayesianModel hb = new HierarchicalBayesianModel(lorStd, degreeOfFreedom, samplingTime,
                          burnIn, majorCount, minorCount, majorBackground, minorBackground);
                  double pVal = hb.testSignificant();
                  double peakOddRatio = Math.exp(hb.quantifyGeneLOR());
                  double peakMAF = Math.min(1.0, peakOddRatio / (peakOddRatio+1));

                  lock.lock();
                  try {
                      ArrayList<String> samePValPeaks = asmPValue.getOrDefault(pVal, new ArrayList<>());
                      samePValPeaks.add(name);
                      asmPValue.put(pVal, samePValPeaks);
                      peakMajorAlleleFrequency.put(name, peakMAF);
                      countDown.countDown();
                  } finally {
                      lock.unlock();
                      majorCount = null;
                      minorCount = null;
                      majorBackground = null;
                      minorBackground = null;
                      hb = null;
                  }
              }
          };
        };

        for (String str: this.statisticForTest.keySet()) {
            this.totalPeakCount++;
            Runnable runnable = task.runTask(str);
            threadPoolExecutor.submit(runnable);
        }
        try {
            countDown.await();
        } catch (InterruptedException ie) {
            this.log.error("analysis interrupted");
        } finally {
            this.statisticForTest = null;
            threadPoolExecutor.shutdown();
            this.lock = null;
        }
        this.log.debug("model test complete");
    }

    /**
     * recalibrate the p value with BH method, get significant q value
     */
    private void bhRecalibrationOfEachPeak() {
        this.log.debug("start recalibrating p values of hierarchical model");
        this.log.debug("sorting test result in order");
        ArrayList<Map.Entry<Double, ArrayList<String>>> sortedPVals = new ArrayList<>(this.asmPValue.entrySet());
        // sort p value from large to small
        Collections.sort(sortedPVals, new Comparator<Map.Entry<Double, ArrayList<String>>>() {
            @Override
            public int compare(Map.Entry<Double, ArrayList<String>> o1, Map.Entry<Double, ArrayList<String>> o2) {
                return o2.getKey().compareTo(o1.getKey());
            }
        });

        int totalPeak = this.totalPeakCount, rankage = totalPeak;
        double prevQValue = 1.0, qValue;
        String pValString, qValString;
        for (Map.Entry<Double, ArrayList<String>> entry: sortedPVals) {
            Double pVal = entry.getKey();
            ArrayList<String> samePValPeaks = entry.getValue();
            // sort items with its SNV number when p value is same
            HashMap<String, Integer> samePValPeaksSNVs = new HashMap<>();
            for (String peak: samePValPeaks)
                samePValPeaksSNVs.put(peak, this.peakSNVNum.get(peak));

            // sort items with its major allele frequency when p value and SNV numbers are same
            HashMap<String, Double> samePValPeakMajorAlleleFrequency = new HashMap<>();
            for (String peak: samePValPeaks)
                samePValPeakMajorAlleleFrequency.put(peak, this.peakMajorAlleleFrequency.get(peak));

            List<Map.Entry<String, Integer>> samePValPeakEntry = new ArrayList<>(samePValPeaksSNVs.entrySet());
            Collections.sort(samePValPeakEntry, new Comparator<Map.Entry<String, Integer>>() {
                @Override
                public int compare(Map.Entry<String, Integer> o1, Map.Entry<String, Integer> o2) {

                    String peak1 = o1.getKey(), peak2 = o2.getKey();
                    Integer peak1SNVs = o1.getValue(), peak2SNVs = o2.getValue();
                    Double peak1MAF = samePValPeakMajorAlleleFrequency.get(peak1), peak2MAF = samePValPeakMajorAlleleFrequency.get(peak2);
                    if (peak1SNVs.equals(peak2SNVs)) {
                        return peak2MAF.compareTo(peak1MAF);
                    } else
                        return peak2SNVs.compareTo(peak1SNVs);
                }
            });

            for (Map.Entry<String, Integer> geneEntry: samePValPeakEntry) {
                String peak = geneEntry.getKey();
                qValue = Math.min(prevQValue, pVal * totalPeak / rankage);
                if (qValue - prevQValue < 0.00001)
                    prevQValue = qValue;
                rankage--;

                pValString = Double.toString(pVal);
                qValString = this.df.format(qValue);
                this.asmQValue.add(String.join("->", new String[]{peak, pValString, qValString}));
            }
        }
        this.log.debug("recalibration complete.");
    }

    /**
     * write the test result into file
     */
    private void outputResult() {
        ArrayList<String[]> outputRecord = new ArrayList<>();
        BufferedWriter bfw = null;
        try {
            bfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(this.asmPeakFile))));
            String line, label, geneId, geneName, chrNum, peakStart, peakEnd, majorAlleleReads, minorAlleleReads,
                   majorAlleleBackground, minorAlleleBackground, pValue, qValue;
            LinkedList<String> majorAlleleNucleotides;
            String[] info, rec, finalInfo;
            String majorAlleleFrequency;
            int snvNum;
            bfw.write("#chr\tpeakStart\tpeakEnd\tgeneId\tgeneName\tp-value\tq-value\tsnvNum\tmajor/minorAlleleReads\tmajor/minorBackground\tMajorAlleleFrequency\tmajorAlleleNC\n");
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
                geneName = this.geneNames.get(chrNum).getOrDefault(geneId, "unknown");
                majorAlleleFrequency = this.df.format(this.peakMajorAlleleFrequency.get(label));
                majorAlleleNucleotides = this.peakMajorAlleleNucleotide.get(label);
                String majorAlleleRecords = this.getString(majorAlleleNucleotides);
                snvNum = this.peakSNVNum.get(label);
                majorAlleleReads = String.valueOf(this.peakMajorMinorAlleleCount.get(label).get("major"));
                minorAlleleReads = String.valueOf(this.peakMajorMinorAlleleCount.get(label).get("minor"));
                majorAlleleBackground = String.valueOf(this.peakMajorMinorBackground.get(label).get("major"));
                minorAlleleBackground = String.valueOf(this.peakMajorMinorBackground.get(label).get("minor"));
                finalInfo = new String[]{chrNum, peakStart, peakEnd, geneId, geneName, pValue, qValue,
                                         Integer.toString(snvNum), majorAlleleReads + "," + minorAlleleReads,
                                         majorAlleleBackground+","+minorAlleleBackground, majorAlleleFrequency,
                                         majorAlleleRecords};
                outputRecord.add(finalInfo);
            }

            Collections.sort(outputRecord, new Comparator<String[]>() {
                @Override
                public int compare(String[] o1, String[] o2) {
                    Double q1 = Double.parseDouble(o1[6]), q2 = Double.parseDouble(o2[6]);
                    if (!q1.equals(q2))
                        return q1.compareTo(q2);
                    // sort records with same q-value after BH recalibration by SNV number
                    Integer snvCount1 = Integer.parseInt(o1[7]), snvCount2 = Integer.parseInt(o2[7]);
                    return snvCount2.compareTo(snvCount1);
                }
            });
            for (String[] record: outputRecord) {
                line = String.join("\t", record);
                bfw.write(line);
                bfw.newLine();
            }
            this.log.debug("result file " + this.asmPeakFile);
        } catch (IOException ie) {
            this.log.error(ie.getMessage());
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

    public void parseGTFFile() {
        BufferedReader bfr = null;
        this.geneNames = new HashMap<>();
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(this.gtfFile)));
            String line = "", chrNum, geneId, geneName;
            String[] info, geneInfo;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    if (!info[2].equals("gene")) {
                        info = null;
                        continue;
                    }
                    chrNum = info[0];
                    geneInfo = this.getGeneInfo(info[8]);
                    geneId = geneInfo[0];
                    geneName = geneInfo[1];
                    HashMap<String, String> chrGenes = this.geneNames.getOrDefault(chrNum, new HashMap<>());
                    chrGenes.put(geneId, geneName);
                    this.geneNames.put(chrNum, chrGenes);
                }
                info = null;
            }
        } catch (IOException ie) {
            ie.printStackTrace();
            System.exit(2);
        } finally {
            if (bfr != null) {
                try {
                    bfr.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    private String[] getGeneInfo(String recordInfo) {
        String[] info = recordInfo.split("; ");
        String geneName = null, geneId = null;
        for (String s: info) {
            if (s.startsWith("gene_id")) {
                String[] name = s.split(" ");
                geneId = name[1].substring(1, name[1].length() -1);
            }
            if (s.startsWith("gene_name")) {
                String[] name = s.split(" ");
                geneName = name[1].substring(1, name[1].length() -1);
            }
        }

        return new String[] {geneId, geneName};
    }

    /**
     * initial log4j Logger instance
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

        option = new Option("g", "gtf", true, "GTF format file");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("wes", "wes_vcf_file", true, "WES SNP calling VCF format file");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("bed", "peak_bed_file", true, "Peak calling output result in BED format");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("db", "dbsnp", true, "big scale SNV annotation data set, like dbsnp, 1000Genome etc. (the file format see https://github.com/Jakob666/allele-specificM6A)");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("o", "output", true, "ASM m6A signal test output file, default ./asmPeak.txt");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("df", "degree_of_freedom", true, "degree of freedom of inverse-Chi-square distribution, default 10");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("rc", "reads_coverage", true, "reads coverage threshold using for filter RNA-seq or MeRIP-seq data SNVs in VCF file (aim for reducing FDR), default 10");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("bc", "bkg_coverage", true, "reads coverage threshold using for filter WES data SNVs in VCF file (aim for reducing FDR), default 30");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("s", "sampling", true, "sampling times, larger than 500, default 10000");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("b", "burn", true, "burn-in times, more than 100 and less than sampling times. Default 2000");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("t", "thread", true, "thread number for running test. Default 2");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("h", "help", false, "help message of AsmPeakDetection");
        option.setRequired(false);
        options.addOption(option);

        CommandLineParser parser = new DefaultParser();

        return parser.parse(options, args);
    }
}
