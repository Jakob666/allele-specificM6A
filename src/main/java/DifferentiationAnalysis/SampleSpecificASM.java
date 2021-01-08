package DifferentiationAnalysis;

import HierarchicalBayesianAnalysis.AsmPeakDetection;
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
import java.util.stream.Collectors;

public class SampleSpecificASM {
    private Logger logger;
    private String outputFile;
    private int samplingTime, burnIn, threadNumber;
    private AsmPeakDetection sample1AsmPeakDetector, sample2AsmPeakDetector;
    private HashMap<String, String> sample1PeakCoveredGene, sample2PeakCoveredGene;
    private HashMap<String, HashSet<String>> mergedPeakCoveredGenes;
    private HashMap<String, Double> peakMajorAlleleFrequency;
    private HashMap<String, HashMap<String, HashMap<String, HashMap<String, Integer>>>> sample1PeakSnpReadsCount, sample2PeakSnpReadsCount;
    private HashMap<String, HashMap<String, String>> geneNames;
    private HashMap<String, HashSet<Integer>> sample1SNVs, sample2SNVs;
    private HashMap<String, HashSet<Integer>> commonSNVs;
    private HashMap<String, Integer> peakSNVNum = new HashMap<>();
    private HashMap<String, ArrayList<int[]>> chrMergedPeaks;
    private HashMap<String, ArrayList<Integer>> sample1Major, sample1Minor, sample2Major, sample2Minor;
    private HashMap<Double, ArrayList<String>> peakSampleSpecificPValue;
    private HashMap<String, int[]> peakMajorMinorReads;
    private ArrayList<String> peakSampleSpecificQValue;
    private DecimalFormat df = new DecimalFormat("0.0000");
    private ReentrantLock lock;

    /**
     * Constructor
     * @param gtfFile GTF annotation file
     * @param sample1PeakBedFile sample1 BED format file via MeRIP-seq IP data
     * @param sample1VcfFile sample1 VCF format file via MeRIP-seq INPUT data
     * @param sample1WesFile sample1 VCF format file via WES data, optional
     * @param sample2PeakBedFile sample2 BED format file via MeRIP-seq IP data
     * @param sample2VcfFile sample2 VCF format file via MeRIP-seq INPUT data
     * @param sample2WesFile sample2 VCF format file via WES data, optional
     * @param dbsnpFile dbsnp file for SNP filtering
     * @param outputFile test result output file
     * @param ipSNPReadInfimum reads coverage threshold when filter INPUT sample SNV sites, default 10
     * @param wesSNPReadInfimum reads coverage threshold when filter WES SNV sites, default 30
     * @param samplingTime sampling time, default 10000
     * @param burnIn burn in time, default 2000
     * @param threadNumber thread number, default 2
     * @param logger log4j instance
     */
    public SampleSpecificASM(String gtfFile, String sample1PeakBedFile, String sample1VcfFile, String sample1WesFile,
                             String sample2PeakBedFile, String sample2VcfFile, String sample2WesFile, String dbsnpFile,
                             String outputFile, int ipSNPReadInfimum, int wesSNPReadInfimum,
                             int samplingTime, int burnIn, int threadNumber, Logger logger) {
        this.logger = logger;
        this.outputFile = outputFile;
        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
        this.threadNumber = threadNumber;

        String peakWithSnvFile = new File(new File(outputFile).getParent(), "sample1_peak_with_snv.txt").getAbsolutePath();
        String peakWithSnvBkgFile = null;
        if (sample1WesFile!=null)
            peakWithSnvBkgFile = new File(new File(outputFile).getParent(), "sample1_peak_with_snv_bkg.txt").getAbsolutePath();
        this.sample1AsmPeakDetector = new AsmPeakDetection(gtfFile, sample1PeakBedFile, sample1VcfFile, sample1WesFile,
                                                           dbsnpFile, peakWithSnvFile, peakWithSnvBkgFile, outputFile,
                                                           ipSNPReadInfimum, wesSNPReadInfimum, samplingTime,
                                                           burnIn, threadNumber, logger);

        peakWithSnvFile = new File(new File(outputFile).getParent(), "sample2_peak_with_snv.txt").getAbsolutePath();
        peakWithSnvBkgFile = new File(new File(outputFile).getParent(), "sample2_peak_with_snv_bkg.txt").getAbsolutePath();
        this.sample2AsmPeakDetector = new AsmPeakDetection(gtfFile, sample2PeakBedFile, sample2VcfFile, sample2WesFile,
                                                           dbsnpFile, peakWithSnvFile, peakWithSnvBkgFile, outputFile,
                                                           ipSNPReadInfimum, wesSNPReadInfimum, samplingTime,
                                                           burnIn, threadNumber, logger);
        this.parseGTFFile(gtfFile);
    }

    public static void main(String[] args) {
        Options options = new Options();
        CommandLine commandLine = null;
        HelpFormatter help = new HelpFormatter();
        String header = "SampleSpecificASM contains following parameters: ";
        String footer = "";

        try {
            commandLine = setCommandLine(args, options);
        } catch (ParseException pe) {
            System.err.println(pe.getMessage());
            help.printHelp("java -jar renlabm6a_allele.jar SampleSpecificASM", header, options, footer, true);
            System.exit(2);
        }

        if (commandLine.hasOption("h")) {
            help.printHelp("java -jar renlabm6a_allele.jar SampleSpecificASM", header, options, footer, true);
            System.exit(0);
        }

        // default parameters
        String gtfFile = null, sample1BedFile = null, sample2BedFile = null, sample1AseVcfFile = null, sample2AseVcfFile = null,
               sample1WesVcfFile = null, sample2WesVcfFile = null, dbsnpFile = null, outputFile, outputDir;
        int ipSNPCoverageInfimum = 10, wesSNPCoverageInfimum = 30, samplingTime = 50000, burn_in = 10000, threadNumber = 2;

        if (!commandLine.hasOption("o"))
            outputFile = new File(System.getProperty("user.dir"), "sampleSpecificASM.txt").getAbsolutePath();
        else
            outputFile = commandLine.getOptionValue("o");

        outputDir = new File(outputFile).getParent();
        Logger logger = initLog(outputDir);

        if (!commandLine.hasOption("s1_bed")) {
            logger.error("sample1 Peak calling BED format file can not be empty");
            help.printHelp("java -jar renlabm6a_allele.jar SampleSpecificASM", header, options, footer, true);
            System.exit(2);
        } else {
            File bed = new File(commandLine.getOptionValue("s1_bed"));
            if (!bed.exists() || !bed.isFile()) {
                logger.error("invalid BED format file path: " + bed.getAbsolutePath());
                System.exit(2);
            }
            sample1BedFile = bed.getAbsolutePath();
        }

        if (!commandLine.hasOption("s2_bed")) {
            logger.error("sample2 Peak calling BED format file can not be empty");
            help.printHelp("java -jar renlabm6a_allele.jar SampleSpecificASM", header, options, footer, true);
            System.exit(2);
        } else {
            File bed = new File(commandLine.getOptionValue("s2_bed"));
            if (!bed.exists() || !bed.isFile()) {
                logger.error("invalid BED format file path: " + bed.getAbsolutePath());
                System.exit(2);
            }
            sample2BedFile = bed.getAbsolutePath();
        }

        if (!commandLine.hasOption("s1_vcf")) {
            logger.error("sample1 SNP calling VCF format file can not be empty");
            help.printHelp("java -jar renlabm6a_allele.jar SampleSpecificASM", header, options, footer, true);
            System.exit(2);
        } else {
            File vcf = new File(commandLine.getOptionValue("s1_vcf"));
            if (!vcf.exists() || !vcf.isFile()) {
                logger.error("invalid file path: " + vcf.getAbsolutePath());
                System.exit(2);
            }
            sample1AseVcfFile = vcf.getAbsolutePath();
        }

        if (!commandLine.hasOption("s2_vcf")) {
            logger.error("sample2 SNP calling VCF format file can not be empty");
            help.printHelp("java -jar renlabm6a_allele.jar SampleSpecificASM", header, options, footer, true);
            System.exit(2);
        } else {
            File vcf = new File(commandLine.getOptionValue("s2_vcf"));
            if (!vcf.exists() || !vcf.isFile()) {
                logger.error("invalid file path: " + vcf.getAbsolutePath());
                System.exit(2);
            }
            sample2AseVcfFile = vcf.getAbsolutePath();
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

        if (!commandLine.hasOption("g")) {
            logger.error("GTF format file can not be empty");
            help.printHelp("java -jar renlabm6a_allele.jar SampleSpecificASM", header, options, footer, true);
            System.exit(2);
        } else {
            File gtf = new File(commandLine.getOptionValue("g"));
            if (!gtf.exists() || !gtf.isFile()) {
                logger.error("invalid file path: " + gtf.getAbsolutePath());
                System.exit(2);
            }
            gtfFile = gtf.getAbsolutePath();
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

        SampleSpecificASM ssasm = new SampleSpecificASM(gtfFile, sample1BedFile, sample1AseVcfFile, sample1WesVcfFile,
                                                        sample2BedFile, sample2AseVcfFile, sample2WesVcfFile, dbsnpFile,
                                                        outputFile, ipSNPCoverageInfimum, wesSNPCoverageInfimum,
                                                        samplingTime, burn_in, threadNumber, logger);
        ssasm.testSampleSpecificAsm();
    }

    public void testSampleSpecificAsm() {
        this.dataPreparation();
        this.asmDifferentiationAnalysis();
        this.bhRecalibrationOfEachPeak();
        this.output();
    }

    /**
     * sample specific ASE detection
     */
    private void dataPreparation() {
        // prepare data for test
        this.sample1AsmPeakDetector.getPeakSNPReadsCount();
        this.sample1AsmPeakDetector.dataPreparation();

        this.sample2AsmPeakDetector.getPeakSNPReadsCount();
        this.sample2AsmPeakDetector.dataPreparation();
        // whether m6A peaks cover common genes between two samples
        // chr:peakStart:peakEnd -> geneId
        this.sample1PeakCoveredGene = this.sample1AsmPeakDetector.getPeakCoveredGenes();
        this.sample2PeakCoveredGene = this.sample2AsmPeakDetector.getPeakCoveredGenes();
        // chr -> [pos1, pos2, ...]
        this.sample1SNVs = this.sample1AsmPeakDetector.getSnvForTest();
        this.sample2SNVs = this.sample2AsmPeakDetector.getSnvForTest();
        this.logger.debug("extract common peak covered SNVs between two samples");
        // chr-> {pos1, pos2,...}
        this.getCommonSNVs();
        if (this.commonSNVs.isEmpty()) {
            this.logger.debug("detect no common peak covered SNVs for test");
            System.exit(0);
        }
        this.sample1PeakSnpReadsCount = this.sample1AsmPeakDetector.getPeakSnpReadsCount();
        this.sample2PeakSnpReadsCount = this.sample2AsmPeakDetector.getPeakSnpReadsCount();
        HashSet<String> sample1Genes = new HashSet<>(), sample2Genes = new HashSet<>();
        sample1PeakCoveredGene.entrySet().forEach(entry -> sample1Genes.add(entry.getValue()));
        sample2PeakCoveredGene.entrySet().forEach(entry -> sample2Genes.add(entry.getValue()));
        sample1Genes.retainAll(sample2Genes);
        if (sample1Genes.isEmpty()) {
            this.logger.debug("detect no common peak covered genes for allele-specific modification differentiation analysis between two samples");
            System.exit(0);
        } else
            this.logger.debug("start data preparation");

        // merging m6A signal peaks on common genes between two samples
        this.sample1Major = new HashMap<>();
        this.sample1Minor = new HashMap<>();
        this.sample2Major = new HashMap<>();
        this.sample2Minor = new HashMap<>();
        this.mergedPeakCoveredGenes = new HashMap<>();
        this.mergingM6ASignalPeaks(sample1Genes);
    }

    /**
     * common peak covered SNVs between two samples
     */
    private void getCommonSNVs() {
        this.commonSNVs = new HashMap<>();
        HashSet<String> sample1Chrs = new HashSet<>(this.sample1SNVs.keySet());
        HashSet<String> sample2Chrs = new HashSet<>(this.sample2SNVs.keySet());
        sample1Chrs.retainAll(sample2Chrs);
        for (String chr: sample1Chrs) {
            HashSet<Integer> sites1 = this.sample1SNVs.get(chr);
            HashSet<Integer> sites2 = this.sample2SNVs.get(chr);
            sites1.retainAll(sites2);
            if (!sites1.isEmpty())
                this.commonSNVs.put(chr, sites1);
        }
    }

    /**
     * merge overlapped peaks no chromosomes
     * @param commonGenes common genes covered by peaks in two samples
     */
    private void mergingM6ASignalPeaks(HashSet<String> commonGenes) {
        this.logger.debug("merging m6A signal peaks between two samples");
        // chr1-> [peakStart:peakEnd-> position-> [majorNC-> count, minorNC-> count]
        HashMap<String, HashMap<String, HashMap<String, HashMap<String, Integer>>>> sample1PeakCoveredSNV = this.sample1AsmPeakDetector.getPeakSnpReadsCount();
        HashMap<String, HashMap<String, HashMap<String, HashMap<String, Integer>>>> sample2PeakCoveredSNV = this.sample2AsmPeakDetector.getPeakSnpReadsCount();
        // common chromosome
        HashSet<String> sample1Chrs = new HashSet<>(sample1PeakCoveredSNV.keySet());
        HashSet<String> sample2Chrs = new HashSet<>(sample2PeakCoveredSNV.keySet());
        sample1Chrs.retainAll(sample2Chrs);
        if (sample1Chrs.isEmpty()) {
            this.logger.debug("peaks in two samples cover no common chromosome");
            System.exit(0);
        }

        // sample1 and sample2 peaks on common covered genes
        this.chrMergedPeaks = new HashMap<>();
        HashMap<String, HashMap<String, HashMap<String, Integer>>> chrRecords;
        ArrayList<int[]> sample1Peaks, sample2Peaks;
        String label;
        String[] startEnd;
        for (String chrNum: sample1Chrs) {
            sample1Peaks = new ArrayList<>();
            sample2Peaks = new ArrayList<>();
            chrRecords = sample1PeakCoveredSNV.get(chrNum);
            for (String peakRange: chrRecords.keySet()) {
                label = chrNum + ":" + peakRange;
                if (!commonGenes.contains(this.sample1PeakCoveredGene.get(label)))
                    continue;
                startEnd = peakRange.split(":");
                sample1Peaks.add(new int[] {Integer.valueOf(startEnd[0]), Integer.valueOf(startEnd[1])});
            }
            sample1Peaks = (ArrayList<int[]>) sample1Peaks.stream().sorted((o1, o2) -> o1[0]-o2[0]).collect(Collectors.toList());   // sort from small to large
            chrRecords = sample2PeakCoveredSNV.get(chrNum);
            for (String peakRange: chrRecords.keySet()) {
                label = chrNum + ":" + peakRange;
                if (!commonGenes.contains(this.sample2PeakCoveredGene.get(label)))
                    continue;
                startEnd = peakRange.split(":");
                sample2Peaks.add(new int[] {Integer.valueOf(startEnd[0]), Integer.valueOf(startEnd[1])});
            }
            sample2Peaks = (ArrayList<int[]>) sample2Peaks.stream().sorted((o1, o2) -> o1[0]-o2[0]).collect(Collectors.toList());   // sort from small to large

            // merge overlapped m6A peaks between two samples
            ArrayList<int[]> mergedPeaks = new ArrayList<>();
            ArrayList<ArrayList<String>> sample1Component = new ArrayList<>(), sample2Component = new ArrayList<>();
            boolean overlap, tailOverlap;
            Iterator<int[]> sample1PeakItr = sample1Peaks.iterator(), sample2PeakItr = sample2Peaks.iterator();
            int[] interval1 = null, interval2 = null;
            while (sample1PeakItr.hasNext() || sample2PeakItr.hasNext()) {
                if (interval1 == null && sample1PeakItr.hasNext())
                    interval1 = sample1PeakItr.next();
                if (interval2 == null && sample2PeakItr.hasNext())
                    interval2 = sample2PeakItr.next();

                if (interval1 != null && interval2 != null) {   // test whether the two intervals are overlapped
                    int start1 = interval1[0], end1 = interval1[1], start2 = interval2[0], end2 = interval2[1];
                    overlap = this.ifOverlapped(start1, end1, start2, end2);
                    if (overlap) {  // two interval overlapped
                        int[] newRange = new int[]{Math.min(start1, start2), Math.max(end1, end2)};
                        if (mergedPeaks.isEmpty()) {
                            ArrayList<String> component1 = new ArrayList<>();
                            component1.add(start1 + ":" + end1);
                            sample1Component.add(component1);
                            ArrayList<String> component2 = new ArrayList<>();
                            component2.add(start2 + ":" + end2);
                            sample2Component.add(component2);
                            mergedPeaks.add(newRange);
                        } else {
                            int[] lastPeak = mergedPeaks.get(mergedPeaks.size() - 1);
                            tailOverlap = this.ifOverlapped(lastPeak[0], lastPeak[1], newRange[0], newRange[1]);
                            if (tailOverlap) {
                                ArrayList<String> component1 = sample1Component.get(sample1Component.size() - 1);
                                component1.add(start1 + ":" + end1);
                                sample1Component.remove(sample1Component.size() - 1);
                                sample1Component.add(component1);
                                ArrayList<String> component2 = sample2Component.get(sample2Component.size() - 1);
                                component2.add(start2 + ":" + end2);
                                sample2Component.remove(sample2Component.size() - 1);
                                sample2Component.add(component2);
                                mergedPeaks.remove(mergedPeaks.size() - 1);
                                lastPeak = new int[]{Math.min(lastPeak[0], newRange[0]), Math.max(lastPeak[1], newRange[1])};
                                mergedPeaks.add(lastPeak);
                            } else {
                                ArrayList<String> component1 = new ArrayList<>();
                                component1.add(start1 + ":" + end1);
                                sample1Component.add(component1);
                                ArrayList<String> component2 = new ArrayList<>();
                                component2.add(start2 + ":" + end2);
                                sample2Component.add(component2);
                                mergedPeaks.add(newRange);
                            }
                        }
                        interval1 = null;
                        interval2 = null;
                    } else {    // two independent interval
                        int[] front = (interval1[0] <= interval2[0])? interval1 : interval2;
                        boolean interval1Front = interval1[0] <= interval2[0];
                        int[] lastPeak = (mergedPeaks.isEmpty()) ? null : mergedPeaks.get(mergedPeaks.size() - 1);
                        tailOverlap = (lastPeak != null) && this.ifOverlapped(lastPeak[0], lastPeak[1], front[0], front[1]);
                        if (tailOverlap) {
                            if (interval1Front) {
                                ArrayList<String> component1 = sample1Component.get(sample1Component.size() - 1);
                                component1.add(start1 + ":" + end1);
                                sample1Component.remove(sample1Component.size() - 1);
                                sample1Component.add(component1);
                            } else {
                                ArrayList<String> component2 = sample2Component.get(sample2Component.size() - 1);
                                component2.add(start2 + ":" + end2);
                                sample2Component.remove(sample2Component.size() - 1);
                                sample2Component.add(component2);
                            }
                            mergedPeaks.remove(mergedPeaks.size() - 1);
                            lastPeak = new int[]{Math.min(lastPeak[0], front[0]), Math.max(lastPeak[1], front[1])};
                            mergedPeaks.add(lastPeak);
                        }
                        if (interval1Front)
                            interval1 = null;
                        else
                            interval2 = null;
                    }
                } else if (interval1 != null || interval2 != null){
                    int[] front = (interval1!=null)? interval1: interval2;
                    boolean interval1Front = interval1!=null;
                    int[] lastPeak = (mergedPeaks.isEmpty()) ? null : mergedPeaks.get(mergedPeaks.size() - 1);
                    tailOverlap = (lastPeak != null) && this.ifOverlapped(lastPeak[0], lastPeak[1], front[0], front[1]);
                    if (tailOverlap) {
                        if (interval1Front) {
                            ArrayList<String> component1 = sample1Component.get(sample1Component.size() - 1);
                            component1.add(front[0] + ":" + front[1]);
                            sample1Component.remove(sample1Component.size() - 1);
                            sample1Component.add(component1);
                        } else {
                            ArrayList<String> component2 = sample2Component.get(sample2Component.size() - 1);
                            component2.add(front[0] + ":" + front[1]);
                            sample2Component.remove(sample2Component.size() - 1);
                            sample2Component.add(component2);
                        }
                        mergedPeaks.remove(mergedPeaks.size() - 1);
                        lastPeak = new int[]{Math.min(lastPeak[0], front[0]), Math.max(lastPeak[1], front[1])};
                        mergedPeaks.add(lastPeak);
                    }
                    if (interval1Front)
                        interval1 = null;
                    else
                        interval2 = null;
                } else
                    break;
            }
            if (!mergedPeaks.isEmpty()) {
                // common SNVs of two samples under merged peaks
                this.commonSNVCoveredByMergedPeaks(chrNum, mergedPeaks, sample1Component, sample2Component);
            }

        }
        if (chrMergedPeaks.isEmpty()) {
            this.logger.debug("contains no overlapped peaks in two samples, analysis interrupt");
            System.exit(0);
        } else
            this.logger.debug("peak merging completed");
    }

    /**
     * whether the two ranges are overlap
     * @param start1 start point of first region
     * @param end1 end point of first region
     * @param start2 start point of second region
     * @param end2 end point of second region
     * @return true, if is overlapped; otherwise, false
     */
    private boolean ifOverlapped(int start1, int end1, int start2, int end2) {
        //   s1 ------- e1       s1 -------- e1         s1 --------- e1     s1 ---- e2
        // s2 ------- e2             s2------- e2           s2 --- e2     s2 --------- e2
        if (Math.max(start1, start2) <= Math.min(end1, end2))
            return true;
        else
            return !(start1 > end2 || end1 < start2);
    }

    /**
     * common SNVs of two samples under merged peaks, record its reads count
     * @param chrNum chromosome number
     * @param mergedPeaks merged peaks
     * @param sample1Component merged peak sample1 sub-peaks
     * @param sample2Component merged peak sample2 sub-peaks
     */
    private void commonSNVCoveredByMergedPeaks(String chrNum, ArrayList<int[]> mergedPeaks,
                                               ArrayList<ArrayList<String>> sample1Component,
                                               ArrayList<ArrayList<String>> sample2Component) {
        HashSet<Integer> commonSites = this.commonSNVs.get(chrNum);
        if (commonSites == null)
            return;
        this.chrMergedPeaks.put(chrNum, mergedPeaks);
        // position-> [majorNC-> count, minorNC-> count]
        HashMap<String, HashMap<String, Integer>> sample1Reads, sample2Reads;
        for (int i=0; i<mergedPeaks.size(); i++) {
            int[] mergePeak = mergedPeaks.get(i);
            ArrayList<String> component1 = sample1Component.get(i);
            String label = chrNum+":"+mergePeak[0]+":"+mergePeak[1];
            for (String peakRange: component1) {
                String geneName = this.sample1PeakCoveredGene.get(chrNum+":"+peakRange);
                HashSet<String> genes = this.mergedPeakCoveredGenes.getOrDefault(label, new HashSet<>());
                genes.add(geneName);
                this.mergedPeakCoveredGenes.put(label, genes);
                sample1Reads = this.sample1PeakSnpReadsCount.get(chrNum).get(peakRange);
                List<Integer> mutPositions = sample1Reads.keySet().stream().map(x->Integer.valueOf(x)).sorted(((o1, o2) -> o1-o2)).collect(Collectors.toList());
                for (Integer pos: mutPositions) {
                    if (!commonSites.contains(pos))
                        continue;
                    HashMap<String, Integer> majorMinorRecords = sample1Reads.get(Integer.toString(pos));
                    List<Map.Entry<String, Integer>> sorted = majorMinorRecords.entrySet().stream().sorted(((o1, o2) -> o1.getValue()-o2.getValue())).collect(Collectors.toList());
                    ArrayList<Integer> snvs = this.sample1Major.getOrDefault(label, new ArrayList<>());
                    snvs.add(sorted.get(1).getValue());
                    this.sample1Major.put(label, snvs);
                    snvs = this.sample1Minor.getOrDefault(label, new ArrayList<>());
                    snvs.add(sorted.get(0).getValue());
                    this.sample1Minor.put(label, snvs);
                }
            }
            ArrayList<String> component2 = sample2Component.get(i);
            for (String peakRange: component2) {
                String geneName = this.sample1PeakCoveredGene.get(chrNum+":"+peakRange);
                HashSet<String> genes = this.mergedPeakCoveredGenes.getOrDefault(label, new HashSet<>());
                genes.add(geneName);
                this.mergedPeakCoveredGenes.put(label, genes);
                sample2Reads = this.sample2PeakSnpReadsCount.get(chrNum).get(peakRange);
                List<Integer> mutPositions = sample2Reads.keySet().stream().map(x->Integer.valueOf(x)).sorted(((o1, o2) -> o1-o2)).collect(Collectors.toList());
                for (Integer pos: mutPositions) {
                    if (!commonSites.contains(pos))
                        continue;
                    HashMap<String, Integer> majorMinorRecords = sample2Reads.get(Integer.toString(pos));
                    List<Map.Entry<String, Integer>> sorted = majorMinorRecords.entrySet().stream().sorted(((o1, o2) -> o1.getValue()-o2.getValue())).collect(Collectors.toList());
                    ArrayList<Integer> snvs = this.sample2Major.getOrDefault(label, new ArrayList<>());
                    snvs.add(sorted.get(1).getValue());
                    this.sample2Major.put(label, snvs);
                    snvs = this.sample2Minor.getOrDefault(label, new ArrayList<>());
                    snvs.add(sorted.get(0).getValue());
                    this.sample2Minor.put(label, snvs);
                }
            }
        }
    }

    /**
     * sample specific ASM test
     */
    private void asmDifferentiationAnalysis() {
        this.logger.debug(this.sample1Major.size() + " merged peaks can be used for sample specific ASM test");
        this.lock = new ReentrantLock();
        this.peakSampleSpecificPValue = new HashMap<>();
        this.peakMajorAlleleFrequency = new HashMap<>();
        this.peakMajorMinorReads = new HashMap<>();
        ExecutorService threadPoolExecutor = Executors.newFixedThreadPool(this.threadNumber);
        CountDownLatch countDown = new CountDownLatch(this.sample1Major.size());
        long tenPercent = Math.round(this.sample1Major.size() * 0.1);

        ArrayList<int[]> chrPeaks;
        for (String chr: this.chrMergedPeaks.keySet()) {
            chrPeaks = this.chrMergedPeaks.get(chr);
            for (int[] peakRange: chrPeaks) {
                Runnable task = () -> {
                    String label = chr+":"+peakRange[0]+":"+peakRange[1];
                    int size = this.sample1Major.get(label).size();
                    ArrayList<Integer> s1Major = this.sample1Major.get(label);
                    ArrayList<Integer> s1Minor = this.sample1Minor.get(label);
                    ArrayList<Integer> s2Major = this.sample2Major.get(label);
                    ArrayList<Integer> s2Minor = this.sample2Minor.get(label);
                    int[] major1 = new int[size], minor1 = new int[size], major2 = new int[size], minor2 = new int[size];
                    for (int i=0; i<size; i++) {
                        major1[i] = s1Major.get(i);
                        minor1[i] = s1Minor.get(i);
                        major2[i] = s2Major.get(i);
                        minor2[i] = s2Minor.get(i);
                    }

                    HierarchicalBayesianModel hb;
                    try {
                        this.peakSNVNum.put(label, s1Major.size());
                        int majorCount1 = s1Major.stream().reduce((x, y) -> x + y).get();
                        int minorCount1 = s1Minor.stream().reduce((x, y) -> x + y).get();
                        int majorCount2 = s2Major.stream().reduce((x, y) -> x + y).get();
                        int minorCount2 = s2Minor.stream().reduce((x, y) -> x + y).get();
                        this.peakMajorMinorReads.put(label, new int[]{majorCount1, minorCount1, majorCount2, minorCount2});
                        double df, aveDepth, scaleParam, epsilon = 0.00000000001;
                        df = Math.max(3, major1.length);
                        aveDepth = Arrays.stream(major1).average().getAsDouble();
                        scaleParam = (aveDepth < 10)? 50: 100;
                        if (minorCount1 == 0 && minorCount2 == 0)
                            hb = new HierarchicalBayesianModel(df, scaleParam, this.samplingTime,
                                    this.burnIn, major1, major1, major2, major2);
                        else
                            hb = new HierarchicalBayesianModel(df, scaleParam, this.samplingTime,
                                    this.burnIn, major1, minor1, major2, minor2);
                        boolean wellDone = false;
                        int maxTrail = 0;
                        double pVal = 0, peakOddRatio, peakMAF = 0;
                        while (!wellDone && maxTrail < 10) {  // avoid MCMC failed
                            pVal = hb.testSignificant();
                            peakOddRatio = Math.exp(hb.quantifyGeneLOR());
                            peakMAF = Math.min(1.0, peakOddRatio / (peakOddRatio + 1));
                            Double maf = new Double(peakMAF);
                            if (!maf.equals(Double.NaN) && !(Math.abs(peakMAF - 1.0) < epsilon) && !(Math.abs(peakMAF) < epsilon))
                                wellDone = true;
                            maxTrail++;
                        }

                        this.lock.lock();
                        ArrayList<String> samePValPeaks = this.peakSampleSpecificPValue.getOrDefault(pVal, new ArrayList<>());
                        samePValPeaks.add(label);
                        this.peakSampleSpecificPValue.put(pVal, samePValPeaks);
                        this.peakMajorAlleleFrequency.put(label, peakMAF);
                    } catch (Exception e) {
                        this.logger.error("error occurs on record " + label);
                        this.logger.error(e.getMessage());
                    } finally {
                        countDown.countDown();
                        if (countDown.getCount() % tenPercent == 0) {
                            double proportion = 100 - 10.0 * countDown.getCount() / tenPercent;
                            logger.debug(proportion + "% completed");
                        }
                        this.lock.unlock();
                        hb = null;
                        s1Major = null;
                        s1Minor = null;
                        s2Major = null;
                        s2Minor = null;
                        major1 = null;
                        minor1 = null;
                        major2 = null;
                        minor2 = null;
                    }
                };
                threadPoolExecutor.submit(task);
            }
        }

        try {
            countDown.await();
        } catch (InterruptedException ie) {
            this.logger.error("analysis interrupted");
            this.logger.error(ie.getMessage());
        } finally {
            this.sample1AsmPeakDetector = null;
            this.sample2AsmPeakDetector = null;
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
    private void bhRecalibrationOfEachPeak() {
        this.logger.debug("start recalibrating sample specific ASM peak p values");
        this.logger.debug("sorting test result in order");
        ArrayList<Map.Entry<Double, ArrayList<String>>> sortedPVals = new ArrayList<>(this.peakSampleSpecificPValue.entrySet());
        // sort p value from large to small
        Collections.sort(sortedPVals, new Comparator<Map.Entry<Double, ArrayList<String>>>() {
            @Override
            public int compare(Map.Entry<Double, ArrayList<String>> o1, Map.Entry<Double, ArrayList<String>> o2) {
                return o2.getKey().compareTo(o1.getKey());
            }
        });

        int totalPeak = this.sample1Major.size(), rankage = totalPeak;
        double prevQValue = 1.0, qValue;
        String pValString, qValString;
        this.peakSampleSpecificQValue = new ArrayList<>(totalPeak);

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

                pValString = this.df.format(pVal);
                qValString = this.df.format(qValue);
                this.peakSampleSpecificQValue.add(String.join("->", new String[]{peak, pValString, qValString}));
            }
        }
        this.logger.debug("recalibration complete.");
    }

    /**
     * test result output
     */
    private void output() {
        ArrayList<String[]> outputRecord = new ArrayList<>();
        BufferedWriter bfw = null;
        try {
            bfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(this.outputFile))));
            String line, label, geneId, geneName, chrNum, peakStart, peakEnd, sample1MajorReads, sample1MinorReads,
                    sample2MajorReads, sample2MinorReads, pValue, qValue;
            String[] info, rec, finalInfo;
            int[] twoSampleReads;
            int snvNum;
            bfw.write("#chr\tpeakStart\tpeakEnd\tgeneId\tgeneName\tp-value\tq-value\tsnvNum\tsample1Major/minorAlleleReads\tsample2Major/minorAlleleReads\n");
            for (String record: this.peakSampleSpecificQValue) {
                rec = record.split("->");
                label = rec[0];
                // info = [chrNum, peakStart, peakEnd]
                info = label.split(":");
                chrNum = info[0];
                peakStart = info[1];
                peakEnd = info[2];
                pValue = rec[1];
                qValue = rec[2];
                HashSet<String> peakGenes = this.mergedPeakCoveredGenes.get(label);
                geneId = peakGenes.stream().filter(x -> x!=null).collect(Collectors.joining());
                geneName = this.getGeneNames(chrNum, peakGenes);
                snvNum = this.peakSNVNum.get(label);
                twoSampleReads = this.peakMajorMinorReads.get(label);
                sample1MajorReads = String.valueOf(twoSampleReads[0]);
                sample1MinorReads = String.valueOf(twoSampleReads[1]);
                sample2MajorReads = String.valueOf(twoSampleReads[2]);
                sample2MinorReads = String.valueOf(twoSampleReads[3]);
                finalInfo = new String[]{chrNum, peakStart, peakEnd, geneId, geneName, pValue, qValue,
                        Integer.toString(snvNum), sample1MajorReads + "," + sample1MinorReads,
                        sample2MajorReads+","+sample2MinorReads};
                outputRecord.add(finalInfo);
            }

            Collections.sort(outputRecord, new Comparator<String[]>() {
                @Override
                public int compare(String[] o1, String[] o2) {
                    Double q1 = Double.parseDouble(o1[6]), q2 = Double.parseDouble(o2[6]);
                    if (!q1.equals(q2))
                        return q1.compareTo(q2);
                    Double p1 = Double.parseDouble(o1[5]), p2 = Double.parseDouble(o2[5]);
                    if (!p1.equals(p2))
                        return p1.compareTo(p2);
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
            this.logger.debug("result file " + this.outputFile);
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

    private void parseGTFFile(String gtfFile) {
        BufferedReader bfr = null;
        this.geneNames = new HashMap<>();
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(gtfFile)));
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

    private String getGeneNames(String chr, HashSet<String> geneIds) {
        ArrayList<String> names = new ArrayList<>();
        for (String id: geneIds) {
            String name = this.geneNames.get(chr).getOrDefault(id, "unknown");
            if (name.equals("unknown"))
                continue;
            names.add(name);
        }
        return names.stream().filter(x -> !x.equals("unknown")).collect(Collectors.joining());
    }

    private static Logger initLog(String logHome) {
        System.setProperty("log_home", logHome);
        return Logger.getLogger(SampleSpecificASM.class);
    }

    private static CommandLine setCommandLine(String[] args, Options options) throws ParseException {
        Option option = new Option("s1_vcf", "sample1_vcf_file", true, "sample1 SNP calling result VCF file, required");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("s1_bed", "sample1_bed_file", true, "sample1 peak calling result BED file, required");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("s1_wes", "sample1_wes_vcf_file", true, "sample1 WES data SNP calling VCF format file, optional");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("s2_vcf", "sample2_vcf_file", true, "sample2 SNP calling result VCF file, required");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("s2_bed", "sample2_bed_file", true, "sample2 peak calling result BED file, required");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("s2_wes", "sample2_wes_vcf_file", true, "sample2 WES data SNP calling VCF format file, optional");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("g", "gtf", true, "GTF annotation file, required");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("db", "dbsnp", true, "big scale SNV annotation data set, like dbsnp, 1000Genome etc. Optional, the file format see https://github.com/Jakob666/allele-specificM6A");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("o", "output", true, "ASM peak test output file. Optional, default ./sampleSpecificASM.txt");
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

        option = new Option("b", "burn", true, "burn-in times, more than 100 and less than sampling times. Optional, default 10000");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("t", "thread", true, "thread number for running test. Optional, Default 2");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("h", "help", false, "help message of SampleSpecificASM");
        option.setRequired(false);
        options.addOption(option);

        CommandLineParser parser = new DefaultParser();

        return parser.parse(options, args);
    }
}
