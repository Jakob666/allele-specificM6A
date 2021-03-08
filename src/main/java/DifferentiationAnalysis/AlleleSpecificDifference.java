package DifferentiationAnalysis;


import HierarchicalBayesianAnalysis.AseGeneDetection;
import HierarchicalBayesianAnalysis.AsmPeakDetection;
import HierarchicalBayesianAnalysis.OddRatioCalc;
import HierarchicalBayesianAnalysis.RunTest;
import org.apache.commons.math3.util.MathArrays;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.locks.ReentrantLock;
import java.util.function.Function;

public class AlleleSpecificDifference {
    private String s1VcfFile, s2VcfFile, s1WesFile, s2WesFile, dbsnpFile, bedFile, gtfFile, outputFile;
    private int readsCoverageThreshold, wesCoverageThreshold, samplingTime, burnIn, threadNum;
    private Logger logger;
    private ReentrantLock lock;
    private ExecutorService executorService;
    private HashMap<String, String> bf;
    private HashMap<String, double[]> quantifiedResult;
    private HashMap<String, HashMap<String, int[]>> statisticForTest;

    private AlleleSpecificDifference(String s1VcfFile, String s2VcfFile, String s1WesFile, String s2WesFile,
                                     String dbsnpFile, String bedFile, String gtfFile, String outputFile,
                                     int readsCoverageThreshold, int wesCoverageThreshold, int samplingTime, int burnIn, int threadNum) {
        this.s1VcfFile = s1VcfFile;
        this.s2VcfFile = s2VcfFile;
        this.s1WesFile = s1WesFile;
        this.s2WesFile = s2WesFile;
        this.dbsnpFile = dbsnpFile;
        this.bedFile = bedFile;
        this.gtfFile = gtfFile;
        this.outputFile = outputFile;
        this.readsCoverageThreshold = readsCoverageThreshold;
        this.wesCoverageThreshold = wesCoverageThreshold;
        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
        this.threadNum = threadNum;
    }

    public void prepare() {
        if (this.bedFile == null) {
            String s1SnvFile = new File(new File(this.outputFile).getParent(), "s1_snv_location.txt").getAbsolutePath();
            AseGeneDetection s1ase = new AseGeneDetection(this.gtfFile, this.s1VcfFile, this.s1WesFile, this.dbsnpFile, s1SnvFile, this.outputFile,
                                                          this.readsCoverageThreshold, this.wesCoverageThreshold, this.samplingTime, this.burnIn, this.threadNum, this.logger);
            s1ase.dataPreparation();

            String s2SnvFile = new File(new File(this.outputFile).getParent(), "s2_snv_location.txt").getAbsolutePath();
            AseGeneDetection s2ase = new AseGeneDetection(this.gtfFile, this.s2VcfFile, this.s2WesFile, this.dbsnpFile, s2SnvFile, this.outputFile,
                    this.readsCoverageThreshold, this.wesCoverageThreshold, this.samplingTime, this.burnIn, this.threadNum, this.logger);
            s2ase.dataPreparation();

            // geneAlleleReads = {"geneId->geneName": {pos1: [majorAllele:count, minorAllele: count1]}, ...}
            HashMap<String, HashMap<Integer, String[]>> s1GeneAlleleReads = s1ase.getGeneAlleleReads();
            HashMap<String, HashMap<Integer, String[]>> s2GeneAlleleReads = s2ase.getGeneAlleleReads();
            HashSet<String> sample1Genes = new HashSet<> (s1GeneAlleleReads.keySet()), sample2Genes = new HashSet<>(s2GeneAlleleReads.keySet());
            // sample1Genes represents common genes
            sample1Genes.retainAll(sample2Genes);
            if (sample1Genes.isEmpty()) {
                this.logger.debug("detect no common genes for allele-specific expression differentiation analysis between two samples");
                System.exit(0);
            }

            this.statisticForTest = new HashMap<>();
            HashMap<String, String> s1GeneMajorNucleotide = s1ase.getGeneMajorNucleotide();
            HashMap<String, String> s2GeneMajorNucleotide = s2ase.getGeneMajorNucleotide();

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
            HashMap<String, int[]> temp;
            for (String gene: sample1Genes) {
                info = gene.split("->");
                geneId = info[0];
                String sample1GeneRecord = s1GeneMajorNucleotide.get(geneId);
                HashSet<Integer> sample1GeneSNVPositions = getGeneSNVPositions.apply(sample1GeneRecord);
                String sample2GeneRecord = s2GeneMajorNucleotide.get(geneId);
                HashSet<Integer> sample2GeneSNVPositions = getGeneSNVPositions.apply(sample2GeneRecord);

                sample1Reads = this.getGeneSNVReadsCount(s1GeneAlleleReads.get(gene), sample1GeneSNVPositions);
                sample1MajorCounts = sample1Reads[0];
                sample1MinorCounts = sample1Reads[1];
                sample2Reads = this.getGeneSNVReadsCount(s2GeneAlleleReads.get(gene), sample2GeneSNVPositions);
                sample2MajorCounts = sample2Reads[0];
                sample2MinorCounts = sample2Reads[1];
                // common SNV sites major allele nucleotide
//                this.formMajorAlleleRecord(gene, sample1GeneSNVPositions);
                temp = new HashMap<>();
                temp.put("s1Major", sample1MajorCounts);
                temp.put("s1Minor", sample1MinorCounts);
                temp.put("s2Major", sample2MajorCounts);
                temp.put("s2Minor", sample2MinorCounts);

                this.statisticForTest.put(gene, temp);
            }
            if (this.statisticForTest.isEmpty()) {
                this.logger.error("contains no common genes with SNV sites between two samples for differentiation test, please check the input data");
                System.exit(2);
            } else
                this.logger.debug("data preparation completed");
        } else {
            String s1PeakWithSnvFile = new File(new File(outputFile).getParent(), "s1_peak_with_snv.txt").getAbsolutePath();
            String s1PeakWithSnvBkgFile = null;
            if (this.s1WesFile!=null)
                s1PeakWithSnvBkgFile = new File(new File(outputFile).getParent(), "s1_peak_with_snv_bkg.txt").getAbsolutePath();
            AsmPeakDetection s1asm = new AsmPeakDetection(this.gtfFile, this.bedFile, this.s1VcfFile, this.s1WesFile, this.dbsnpFile,
                    s1PeakWithSnvFile, s1PeakWithSnvBkgFile, this.outputFile, this.readsCoverageThreshold, this.wesCoverageThreshold,
                    this.samplingTime, this.burnIn, this.threadNum, this.logger);
            s1asm.getPeakSNPReadsCount();
            s1asm.dataPreparation();

            String s2PeakWithSnvFile = new File(new File(outputFile).getParent(), "s2_peak_with_snv.txt").getAbsolutePath();
            String s2PeakWithSnvBkgFile = null;
            if (this.s1WesFile!=null)
                s2PeakWithSnvBkgFile = new File(new File(outputFile).getParent(), "s2_peak_with_snv_bkg.txt").getAbsolutePath();
            AsmPeakDetection s2asm = new AsmPeakDetection(this.gtfFile, this.bedFile, this.s1VcfFile, this.s1WesFile, this.dbsnpFile,
                    s2PeakWithSnvFile, s2PeakWithSnvBkgFile, this.outputFile, this.readsCoverageThreshold, this.wesCoverageThreshold,
                    this.samplingTime, this.burnIn, this.threadNum, this.logger);
            s2asm.getPeakSNPReadsCount();
            s2asm.dataPreparation();
            // whether m6A peaks cover common genes between two samples
            // chr:peakStart:peakEnd -> geneId
            HashMap<String, String> s1PeakCoveredGene = s1asm.getPeakCoveredGenes();
            HashMap<String, String> s2PeakCoveredGene = s2asm.getPeakCoveredGenes();
            // [chr1->{peak1->{pos1->{major: count, minor: count}, pos2:{major: count, minor:count},...},...}, chr2:....]
            // peak = peakStart:peakEnd
            HashMap<String, HashMap<String, HashMap<String, HashMap<String, Integer>>>> s1PeakSnpReadsCount = s1asm.getPeakSnpReadsCount();
            HashMap<String, HashMap<String, HashMap<String, HashMap<String, Integer>>>> s2PeakSnpReadsCount = s2asm.getPeakSnpReadsCount();

            HashSet<String> sample1Genes = new HashSet<>(), sample2Genes = new HashSet<>();
            s1PeakCoveredGene.entrySet().forEach(entry -> sample1Genes.add(entry.getValue()));
            s2PeakCoveredGene.entrySet().forEach(entry -> sample2Genes.add(entry.getValue()));
            sample1Genes.retainAll(sample2Genes);
            if (sample1Genes.isEmpty()) {
                this.logger.debug("detect no common peak covered genes for allele-specific modification differentiation analysis between two samples");
                System.exit(0);
            }

            this.statisticForTest = new HashMap<>();
            int idx;
            int[] sample1MajorCounts, sample1MinorCounts, sample2MajorCounts, sample2MinorCounts;
            Integer[] temp, readsCount;
            HashMap<String, int[]> record;
            for (String chr: s1PeakSnpReadsCount.keySet()) {
                if (!s2PeakSnpReadsCount.containsKey(chr))
                    continue;
                HashMap<String, HashMap<String, HashMap<String, Integer>>> s1PeaksRecord = s1PeakSnpReadsCount.get(chr);
                HashMap<String, HashMap<String, HashMap<String, Integer>>> s2PeaksRecord = s2PeakSnpReadsCount.get(chr);
                for (String peak: s1PeaksRecord.keySet()) {
                    String peakGene = s1PeakCoveredGene.get(peak);
                    if (!s2PeaksRecord.containsKey(peak))
                        continue;
                    HashMap<String, HashMap<String, Integer>> s1SnpsRecord = s1PeaksRecord.get(peak);
                    sample1MajorCounts = new int[s1SnpsRecord.size()];
                    sample1MinorCounts = new int[s1SnpsRecord.size()];
                    idx = 0;
                    for (String pos: s1SnpsRecord.keySet()) {
                        temp = new Integer[2];
                        readsCount = s1SnpsRecord.get(pos).values().toArray(temp);
                        if (readsCount[0] > readsCount[1]) {
                            sample1MajorCounts[idx] = readsCount[0];
                            sample1MinorCounts[idx] = readsCount[1];
                        } else {
                            sample1MajorCounts[idx] = readsCount[1];
                            sample1MinorCounts[idx] = readsCount[0];
                        }
                        idx++;
                    }

                    HashMap<String, HashMap<String, Integer>> s2SnpsRecord = s2PeaksRecord.get(peak);
                    sample2MajorCounts = new int[s2SnpsRecord.size()];
                    sample2MinorCounts = new int[s2SnpsRecord.size()];
                    idx = 0;
                    for (String pos: s2SnpsRecord.keySet()) {
                        temp = new Integer[2];
                        readsCount = s2SnpsRecord.get(pos).values().toArray(temp);
                        if (readsCount[0] > readsCount[1]) {
                            sample2MajorCounts[idx] = readsCount[0];
                            sample2MinorCounts[idx] = readsCount[1];
                        } else {
                            sample2MajorCounts[idx] = readsCount[1];
                            sample2MinorCounts[idx] = readsCount[0];
                        }
                        idx++;
                    }

                    String label = peakGene + ":" + peak;
                    record = new HashMap<>();
                    record.put("s1Major", sample1MajorCounts);
                    record.put("s1Minor", sample1MinorCounts);
                    record.put("s2Major", sample2MajorCounts);
                    record.put("s2Minor", sample2MinorCounts);
                    this.statisticForTest.put(label, record);
                }
            }
            if (this.statisticForTest.isEmpty()) {
                this.logger.error("contains no common genes with SNV sites between two samples for differentiation test, please check the input data");
                System.exit(2);
            } else
                this.logger.debug("data preparation completed");
        }
    }

    public void testDifference() {
        // test ASE gene with Hierarchical model
        this.logger.debug("hierarchical Bayesian model test start");
        // init thread pool and lock
        ExecutorService threadPoolExecutor = Executors.newFixedThreadPool(this.threadNum);
        this.lock = new ReentrantLock();
        CountDownLatch countDown = new CountDownLatch(this.statisticForTest.size());
        long tenPercent = Math.round(this.statisticForTest.size() * 0.1);

        this.bf = new HashMap<>();
        this.quantifiedResult = new HashMap<>();
        RunTest task = (String name) -> {
            return new Runnable() {
                @Override
                public void run() {
                    HashMap<String, int[]> reads;
                    int[] s1Major, s1Minor, s2Major, s2Minor, majorCount, minorCount;
                    double df, scaleParam, aveDepth, epsilon = 0.00000000001, methDiffThreshold=0.1;
                    ModelSelection sameModel, diffModel;
                    try {
                        reads = statisticForTest.get(name);
                        s1Major = reads.get("s1Major");
                        s1Minor = reads.get("s1Minor");
                        s2Major = reads.get("s2Major");
                        s2Minor = reads.get("s2Minor");
                        if (s1Major.length > s2Major.length) {
                            majorCount = s2Major;
                            minorCount = s2Minor;
                        } else {
                            majorCount = s1Major;
                            minorCount = s1Minor;
                        }
                        df = Math.max(2, majorCount.length);
                        aveDepth = Math.min(Arrays.stream(s1Major).average().getAsDouble(), Arrays.stream(s2Major).average().getAsDouble());
                        scaleParam = (aveDepth - 15 < epsilon)? 50: 100;
                        HashMap<String, double[]> s1LORAndVar = getObserveData(s1Major, s1Minor, null, null);
                        double[] s1ObserveLOR = s1LORAndVar.get("LOR");
                        double[] s1Variances = s1LORAndVar.get("VAR");
                        HashMap<String, double[]> s2LORAndVar = getObserveData(s2Major, s2Minor, null, null);
                        double[] s2ObserveLOR = s2LORAndVar.get("LOR");
                        double[] s2Variances = s2LORAndVar.get("VAR");

                        boolean welldone = false;
                        int time = 0, maxTrial = 10;
                        String choice = null, evidence;
                        double sameModelLogMarginProb = 0, diffModelLogMarginProb = 0, bayesianFactor = 0, methDiff;
                        double[] sameModelRes = null, diffModelRes = null;
                        while (!welldone && time < maxTrial) {
                            time += 1;
                            sameModel = new ModelSelection(s1ObserveLOR, s1Variances, s2ObserveLOR, s2Variances, true, 50000, 10000, df, scaleParam);
                            sameModel.initModel();
                            sameModel.sampling();
                            sameModelLogMarginProb = sameModel.calcModelMarginProb();
                            diffModel = new ModelSelection(s1ObserveLOR, s1Variances, s2ObserveLOR, s2Variances, false, 50000, 10000, df, scaleParam);
                            diffModel.initModel();
                            diffModel.sampling();
                            diffModelLogMarginProb = diffModel.calcModelMarginProb();
                            bayesianFactor = Math.exp(diffModelLogMarginProb-sameModelLogMarginProb);
                            double s1OR, s2OR, s1Maf, s2Maf;
                            sameModelRes = modelChoice(sameModel);
                            diffModelRes = modelChoice(diffModel);
                            methDiff = Math.abs(diffModelRes[0] - diffModelRes[1]);
                            if (bayesianFactor - 3 > epsilon) {
                                choice = "diffModel";
                                if (methDiff - methDiffThreshold > epsilon)
                                    welldone = true;
                            } else {
                                choice = "sameModel";
                                if (methDiff - methDiffThreshold < epsilon)
                                    welldone = true;
                            }

                        }
                        if (bayesianFactor - 1 < epsilon)
                            evidence = "No evidence support differential allele-specific";
                        else if (bayesianFactor -3 < epsilon)
                            evidence = "Anecdotal evidence support differential allele-specific";
                        else if (bayesianFactor - 10 < epsilon)
                            evidence = "Morderate evidence support differential allele-specific";
                        else
                            evidence = "Strong evidence support differential allele-specific";

                        lock.lock();
                        bf.put(name, evidence);
                        double[] res;
                        if (choice.equals("diffModel"))
                            quantifiedResult.put(name, diffModelRes);
                        else
                            quantifiedResult.put(name, sameModelRes);
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
                        sameModel = null;
                        diffModel = null;
                    }

                }
            };
        };

        this.logger.debug(this.statisticForTest.size() + " item to be tested");
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

    public void output() {
        BufferedWriter bfw = null;
        try {
            bfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(this.outputFile))));
            if (this.bedFile == null)
                bfw.write("#geneID\tgeneName\ts1LOR\ts2LOR\tevidence");
            else
                bfw.write("#geneID\tgeneName\tpeakStart\tpeakEnd\ts1LOR\ts2LOR\tevidence");
            bfw.newLine();
            String line;
            for (String label: this.bf.keySet()) {
                String bf = this.bf.get(label);
                double s1LOR = this.quantifiedResult.get(label)[0];
                double s2LOR = this.quantifiedResult.get(label)[1];
                line = String.join("\t", label.split(":")) + "\t" + String.join("\t", new String[]{Double.toString(s1LOR), Double.toString(s2LOR), bf});
                bfw.write(line);
                bfw.newLine();
            }
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

    private HashMap<String, double[]> getObserveData(int[] majorAlleleReads, int[] minorAlleleReads,
                                                     int[] majorAlleleBackground, int[] minorAlleleBackground) {
        OddRatioCalc orc = new OddRatioCalc(majorAlleleReads, minorAlleleReads, majorAlleleBackground, minorAlleleBackground);
        return orc.getLogOddRatio();
    }
    private double[] modelChoice(ModelSelection model) {
        double s1OR, s2OR, s1Maf, s2Maf;
        model.quantify();
        s1OR = Math.exp(model.getQuantifiedS1LOR());
        s2OR = Math.exp(model.getQuantifiedS2LOR());
        s1Maf = Math.min(1.0, s1OR / (s1OR + 1));
        s2Maf = Math.min(1.0, s2OR / (s2OR + 1));
        return new double[] {s1Maf, s2Maf};
    }

}
