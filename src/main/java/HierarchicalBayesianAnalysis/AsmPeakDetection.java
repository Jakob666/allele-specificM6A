package HierarchicalBayesianAnalysis;


import AseM6aPeakDetector.HeterozygoteReadsCount;
import heterozygoteSiteAnalysis.DbsnpAnnotation;
import org.apache.commons.cli.*;
import org.apache.log4j.Logger;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

/**
 * Test ASM peaks with VCF format file and BED format file
 */
public class AsmPeakDetection {
    private String peakBedFile, vcfFile, wesFile, asmPeakFile, peakCoveredSnpFile, peakCoveredWesSnpFile;
    private int ipSNPReadInfimum, wesSNPReadInfimum, samplingTime, burnIn, totalPeakCount = 0;
    private double degreeOfFreedom;
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
    private ArrayList<String> asmQValue = new ArrayList<>();
    private DecimalFormat df = new DecimalFormat("0.0000");

    /**
     * Constructor
     * @param peakBedFile BED format file via MeRIP-seq IP data
     * @param vcfFile VCF format file via MeRIP-seq INPUT data
     * @param wesFile VCF format file via WES data, optional
     * @param dbsnpFile dbsnp file for SNP filtering
     * @param peakCoveredSnpFile output file which record MeRIP-seq INPUT data SNV sites covered by m6A signal
     * @param peakCoveredWesSnpFile output file which record WES data SNV sites covered by m6A signal
     * @param asmPeakFile test result output file
     * @param degreeOfFreedom the degree of freedom of inverse-Chi-square distribution, default 10
     * @param ipSNPReadInfimum reads coverage threshold when filter INPUT sample SNV sites, default 10
     * @param wesSNPReadInfimum reads coverage threshold when filter WES SNV sites, default 30
     * @param samplingTime sampling time, default 5000
     * @param burnIn burn in time, default 200
     */
    public AsmPeakDetection(String peakBedFile, String vcfFile, String wesFile, String dbsnpFile, String peakCoveredSnpFile,
                            String peakCoveredWesSnpFile, String asmPeakFile, double degreeOfFreedom,
                            int ipSNPReadInfimum, int wesSNPReadInfimum, int samplingTime, int burnIn) {
        this.peakBedFile = peakBedFile;
        this.vcfFile = vcfFile;
        this.wesFile = wesFile;
        this.peakCoveredSnpFile = peakCoveredSnpFile;
        this.peakCoveredWesSnpFile = peakCoveredWesSnpFile;
        this.log = initLog(new File(asmPeakFile).getParent());
        if (this.wesFile != null)
            assert this.peakCoveredWesSnpFile != null;
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
    }

    public static void main(String[] args) throws ParseException {
        Options options = new Options();
        CommandLine commandLine = setCommandLine(args, options);

        String bedFile = null, aseVcfFile = null, wesVcfFile = null, dbsnpFile = null, outputFile, outputDir,
               peakCoveredSnpFile, peakCoveredSnpBackgroundFile;
        int ipSNPCoverageInfimum = 10, wesSNPCoverageInfimum = 30, samplingTime = 5000, burn_in = 200;
        double degreeOfFreedom = 10;

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

        if (commandLine.hasOption("db")) {
            File dbsnp = new File(commandLine.getOptionValue("db"));
            if (!dbsnp.exists() || !dbsnp.isFile()) {
                logger.error("invalid file path: " + dbsnp.getAbsolutePath());
                System.exit(2);
            }
            dbsnpFile = dbsnp.getAbsolutePath();
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
        if (commandLine.hasOption("df"))
            degreeOfFreedom = Double.parseDouble(commandLine.getOptionValue("df"));
        if (degreeOfFreedom <= 1) {
            System.out.println("invalid inverse-Chi-square distribution parameter for tau sampling. Must larger than 1.0");
            System.exit(2);
        }
        if (commandLine.hasOption("rc"))
            ipSNPCoverageInfimum = Integer.parseInt(commandLine.getOptionValue("rc"));
        if (commandLine.hasOption("wc"))
            wesSNPCoverageInfimum = Integer.parseInt(commandLine.getOptionValue("wc"));

        AsmPeakDetection apd = new AsmPeakDetection(bedFile, aseVcfFile, wesVcfFile, dbsnpFile, peakCoveredSnpFile,
                                                    peakCoveredSnpBackgroundFile, outputFile, degreeOfFreedom, ipSNPCoverageInfimum,
                                                    wesSNPCoverageInfimum, samplingTime, burn_in);
        apd.getTestResult();
    }

    public void getTestResult() {
        this.getPeakCoveredGene();
        this.getPeakCoveredSnpResult();
        if (this.wesFile != null)
            this.getPeakCoveredSnpBackground();
        this.asmPeakTest();
        this.bhRecalibrationOfEachPeak();
        this.outputResult();
    }

    /**
     * locate m6A signal in BED format file to corresponding gene, peakCoveredGene = {chr:peakStart:peakEnd -> geneId, ...}
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
     * get m6A signal covered MeRIP-seq INPUT sample SNV sites
     * #chr	SNP strand	position	peakStart	peakEnd	majorAlleleStrand	majorCount	minorCount
     */
    private void getPeakCoveredSnpResult() {
        PeakCoveredSNPRecord pcsr = new PeakCoveredSNPRecord(this.vcfFile, this.peakBedFile,
                                                             this.peakCoveredSnpFile, this.ipSNPReadInfimum);
        pcsr.getPeakCoveredSNP();
    }

    /**
     * get m6A signal covered WES SNV sites
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
     * get major allele and minor allele nucleotide and reads counts of each MeRIP-seq INPUT sample SNV sites covered by m6A peak
     * [chr1: [peak1: position1: [major: count, minor: count], position2:[major: count, minor:count]], chr2:....]
     */
    private void getPeakSNPReadsCount() {
        HeterozygoteReadsCount hrc = new HeterozygoteReadsCount(this.peakCoveredSnpFile, this.log);
        this.peakSnpReadsCount = hrc.getMajorMinorHaplotype();
        this.peakMajorAlleleNucleotide = hrc.getPeakMajorAlleleNucleotides();
    }

    /**
     * get major allele and minor allele nucleotide and reads counts of each WES SNV sites covered by m6A peak
     * @return [chr1: [peak1: position1: [major: count, minor: count], position2:[major: count, minor:count]], chr2:....]
     */
    private HashMap<String, HashMap<String, HashMap<String, HashMap<String, Integer>>>> getPeakSNPBackground() {
        HeterozygoteReadsCount hrc = new HeterozygoteReadsCount(this.peakCoveredWesSnpFile, this.log);
        return hrc.getMajorMinorHaplotype();
    }

    /**
     * test the ASM significant p value of a m6A peak
     */
    private void asmPeakTest() {
        // [chr1: [peak1: position1: [major: count, minor: count], position2:[major: count, minor:count]], chr2:....]
        this.getPeakSNPReadsCount();
        if (this.wesFile != null)
            this.peakSnpBackground = this.getPeakSNPBackground();
        else
            this.peakSnpBackground = null;
        this.hierarchicalModelTest();
    }

    /**
     * calculate m6A signal p value via hierarchical model
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
            // m6A signal on a particular chromosome
            rnaSeqPeakSnvAlleleReads = this.peakSnpReadsCount.get(chrNum);

            if (this.peakSnpBackground != null)
                wesPeakSnvAlleleReads = this.peakSnpBackground.getOrDefault(chrNum, null);
            else
                wesPeakSnvAlleleReads = null;

            for (String peakRange: rnaSeqPeakSnvAlleleReads.keySet()) {
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
                maf = this.calcMajorAlleleFrequency(majorCount, minorCount);
                // MAF threshold for reducing false positive ASM peaks. Peaks with MAF < 0.55 are seen as normal peak
                // caused by alignment error
                if (maf - 0.55 < 0.00001)
                    continue;
                // MAF threshold for reducing false positive ASE genes. Genes with MAF > 0.95 are seen as homo-zygote gene
                // cause by alignment error
                if (maf - 0.95 > 0.00001)
                    continue;
                this.peakMajorAlleleFrequency.put(label, maf);
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

        double lorStd = this.calcLorStd();
        for (String str: this.statisticForTest.keySet()) {
            this.totalPeakCount++;
            ArrayList<int[]> statistic = this.statisticForTest.get(str);
            majorCount = statistic.get(0);
            minorCount = statistic.get(1);
            majorBackground = statistic.get(2);
            minorBackground = statistic.get(3);
            // get p value via hierarchical model
            HierarchicalBayesianModel hb = new HierarchicalBayesianModel(lorStd, this.degreeOfFreedom, this.samplingTime,
                    this.burnIn, majorCount, minorCount, majorBackground, minorBackground);
            pVal = hb.testSignificant();
            ArrayList<String> samePValPeaks = this.asmPValue.getOrDefault(pVal, new ArrayList<>());
            samePValPeaks.add(str);
            this.asmPValue.put(pVal, samePValPeaks);
            hb = null;
        }
    }

    /**
     * calculate the standard deviation of LOR of all SNV sites covered by m6A signals on genome
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
     * calculate major allele frequency of SNV site
     * @param majorCounts major allele reads count
     * @param minorCounts minor allele reads count
     * @return major allele frequency
     */
    private double calcMajorAlleleFrequency(int[] majorCounts, int[] minorCounts) {
        int major = 0, minor = 0;
        for (int i=0; i<majorCounts.length; i++) {
            major += majorCounts[i];
            minor += minorCounts[i];
        }

        return (double) major / (double) (major + minor);
    }

    /**
     * recalibrate the p value with BH method, get significant q value
     */
    private void bhRecalibrationOfEachPeak() {
        ArrayList<Map.Entry<Double, ArrayList<String>>> sortedPVals = new ArrayList<>(this.asmPValue.entrySet());
        // sort p value from small to large
        Collections.sort(sortedPVals, new Comparator<Map.Entry<Double, ArrayList<String>>>() {
            @Override
            public int compare(Map.Entry<Double, ArrayList<String>> o1, Map.Entry<Double, ArrayList<String>> o2) {
                return o2.getKey().compareTo(o1.getKey());
            }
        });

        int totalPeak = this.totalPeakCount, rankage = totalPeak;
        double qValue;
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
                rankage--;

                pValString = Double.toString(pVal);
                qValString = this.df.format(qValue);
                this.asmQValue.add(String.join("->", new String[]{peak, pValString, qValString}));
            }
        }
    }

    /**
     * write the test result into file
     */
    private void outputResult() {
        ArrayList<String[]> outputRecord = new ArrayList<>();
        BufferedWriter bfw = null;
        try {
            bfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(this.asmPeakFile))));
            String line, label, geneId, chrNum, peakStart, peakEnd, majorAlleleReads, minorAlleleReads,
                   majorAlleleBackground, minorAlleleBackground, pValue, qValue;
            LinkedList<String> majorAlleleNucleotides;
            String[] info, rec, finalInfo;
            String majorAlleleFrequency;
            int snvNum;
            bfw.write("#chr\tpeakStart\tpeakEnd\tgeneId\tp-value\tq-value\tsnvNum\tmajor/minorAlleleReads\tmajor/minorBackground\tMajorAlleleFrequency\tmajorAlleleNC\n");
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
                majorAlleleFrequency = this.df.format(this.peakMajorAlleleFrequency.get(label));
                majorAlleleNucleotides = this.peakMajorAlleleNucleotide.get(label);
                String majorAlleleRecords = this.getString(majorAlleleNucleotides);
                snvNum = this.peakSNVNum.get(label);
                majorAlleleReads = String.valueOf(this.peakMajorMinorAlleleCount.get(label).get("major"));
                minorAlleleReads = String.valueOf(this.peakMajorMinorAlleleCount.get(label).get("minor"));
                majorAlleleBackground = String.valueOf(this.peakMajorMinorBackground.get(label).get("major"));
                minorAlleleBackground = String.valueOf(this.peakMajorMinorBackground.get(label).get("minor"));
                finalInfo = new String[]{chrNum, peakStart, peakEnd, geneId, pValue, qValue,
                                         Integer.toString(snvNum), majorAlleleReads + "," + minorAlleleReads,
                                         majorAlleleBackground+","+minorAlleleBackground, majorAlleleFrequency,
                                         majorAlleleRecords};
                outputRecord.add(finalInfo);
            }

            Collections.sort(outputRecord, new Comparator<String[]>() {
                @Override
                public int compare(String[] o1, String[] o2) {
                    Double q1 = Double.parseDouble(o1[5]), q2 = Double.parseDouble(o2[5]);
                    if (!q1.equals(q2))
                        return q1.compareTo(q2);
                    // sort records with same q-value after BH recalibration by SNV number
                    Integer snvCount1 = Integer.parseInt(o1[6]), snvCount2 = Integer.parseInt(o2[6]);
//                    if (snvCount1 - snvCount2 != 0)
                    return snvCount2.compareTo(snvCount1);
                    // sort records with same q-value and SNV number with MAF
//                    Double maf1 = Double.parseDouble(o1[9]), maf2 = Double.parseDouble(o2[9]);
//                    return maf2.compareTo(maf1);
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

        option = new Option("wes", "wes_vcf_file", true, "WES SNP calling VCF format file");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("bed", "peak_bed_file", true, "Peak calling output result in BED format");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("db", "dbsnp", true, "dbsnp file for SNP filtering");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("o", "output", true, "ASM peak test output file, default ./asmPeak.txt");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("peak_cover_snp", "peak_cover_snp_record_file", true, "Peak covered INPUT sample SNP record");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("peak_cover_snp_bkg", "peak_cover_snp_background_record_file", true, "Peak covered WES sample SNP record");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("df", "degree_of_freedom", true, "degree of freedom of inverse-Chi-square distribution, default 10");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("rc", "reads_coverage", true, "IP sample SNP site coverage infimum, default 10");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("wc", "wes_coverage", true, "WES sample SNP site coverage infimum, default 30");
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
