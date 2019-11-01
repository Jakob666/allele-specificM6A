package HierarchicalBayesianAnalysis;

import org.apache.commons.cli.*;
import org.apache.log4j.Logger;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

public class AlleleSpecificM6aSignal {
    private File aseTestFile, asmTestFile, outputFile;
    private HashMap<String, String[]> aseSignificantGene = new HashMap<>(), aseNonsignificantGene = new HashMap<>(),

                                      significantResult = new HashMap<>(), nonSignificantResult = new HashMap<>();
    private HashMap<String, ArrayList<String[]>> asmSignificantPeak = new HashMap<>(), asmNonsignificantPeak = new HashMap<>();
    private Logger log;

    public static void main(String[] args) throws ParseException{
        Options options = new Options();
        CommandLine commandLine = setCommandLine(args, options);
        String aseGeneFile = null, asmPeakFile = null, outputFile = null;
        Logger logger;

        if (commandLine.hasOption("ase"))
            aseGeneFile = commandLine.getOptionValue("ase");
        if (commandLine.hasOption("asm"))
            asmPeakFile = commandLine.getOptionValue("asm");
        if (aseGeneFile==null) {
            System.out.println("ASE test result should be specified with -ase(--ase_gene_result)");
            System.exit(2);
        }
        if (asmPeakFile==null) {
            System.out.println("ASM test result should be specified with -asm(--asm_peak_result)");
            System.exit(2);
        }
        if (commandLine.hasOption("o"))
            outputFile = commandLine.getOptionValue("o");
        else
            outputFile = new File(new File(asmPeakFile).getParent(), "allele_specific_m6A_peak.txt").getAbsolutePath();

        logger = initLog(new File(asmPeakFile).getParent());
        AlleleSpecificM6aSignal asms = new AlleleSpecificM6aSignal(aseGeneFile, asmPeakFile, outputFile, logger);
        asms.detect();
    }

    /**
     * Constructor
     * @param aseTestFile ASE gene test result file
     * @param asmTestFile ASM peak result file
     * @param outputFile output file record allele-specific m6A peak
     * @param log Log4j Logger instance
     */
    public AlleleSpecificM6aSignal(String aseTestFile, String asmTestFile, String outputFile, Logger log) {
        this.aseTestFile = new File(aseTestFile);
        this.asmTestFile = new File(asmTestFile);
        this.outputFile = new File(outputFile);
        this.log = log;
    }

    /**
     * test if m6A signal is allele specific, and record it into file
     */
    public void detect() {
        this.parseAseTestFile();
        this.parseAsmTestFile();
        this.mergeAseAsmTestResult();
        this.outputResult();
    }

    /**
     * parse ASE Gene test result
     */
    private void parseAseTestFile() {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(this.aseTestFile)));
            String line = "", geneId, geneName, majorAlleleRecord;
            String[] info, record;
            double qVal;
            while (line != null) {
                line = bfr.readLine();
                if (line != null && !line.startsWith("#")) {
                    info = line.split("\t");
                    geneId = info[0];
                    geneName = info[1];
                    qVal = Double.parseDouble(info[3]);
                    majorAlleleRecord = info[8];
                    record = new String[] {geneName, info[3], majorAlleleRecord};
                    if (qVal - 0.05 < 0.000001)
                        this.aseSignificantGene.put(geneId, record);
                    else
                        this.aseNonsignificantGene.put(geneId, record);
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
     * parse ASM Peak test result
     */
    private void parseAsmTestFile() {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(this.asmTestFile)));
            String line = "", chrNum, peakStart, peakEnd, geneId, majorAlleleStrand;
            String[] info, record;
            ArrayList<String[]> recordList;
            double qVal;
            while (line != null) {
                line = bfr.readLine();
                if (line != null && !line.startsWith("#")) {
                    info = line.split("\t");
                    chrNum = info[0];
                    peakStart = info[1];
                    peakEnd = info[2];
                    geneId = info[3];
                    qVal = Double.parseDouble(info[5]);
                    majorAlleleStrand = info[10];

                    record = new String[] {chrNum, peakStart, peakEnd, info[4], majorAlleleStrand};

                    if (qVal - 0.05 < 0.00001) {
                        recordList = this.asmSignificantPeak.getOrDefault(geneId, new ArrayList<>());
                        recordList.add(record);
                        this.asmSignificantPeak.put(geneId, recordList);
                    } else {
                        recordList = this.asmNonsignificantPeak.getOrDefault(geneId, new ArrayList<>());
                        recordList.add(record);
                        this.asmNonsignificantPeak.put(geneId, recordList);
                    }
                }
            }
            bfr.close();
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
     * merge ASE and ASM test result get m6A peak with real allele specific
     */
    private void mergeAseAsmTestResult() {
        ArrayList<String[]> asmRecordList;
        String[] aseRecord, combineRecord;
        Set<String> aseGeneId = this.aseSignificantGene.keySet();
        Set<String> asmGeneId = this.asmSignificantPeak.keySet();
        boolean judge;
        // both ASE and ASM significant
        aseGeneId.retainAll(asmGeneId);
        for (String geneId: aseGeneId) {
            aseRecord = this.aseSignificantGene.get(geneId);
            asmRecordList = this.asmSignificantPeak.get(geneId);
            for (String[] asmRecord: asmRecordList) {
                judge = this.compareMajorAllele(aseRecord[aseRecord.length-1], asmRecord[asmRecord.length-1]);
                // chrNum, geneName, peakStart, peakEnd, ASE q-value, ASM q-value
                combineRecord = new String[] {asmRecord[0], aseRecord[0], asmRecord[1], asmRecord[2], aseRecord[1], asmRecord[3]};
                if (judge)
                    this.significantResult.put(geneId, combineRecord);
                else
                    this.nonSignificantResult.put(geneId, combineRecord);
            }
        }

        // ASE significant but ASM not
        aseGeneId = this.aseSignificantGene.keySet();
        Set<String> nonAsmGeneId = this.asmNonsignificantPeak.keySet();
        nonAsmGeneId.retainAll(aseGeneId);
        for (String geneId: nonAsmGeneId) {
            aseRecord = this.aseSignificantGene.get(geneId);
            asmRecordList = this.asmNonsignificantPeak.get(geneId);
            for (String[] asmRecord: asmRecordList) {
                combineRecord = new String[] {asmRecord[0], aseRecord[0], asmRecord[1], asmRecord[2], aseRecord[1], asmRecord[3]};
                this.nonSignificantResult.put(geneId, combineRecord);
            }
        }

        // ASM significant but ASE not
        Set<String> nonAseGeneId = this.aseNonsignificantGene.keySet();
        asmGeneId.retainAll(nonAseGeneId);
        for (String geneId: asmGeneId) {
            aseRecord = this.aseNonsignificantGene.get(geneId);
            asmRecordList = this.asmSignificantPeak.get(geneId);
            for (String[] asmRecord: asmRecordList) {
                combineRecord = new String[] {asmRecord[0], aseRecord[0], asmRecord[1], asmRecord[2], aseRecord[1], asmRecord[3]};
                this.significantResult.put(geneId, combineRecord);
            }
        }

        // both insignificant
        nonAseGeneId = this.aseNonsignificantGene.keySet();
        nonAsmGeneId = this.asmNonsignificantPeak.keySet();
        nonAseGeneId.retainAll(nonAsmGeneId);
        for (String geneId: nonAseGeneId) {
            aseRecord = this.aseNonsignificantGene.get(geneId);
            asmRecordList = this.asmNonsignificantPeak.get(geneId);
            for (String[] asmRecord: asmRecordList) {
                combineRecord = new String[] {asmRecord[0], aseRecord[0], asmRecord[1], asmRecord[2], aseRecord[1], asmRecord[3]};
                this.nonSignificantResult.put(geneId, combineRecord);
            }
        }
    }

    /**
     * when both significant in ASE and ASM test, judge whether the major allele is difference between gene and peak
     * @param aseMajorAllele major allele on gene
     * @param asmMajorAllele major allele on peak
     * @return true, if difference exists; otherwise false
     */
    private boolean compareMajorAllele(String aseMajorAllele, String asmMajorAllele) {
        String[] aseInfo = aseMajorAllele.split(";");
        String[] asmInfo = asmMajorAllele.split(";");
        HashMap<String, String> aseMajorNC = this.formMap(aseInfo);
        HashMap<String, String> asmMajorNC = this.formMap(asmInfo);
        Set<String> aseSNVSite = aseMajorNC.keySet();
        Set<String> asmSNVSite = asmMajorNC.keySet();

        // common sites
        aseSNVSite.retainAll(asmSNVSite);
        for (String site: aseSNVSite) {
            if (!aseMajorNC.get(site).equals(asmMajorNC.get(site)))
                return true;
        }

        return false;
    }

    private HashMap<String, String> formMap(String[] records) {
        HashMap<String, String> result = new HashMap<>();
        for (String rec: records) {
            String[] info = rec.split(":");
            result.put(info[0], info[1]);
        }

        return result;
    }

    /**
     * output final result
     */
    private void outputResult() {
        BufferedWriter bfw = null;
        try {
            bfw = new BufferedWriter(
                    new OutputStreamWriter(new FileOutputStream(this.outputFile))
            );

            // sort by ASE q value
            ArrayList<Map.Entry<String, String[]>> sigResult = new ArrayList<>(this.significantResult.entrySet());
            Collections.sort(sigResult, new Comparator<Map.Entry<String, String[]>>() {
                public int compare(Map.Entry<String, String[]> o1,
                                   Map.Entry<String, String[]> o2) {
                    Double val1 = Double.parseDouble(o1.getValue()[4]);
                    Double val2 = Double.parseDouble(o2.getValue()[4]);
                    return val1.compareTo(val2);
                }
            });

            ArrayList<Map.Entry<String, String[]>> nonsigResult = new ArrayList<>(this.nonSignificantResult.entrySet());
            Collections.sort(nonsigResult, new Comparator<Map.Entry<String, String[]>>() {
                public int compare(Map.Entry<String, String[]> o1,
                                   Map.Entry<String, String[]> o2) {
                    Double val1 = Double.parseDouble(o1.getValue()[4]);
                    Double val2 = Double.parseDouble(o2.getValue()[4]);
                    return val1.compareTo(val2);
                }
            });

            String geneId, chrNum, geneName, peakStart, peakEnd, aseQVal, asmQVal, line;
            String[] record;
            bfw.write("#chr\tgeneName\tgeneId\tpeakStart\tpeakEnd\tASE q-value\tASM q-value\tsignificant\n");
            for (Map.Entry<String, String[]> result: sigResult) {
                geneId = result.getKey();
                // chrNum, geneName, peakStart, peakEnd, ASE q-value, ASM q-value
                record = result.getValue();
                chrNum = record[0];
                geneName = record[1];
                peakStart = record[2];
                peakEnd = record[3];
                aseQVal = record[4];
                asmQVal = record[5];
                line = String.join("\t", new String[] {chrNum, geneName, geneId, peakStart, peakEnd, aseQVal,
                                                                  asmQVal, String.valueOf(true)});
                bfw.write(line);
                bfw.newLine();
            }
            for (Map.Entry<String, String[]> result: nonsigResult) {
                geneId = result.getKey();
                // chrNum, geneName, peakStart, peakEnd, ASE q-value, ASM q-value
                record = result.getValue();
                chrNum = record[0];
                geneName = record[1];
                peakStart = record[2];
                peakEnd = record[3];
                aseQVal = record[4];
                asmQVal = record[5];
                line = String.join("\t", new String[] {chrNum, geneName, geneId, peakStart, peakEnd, aseQVal,
                                                                 asmQVal, String.valueOf(false)});
                bfw.write(line);
                bfw.newLine();
            }
        } catch (IOException ie) {
            this.log.error("can not write ASE m6a peak result file");
            this.log.error(ie.getMessage());
        } finally {
            if (bfw != null) {
                try {
                    bfw.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    private static Logger initLog(String logHome) {
        System.setProperty("log_home", logHome);
        return Logger.getLogger(AlleleSpecificM6aSignal.class);
    }

    private static CommandLine setCommandLine(String[] args, Options options) throws ParseException {
        Option option = new Option("ase", "ase_gene_result", true, "ASE gene test result");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("asm", "asm_peak_result", true, "ASM m6A peak test result");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("o", "output_file", true, "output file path, default allele_specific_m6A_peak.txt");
        option.setRequired(false);
        options.addOption(option);

        CommandLineParser parser = new DefaultParser();

        return parser.parse(options, args);
    }
}
