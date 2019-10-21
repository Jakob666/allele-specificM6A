package HierarchicalBayesianAnalysis;

import org.apache.commons.cli.*;
import org.apache.log4j.Logger;

import java.io.File;

public class HierarchicalTest {
    private String aseVcfFile, asmVcfFile, wesFile, gtfFile, bedFile, peakCoveredSnpFile, peakCoveredSnpBackground, aseGeneFile, asmPeakFile, finalOutput;
    private int ipSNPReadInfimum, wesSNPReadInfimum, readsCoverageThreshold, samplingTime, burnInTime;
    private double degreeOfFreedom;
    private Logger logger;

    public static void main(String[] args) throws ParseException {
        Options options = new Options();
        CommandLine commandLine = setCommandLine(args, options);

        String aseVcfFile = null, asmVcfFile = null, wesFile = null, gtfFile = null, bedFile = null, outputDir;
        int ipSNPReadInfimum = 0, wesSNPReadInfimum = 0, readsCoverageThreshold = 10, samplingTime = 5000, burnInTime = 200;
        double tauInfimum = 0, tauSupremum = 2;

        if (commandLine.hasOption("o")) {
            outputDir = commandLine.getOptionValue("o");
            File output = new File(outputDir);
            if (output.exists() && output.isFile()) {
                System.out.println("-o parameter should be a path of a directory, not a file");
                System.exit(2);
            } else if (!output.exists()) {
                boolean res = output.mkdirs();
                if (!res) {
                    System.out.println("make directory " + outputDir + "failed");
                    System.exit(2);
                }
            }
        } else {
            File output = new File(System.getProperty("user.dir"), "outputResult");
            if (!output.exists()) {
                boolean res = output.mkdir();
                if (!res) {
                    System.out.println("make directory " + output +" failed");
                    System.exit(2);
                }
            }
            outputDir = output.getAbsolutePath();
        }
        Logger logger = initLog(outputDir);

        if (commandLine.hasOption("ase"))
            aseVcfFile = commandLine.getOptionValue("ase");
        if (aseVcfFile != null && !new File(aseVcfFile).exists()) {
            logger.error("file " + aseVcfFile + "doesn't exist");
            System.exit(2);
        } if (aseVcfFile != null && !new File(aseVcfFile).isFile()) {
            logger.error(aseVcfFile + " is not a file");
            System.exit(2);
        }

        if (commandLine.hasOption("asm"))
            asmVcfFile = commandLine.getOptionValue("ase");
        if (asmVcfFile != null && !new File(asmVcfFile).exists()) {
            logger.error("file " + asmVcfFile + "doesn't exist");
            System.exit(2);
        } if (asmVcfFile != null && !new File(asmVcfFile).isFile()) {
            logger.error(asmVcfFile + " is not a file");
            System.exit(2);
        }

        if (aseVcfFile == null || asmVcfFile == null) {
            logger.error("both INPUT and IP sample SNP calling VCF result files are needed");
            System.exit(2);
        }

        if (commandLine.hasOption("wes"))
            wesFile = commandLine.getOptionValue("wes");
        if (commandLine.hasOption("g"))
            gtfFile = commandLine.getOptionValue("g");
        if (commandLine.hasOption("b"))
            bedFile = commandLine.getOptionValue("b");
        if (gtfFile == null) {
            logger.error("GTF annotation file is required");
            System.exit(2);
        }
        if (bedFile == null) {
            logger.error("peak calling BED format file is required");
            System.exit(2);
        }
        if (commandLine.hasOption("tl"))
            tauInfimum = Double.parseDouble(commandLine.getOptionValue("tl"));
        if (commandLine.hasOption("th"))
            tauSupremum = Double.parseDouble(commandLine.getOptionValue("th"));
        if (tauInfimum >= tauSupremum) {
            logger.error("invalid uniform distribution parameter for tau sampling.");
            System.exit(2);
        }
        if (commandLine.hasOption("ip_cov_infimum"))
            ipSNPReadInfimum = Integer.parseInt(commandLine.getOptionValue("ip_cov_infimum"));
        if (commandLine.hasOption("wes_cov_infimum"))
            wesSNPReadInfimum = Integer.parseInt(commandLine.getOptionValue("wes_cov_infimum"));
        if (commandLine.hasOption("rc"))
            readsCoverageThreshold = Integer.valueOf(commandLine.getOptionValue("rc"));
        if (commandLine.hasOption("st"))
            samplingTime = Integer.parseInt(commandLine.getOptionValue("st"));
        if (commandLine.hasOption("bt"))
            burnInTime = Integer.parseInt(commandLine.getOptionValue("bt"));

        HierarchicalTest ht = new HierarchicalTest(aseVcfFile, asmVcfFile, wesFile, gtfFile, bedFile, outputDir,
                                                   tauSupremum, ipSNPReadInfimum, wesSNPReadInfimum,
                                                   readsCoverageThreshold, samplingTime, burnInTime, logger);
        ht.getResult();
    }


    public HierarchicalTest(String aseVcfFile, String asmVcfFile, String wesFile, String gtfFile, String bedFile, String outputDir,
                            double degreeOfFreedom, int ipSNPReadInfimum, int wesSNPReadInfimum,
                            int readsCoverageThreshold, int samplingTime, int burnInTime, Logger logger) {
        this.aseVcfFile = aseVcfFile;
        this.asmVcfFile = asmVcfFile;
        this.wesFile = wesFile;
        this.gtfFile = gtfFile;
        this.bedFile = bedFile;
        this.degreeOfFreedom = degreeOfFreedom;
        this.ipSNPReadInfimum = ipSNPReadInfimum;
        this.wesSNPReadInfimum = wesSNPReadInfimum;
        this.readsCoverageThreshold = readsCoverageThreshold;
        this.samplingTime = samplingTime;
        this.burnInTime = burnInTime;
        this.logger = logger;
        this.peakCoveredSnpFile = new File(outputDir, "asm_peakCoveredSNP.txt").getAbsolutePath();
        if (this.wesFile != null)
            this.peakCoveredSnpBackground = new File(outputDir, "asm_peakCoveredSNPBackground.txt").getAbsolutePath();
        else
            this.peakCoveredSnpBackground = null;
        this.asmPeakFile = new File(outputDir, "asmPeak.txt").getAbsolutePath();
        this.aseGeneFile = new File(outputDir, "aseGene.txt").getAbsolutePath();
        this.finalOutput = new File(outputDir, "finalResult.txt").getAbsolutePath();
    }

    public void getResult() {
        this.aseGeneDetected();
        this.asmPeakDetected();
        this.finalOutput();
    }

    /**
     * 检验m6A信号的ASM显著性
     */
    private void asmPeakDetected() {
        this.logger.debug("detect ASM m6A signal, m6A coveredSNP sites are shown in " + this.peakCoveredSnpFile);
        AsmPeakDetection apd = new AsmPeakDetection(this.bedFile, this.asmVcfFile,this.wesFile, null, peakCoveredSnpFile,
                                                    this.peakCoveredSnpBackground, this.asmPeakFile,
                                                    this.degreeOfFreedom, this.ipSNPReadInfimum, this.wesSNPReadInfimum,
                                                    this.samplingTime, this.burnInTime);
        apd.getTestResult();
        this.logger.debug("Hierarchical test result output in " + this.asmPeakFile + ", ASM specific m6A signals with q-value less than 0.05");
    }

    /**
     * 检验等位基因的ASE显著性
     */
    private void aseGeneDetected() {
        this.logger.debug("detect ASE Gene");
        AseGeneDetection agd = new AseGeneDetection(this.gtfFile, this.aseVcfFile, this.wesFile, null, this.aseGeneFile,
                                                    this.degreeOfFreedom, this.readsCoverageThreshold,
                                                    this.wesSNPReadInfimum, this.samplingTime, this.burnInTime);
        agd.getTestResult();
        this.logger.debug("Hierarchical test result output in " + this.aseGeneFile + ", ASE specific Genes with q-value less than 0.05");
    }

    /**
     * 依据ASE和ASM的检验结果，将IP样本的m6A信号分为4类
     */
    private void finalOutput() {
        AseSpecificM6aSignal asms = new AseSpecificM6aSignal(this.aseGeneFile, this.asmPeakFile, this.finalOutput, this.logger);
        asms.detect();
    }

    private static Logger initLog(String logHome) {
        System.setProperty("log_home", logHome);
        return Logger.getLogger(HierarchicalTest.class);
    }

    private static CommandLine setCommandLine(String[] args, Options options) throws ParseException {
        Option option = new Option("ase", "ase_vcf_file", true, "INPUT sample SNP calling result VCF file");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("asm", "asm_vcf_file", true, "IP sample SNP calling result VCF file");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("wes", "wes_vcf_file", true, "WES SNP calling result VCF file");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("g", "gft_file", true,"GTF annotation file");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("b", "bed_file", true, "Peak calling BED format file");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("o", "output_dir", true, "result output directory");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("df", "degree_of_freedom", true, "degree of freedom of inverse-Chi-square distribution, default 10");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("ip_cov_infimum", "ip_snp_coverage_infimum", true, "IP sample SNP site coverage infimum, default 0");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("wes_cov_infimum", "wes_snp_coverage_infimum", true, "WES sample SNP site coverage infimum, default 0");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("rc", "readsCoverage", true, "reads coverage threshold using for filter SNV records, default 10");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("st", "sampling_time", true, "sampling time of MH sampling and Gibbs sampling");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("bt", "burn_in", true, "burn-in time of MH sampling and Gibbs sampling");
        option.setRequired(false);
        options.addOption(option);

        CommandLineParser parser = new DefaultParser();

        return parser.parse(options, args);
    }
}
