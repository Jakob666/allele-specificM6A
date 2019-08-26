package AseM6aPeakDetector;

import CommonThreadClass.RnaSeqReadsMapping;
import CommonThreadClass.RunSnpCalling;
import ExomeSeqSnpCalling.ExomeSeqSnpCalling;
import HierarchicalBayesianAnalysis.HierarchicalTest;
import RnaSeqSnpCalling.RnaSeqSnpCalling;
import meripSeqPeakCalling.PeakCaller;
import org.apache.commons.cli.*;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.*;

public class AseM6aPeakDetector {
    private String ipParamVal, inputParamVal, wesParamVal, genomeFile, gtfFile, fastqTempDir, outputDir, direct, gatk,
                   picard, samtools, bcftools, bwa, inputFormat, wesInputFormat, experimentName;
    private int execThread, samplingTime, burnInTime;
    private boolean compress, wesCompress;
    private Logger log;

    public static void main(String[] args) throws ParseException{
        String ipParamVal = null, inputParamVal = null, wesParamVal = null, genomeFile = null, gtfFile = null, fastqTempDir = null,
                outputDir, direct = "se", gatk = "gatk", picard = "picard", samtools = "samtools", bcftools = "bcftools", bwa = "bwa",
                inputFormat = "sra", wesInputFormat = null, experimentName = "m6aPeak";
        boolean compress = false, wesCompress = false;
        int samplingTime = 5000, burnInTime = 200, execThread = 2;

        Options options = new Options();
        CommandLine commandLine = setCommand(args, options);

        if (commandLine.hasOption("r"))
            genomeFile = commandLine.getOptionValue("r");
        if (genomeFile == null) {
            System.out.println("reference genome file is required");
            System.exit(2);
        } else {
            File genome = new File(genomeFile);
            if (!genome.exists()) {
                System.out.println("file " + genomeFile + " doesn't exist");
                System.exit(2);
            } else if (!genome.isFile()) {
                System.out.println(genomeFile + " is not a file");
                System.exit(2);
            }
        }
        if (commandLine.hasOption("g"))
            gtfFile = commandLine.getOptionValue("g");
        if (gtfFile == null) {
            System.out.println("GTF annotation file is required");
            System.exit(2);
        } else {
            File gtf = new File(gtfFile);
            if (!gtf.exists()) {
                System.out.println("file " + gtfFile + " doesn't exist");
                System.exit(2);
            } else if (!gtf.isFile()) {
                System.out.println(gtfFile + " is not a file");
                System.exit(2);
            }
        }

        String genomeFileDir = new File(genomeFile).getParent();
        if (commandLine.hasOption('o'))
            outputDir = commandLine.getOptionValue('o');
        else
            outputDir = new File(genomeFileDir, "outputResult").getAbsolutePath();
        File output = new File(outputDir);
        if (!output.exists()) {
            boolean res = output.mkdirs();
            if (!res) {
                System.out.println("can not make directory " + outputDir);
                System.exit(2);
            }
        } else {
            if (output.isFile()) {
                boolean res = output.mkdirs();
                if (!res) {
                    System.out.println("can not make directory " + outputDir);
                    System.exit(2);
                }
            }
        }

        Logger logger = initLog(outputDir);

        if (commandLine.hasOption("ip"))
            ipParamVal = commandLine.getOptionValue("ip");
        if (ipParamVal == null) {
            logger.error("-ip parameter is required for peak calling");
            System.exit(2);
        }
        if (commandLine.hasOption("input"))
            inputParamVal = commandLine.getOptionValue("input");
        if (inputParamVal == null) {
            logger.error("-input parameter is required for peak calling");
            System.exit(2);
        }
        if (commandLine.hasOption("fmt"))
            inputFormat = commandLine.getOptionValue("fmt").toLowerCase();
        if (!inputFormat.equals("sra") && !inputFormat.equals("fastq") && !inputFormat.equals("bam")) {
            logger.error("invalid value of -fmt parameter, should be fastq, sra or bam");
            System.exit(2);
        }
        if (commandLine.hasOption("exp"))
            experimentName = commandLine.getOptionValue("exp");
        if (commandLine.hasOption("wes"))
            wesParamVal = commandLine.getOptionValue("wes");
        if (wesParamVal != null) {
            if (commandLine.hasOption("wes_fmt")) {
                wesInputFormat = commandLine.getOptionValue("wes_fmt").toLowerCase();
                if (!wesInputFormat.equals("sra") && !wesInputFormat.equals("fastq") && !wesInputFormat.equals("bam")) {
                    logger.error("invalid value of -wes_fmt parameter, should be fastq, sra or bam");
                    System.exit(2);
                }
            } else
                wesInputFormat = "sra";
        }

        if (commandLine.hasOption("c")) {
            String val = commandLine.getOptionValue("c");
            if (!val.equals("true") && !val.equals("false")) {
                logger.error("invalid value. -c parameter should be true or false");
                System.exit(2);
            }
            compress = Boolean.parseBoolean(val);
        }
        if (commandLine.hasOption("wc")) {
            String val = commandLine.getOptionValue("wc");
            if (!val.equals("true") && !val.equals("false")) {
                logger.error("invalid value. -wc parameter should be true or false");
                System.exit(2);
            }
            wesCompress = Boolean.parseBoolean(val);
        }
        if (inputFormat.equals("sra")) {
            if (commandLine.hasOption("tmp"))
                fastqTempDir = commandLine.getOptionValue("tmp");
        }
        if (fastqTempDir != null) {
            File tmp = new File(fastqTempDir);
            if (!tmp.exists()) {
                boolean res = tmp.mkdirs();
                if (!res) {
                    System.out.println("can not make directory " + tmp);
                    System.exit(2);
                }
            } else {
                if (tmp.isFile()) {
                    boolean res = tmp.mkdirs();
                    if (!res) {
                        System.out.println("can not make directory " + tmp);
                        System.exit(2);
                    }
                }
            }
        }

        if (commandLine.hasOption("d")) {
            direct = commandLine.getOptionValue("d");
            if (!direct.toLowerCase().equals("se") && !direct.toLowerCase().equals("pe")) {
                logger.error("invalid input sequencing direction, must be SE(single-end) or PE(pair-end)");
                System.exit(2);
            }
        }
        if (commandLine.hasOption("smt"))
            samtools = commandLine.getOptionValue("smt");
        if (commandLine.hasOption("bft"))
            bcftools = commandLine.getOptionValue("bft");
        if (commandLine.hasOption("gatk"))
            gatk = commandLine.getOptionValue("gatk");
        if (commandLine.hasOption("picard"))
            picard = commandLine.getOptionValue("picard");
        if (commandLine.hasOption("bwa"))
            bwa = commandLine.getOptionValue("bwa");
        if (commandLine.hasOption('t'))
            execThread = Integer.parseInt(commandLine.getOptionValue('t'));
        if (commandLine.hasOption("st"))
            samplingTime = Integer.parseInt(commandLine.getOptionValue("st"));
        if (commandLine.hasOption("bt"))
            burnInTime = Integer.parseInt(commandLine.getOptionValue("bt"));

        AseM6aPeakDetector ampd = new AseM6aPeakDetector(ipParamVal, inputParamVal, wesParamVal, genomeFile, gtfFile,
                fastqTempDir, outputDir, direct, gatk, picard, samtools, bcftools, bwa, inputFormat, wesInputFormat,
                experimentName, execThread, samplingTime, burnInTime, compress, wesCompress, logger);
        ampd.getResult();

    }

    public AseM6aPeakDetector(String ipParamVal, String inputParamVal, String wesParamVal, String genomeFile, String gtfFile,
                              String fastqTempDir, String outputDir, String direct, String gatk, String picard, String samtools,
                              String bcftools, String bwa, String inputFormat, String wesInputFormat, String experimentName,
                              int execThread, int samplingTime, int burnInTime, boolean compress, boolean wesCompress, Logger logger) {
        this.ipParamVal = ipParamVal;
        this.inputParamVal = inputParamVal;
        this.wesParamVal = wesParamVal;
        this.genomeFile = new File(genomeFile).getAbsolutePath();
        this.gtfFile = new File(gtfFile).getAbsolutePath();
        this.fastqTempDir = fastqTempDir == null? null: new File(fastqTempDir).getAbsolutePath();
        this.outputDir = new File(outputDir).getAbsolutePath();
        this.direct = direct;
        this.gatk = gatk;
        this.picard = picard;
        this.samtools = samtools;
        this.bcftools = bcftools;
        this.bwa = bwa;
        this.inputFormat = inputFormat;
        this.wesInputFormat = wesInputFormat;
        this.experimentName = experimentName;
        this.execThread = execThread;
        this.samplingTime = samplingTime;
        this.burnInTime = burnInTime;
        this.compress = compress;
        this.wesCompress = wesCompress;
        this.log = logger;
    }

    private void getResult() {
        // 如果只给出了IP和INPUT的MeRIP-seq数据
        ArrayList<String> inputFastqFileNames, ipFastqFileNames;
        File[] inputFastqFiles, ipFastqFiles;
        RnaSeqReadsMapping inputRsrm, ipRsrm;
        String ipOutputDir = new File(this.outputDir, "IP").getAbsolutePath();
        String ipBamFile = new File(ipOutputDir, "asm_alignment.bam").getAbsolutePath();
        String inputOutputDir = new File(this.outputDir, "INPUT").getAbsolutePath();
        String inputBamFile = new File(inputOutputDir, "ase_alignment.bam").getAbsolutePath();
        String aseSnpFile = new File(this.outputDir, "ase_filtered.vcf").getAbsolutePath();
        String asmSnpFile = new File(this.outputDir, "asm_filtered.vcf").getAbsolutePath();
        String wesSnpFile = null;
        String peakBedDir = new File(this.outputDir, this.experimentName).getAbsolutePath();
        String peakBedFile = new File(peakBedDir, "peak.bed").getAbsolutePath();

        if (wesParamVal == null) { // 仅使用MeRIP-seq数据进行SNP calling和peak calling
            RnaSeqSnpCalling rssc = new RnaSeqSnpCalling(this.ipParamVal, this.inputParamVal, this.genomeFile, this.gtfFile,
                    this.fastqTempDir, this.outputDir, this.direct, this.gatk, this.picard, this.samtools, this.bcftools,
                    this.inputFormat, this.execThread, this.compress, this.log);
            rssc.callSnp();
            PeakCaller pc = new PeakCaller(this.gtfFile, ipBamFile, inputBamFile, this.outputDir, this.experimentName, this.log);
            pc.peakCalling();
        } else { // 如果同时提供MeRIP-seq数据，同时提供WES数据
            RnaSeqSnpCalling rssc = new RnaSeqSnpCalling(this.ipParamVal, this.inputParamVal, this.genomeFile, this.gtfFile,
                    this.fastqTempDir, this.outputDir, this.direct, this.gatk, this.picard, this.samtools, this.bcftools,
                    this.inputFormat, this.execThread, this.compress, this.log);
            ArrayList<File> ipFiles = rssc.getUserFiles(this.direct, this.ipParamVal, this.inputFormat);
            ArrayList<File> inputFiles = rssc.getUserFiles(this.direct, this.inputParamVal, this.inputFormat);

            // 如果是非BAM文件则先对数据进行比对
            if (!this.inputFormat.toLowerCase().equals("bam")) {
                if (inputFormat.toLowerCase().equals("sra")) {
                    inputFastqFileNames = rssc.transToFormat(inputFiles, this.fastqTempDir, this.direct);
                    ipFastqFileNames = rssc.transToFormat(ipFiles, this.fastqTempDir, this.direct);
                } else {
                    inputFastqFileNames = new ArrayList<>();
                    for (File f: inputFiles)
                        inputFastqFileNames.add(f.getAbsolutePath());
                    ipFastqFileNames = new ArrayList<>();
                    for (File f: ipFiles)
                        ipFastqFileNames.add(f.getAbsolutePath());
                }
                this.log.debug(inputFastqFileNames);
                this.log.debug(ipFastqFileNames);
                inputFastqFiles = new File[inputFastqFileNames.size()];
                for (int i=0; i< inputFastqFileNames.size(); i++) {
                    File f = new File(inputFastqFileNames.get(i));
                    inputFastqFiles[i] = f;
                }
                ipFastqFiles = new File[ipFastqFileNames.size()];
                for (int i=0; i< ipFastqFileNames.size(); i++) {
                    File f = new File(ipFastqFileNames.get(i));
                    ipFastqFiles[i] = f;
                }

                inputRsrm = new RnaSeqReadsMapping(inputFastqFiles, this.genomeFile, this.gtfFile, this.picard, this.gatk, this.samtools,
                                                   inputOutputDir, "ase", this.direct, this.execThread, this.compress, this.log);
                Thread inputThread = new Thread(inputRsrm);
                ipRsrm = new RnaSeqReadsMapping(ipFastqFiles, this.genomeFile, this.gtfFile, this.picard, this.gatk, this.samtools,
                                                inputOutputDir, "asm", this.direct, this.execThread, this.compress, this.log);
                Thread ipThread = new Thread(ipRsrm);
                inputThread.start();
                try {
                    inputThread.join();
                    ipThread.start();
                    ipThread.join();
                } catch (InterruptedException ie) {
                    this.log.error(ie.getMessage());
                    System.exit(2);
                }
            } else { // BAM格式文件如果有多个，则需要进行合并
                File mergedBamFile = new File(inputOutputDir, "ase_alignment.bam");
                rssc.mergeBamFiles(inputFiles, mergedBamFile.getAbsolutePath());
            }
            // MeRIP-seq的比对结果用于peak calling
            PeakCaller pc = new PeakCaller(this.gtfFile, ipBamFile, inputBamFile, this.outputDir, this.experimentName, this.log);
            pc.peakCalling();
            // 对IP样本比对结果进行SNP calling
            RunSnpCalling ipRunSnp = new RunSnpCalling(this.genomeFile, ipBamFile, this.samtools, this.bcftools, this.log);
            Thread ipSnpThread = new Thread(ipRunSnp);
            ipSnpThread.start();
            // 传入的WES数据用于INPUT样本的SNP calling得到ASE gene
            ExomeSeqSnpCalling essc = new ExomeSeqSnpCalling(this.wesParamVal, this.genomeFile, this.fastqTempDir, this.outputDir,
                    "ase", this.direct, this.bwa, "", this.picard, this.samtools, this.bcftools, this.wesInputFormat,
                    this.execThread, this.wesCompress);
            essc.callSnp();
            try {
                ipSnpThread.join();
            } catch (InterruptedException ie) {
                this.log.error(ie.getMessage());
                System.exit(2);
            }
        }
        HierarchicalTest ht = new HierarchicalTest(aseSnpFile, asmSnpFile, wesSnpFile, this.gtfFile, peakBedFile, this.outputDir,
                                                   this.samplingTime, this.burnInTime, this.log);
        ht.getResult();
    }

    private static CommandLine setCommand(String[] args, Options options) throws ParseException {
        Option option = new Option("r", "reference_genome", true, "reference genome file");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("g", "gtf_file", true, "GTF annotation file");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("ip", "ip_files", true,
                "when single-end sequencing parameter is the IP data directory or file paths separate by ','\n" +
                        "-ip IP_Data_Dir  or  -ip f1,f2..\n" +
                        "; when pair-end sequencing the format like\n -ip f1_1,f1_2;f2_1,f2-2;...");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("input", "input_files", true,
                "when single-end sequencing parameter is the INPUT data directory or file paths separate by ','\n" +
                        "-input IP_Data_Dir  or  -input f1,f2..\n" +
                        "; when pair-end sequencing the format like\n -input f1_1,f1_2;f2_1,f2-2;...");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("wes", "wes_data", true, "Whole exome sequencing for SNP calling.\n" +
                "When single-end sequencing, parameter value in 2 format: \n\t1) source exome-seq data directory: /path/to/source_data_dir" +
                "\n\t2) sequencing data separate by comma: file1,file2,...\nWhen pair-end sequencing, parameter value like: mate1_1,mate1_2;mate2_1,mate2_2;... set parameter -d or --direct to PE");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("fmt", "input_format", true, "IP and INPUT MeRIP-seq file format, SRA, FastQ or BAM, default SRA");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("wes_fmt", "wes_input_format", true, "WES input file format, SRA, FastQ or BAM, default SRA");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("c", "compress", true, "MeRIP-seq data in GZIP compressed format. default false");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("wc", "wes_compress", true, "WES data in GZIP compressed format. default false");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("d", "direct", true, "SE or PE, represent single-end and pair-end, default SE");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("tmp", "temp_dir", true, "when input files in sra format, a temporary directory store fastq file generate by sra file, default PWD");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("o", "output", true, "absolute path of the output directory, default PWD");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("exp", "experimentName", true, "experimentName, directory for storing peak calling result");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("smt", "samtool", true, "samtools executive file path, default samtools");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("bft", "bcftool", true, "bcftools executive file path, default bcftools");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("gatk", "gatk_tool", true, "GATK executive file path, default gatk");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("picard", "picard_tool", true, "PICARD executive file path, default picard.jar");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("bwa", "bwa_tool", true, "BWA alignment tool if there has WES data");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("t", "threads", true, "number of working threads, default 2");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("st", "sampling_time", true, "sampling times of MH sampling and Gibbs sampling");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("bt", "burn_in_time", true, "burn-in times of MH sampling and Gibbs sampling");
        option.setRequired(false);
        options.addOption(option);

        CommandLineParser parser = new DefaultParser();

        return parser.parse(options, args);
    }

    /**
     * 初始化log4j Logger 对象
     * @param logHome output directory of log file
     * @return Logger instance
     */
    private static Logger initLog(String logHome) {
        System.setProperty("log_home", logHome);
        return Logger.getLogger(AseM6aPeakDetector.class);
    }
}
