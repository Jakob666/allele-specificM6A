package RnaSeqSnpCalling;

import ReadsMapping.ReadsMapping;
import ReadsMapping.sra2fastq;
import SamtoolsPileupSNPCalling.SamtoolsPileupSNPCalling;
import meripSeqPeakCalling.peakCaller;
import org.apache.commons.cli.*;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;

public class RnaSeqSnpCalling {

    private static Logger logger;

    public static void main(String[] args) throws ParseException {
        Options options = new Options();
        CommandLine commandLine = setCommand(args, options);

        // initialize the default params
        String ipDataDir = null, inputDataDir = null, genomeFile, gtfFile, fastqTempDir=null, outputDir, prefix;
        String gatkJar = "gatk";
        String picardJar = "picard";
        String samtools = "samtools";
        String bcftools = "bcftools";
        String inputFormat = "sra";
        int execThread = 2;

        if (commandLine.hasOption('h')) {
            new HelpFormatter().printHelp("flume-ng agent", options, true);
            return;
        }
        // get reference genome file and GTF file
        genomeFile = commandLine.getOptionValue('r');
        gtfFile = commandLine.getOptionValue('g');
        String genomeFileDir = new File(genomeFile).getParent();

        // output result directory, default a new directory name "outputResult" in genome file directory
        if (commandLine.hasOption('o')) {
            outputDir = commandLine.getOptionValue('o');
        } else {
            outputDir = new File(genomeFileDir, "outputResult").getAbsolutePath();
            mkDir(outputDir);
        }
        logger = initLog(outputDir);

        // get input file format(support sra, fastq, fasta)
        if (commandLine.hasOption("fmt")) {
            inputFormat = commandLine.getOptionValue("fmt").toLowerCase();
            if (!inputFormat.equals("sra") && !inputFormat.equals("fastq") && !inputFormat.equals("fasta")) {
                logger.error("invalid input file format, the format should be sra or fastq");
                System.exit(2);
            }
        }

        // source data directory
        if (commandLine.hasOption("ip")) {
            ipDataDir = commandLine.getOptionValue("ip");
        }
        if (commandLine.hasOption("input")) {
            inputDataDir = commandLine.getOptionValue("input");
        }
        if (ipDataDir == null & inputDataDir == null) {
            logger.error("Specifies the value of at least one of the parameters -ip(--ip_dir) and -input(--input_dir)");
            System.exit(2);
        }

        if (inputFormat.equals("sra")) {
            if (commandLine.hasOption("tmp")) {
                fastqTempDir = commandLine.getOptionValue("tmp");
            } else {
                fastqTempDir = new File(genomeFileDir, "temp").getAbsolutePath();
                mkDir(fastqTempDir);
            }
        } else {
            fastqTempDir = ipDataDir;
        }

        if (commandLine.hasOption('p')) {
            prefix = commandLine.getOptionValue('p');
        } else {
            prefix = "AseM6a";
        }
        if (commandLine.hasOption("smt")) {
            samtools = commandLine.getOptionValue("smt");
        }
        if (commandLine.hasOption("bft")) {
            bcftools = commandLine.getOptionValue("bft");
        }
        if (commandLine.hasOption("gatk")) {
            gatkJar = commandLine.getOptionValue("gatk");
        }
        if (commandLine.hasOption("picard")) {
            picardJar = commandLine.getOptionValue("picard");
        }
        if (commandLine.hasOption('t')) {
            execThread = Integer.parseInt(commandLine.getOptionValue('t'));
        }


        String inputOutputDir = new File(outputDir, "INPUT").getAbsolutePath();
        mkDir(inputOutputDir);
        File[] fastqFiles = getFastqFiles(inputDataDir, inputFormat, fastqTempDir);
        String inputBamFile = readsMapping(fastqFiles, genomeFile, gtfFile, picardJar, gatkJar, samtools,
                inputOutputDir, prefix, execThread);
        if (inputFormat.equals("sra"))
            cleanUp(fastqTempDir);
        // 创建线程用于SNP calling
        RunSnpCalling runSnp = new RunSnpCalling(genomeFile, inputBamFile, samtools, bcftools, logger);
        Thread snpThread = new Thread(runSnp);
        snpThread.start();

        // make directories for fastq files and alignment result
        String ipOutputDir = new File(outputDir, "IP").getAbsolutePath();
        mkDir(ipOutputDir);
        fastqFiles = getFastqFiles(ipDataDir, inputFormat, fastqTempDir);
        if (fastqFiles == null) {
            logger.error("empty fastq data directory");
            System.exit(2);
        }
        String ipBamFile = readsMapping(fastqFiles, genomeFile, gtfFile, picardJar, gatkJar, samtools,
                ipOutputDir, prefix, execThread);
        if (inputFormat.equals("sra"))
            cleanUp(fastqTempDir);




//        HashMap<String, String> alignRes = alignment(ipDataDir, inputDataDir, inputFormat, fastqTempDir, genomeFile,
//                                                     gtfFile, gatkJar, picardJar, samtools, outputDir, prefix, execThread);
//
//        String ipBamFile = alignRes.get("ip");
//        String inputBamFile = alignRes.get("input");

        // 创建两个线程，分别用于执行SNP calling和peak calling
        RunPeakCalling rpc = new RunPeakCalling(ipBamFile, inputBamFile, gtfFile, outputDir, logger);
        Thread peakThread = new Thread(rpc);
        peakThread.start();
        // 主线程等待两个子线程运行完才结束
        try {
            snpThread.join();
            peakThread.join();
        } catch (InterruptedException ie) {
            logger.error(ie.getMessage());
            System.exit(2);
        }

        if (inputFormat.equals("sra")) {
            try {
                Process p = Runtime.getRuntime().exec("rm -rf " + fastqTempDir);
                int exitVal = p.waitFor();
                if (exitVal != 0) {
                    logger.error("remove redundant fastq file failed");
                    System.exit(2);
                }
            } catch (IOException | InterruptedException ie) {
                logger.error(ie.getMessage());
            }
        }
    }

    /**
     * Deprecated!
     * 比对过程+预处理
     * @param inputFormat input data format, 'sra', 'fastq' or 'fasta'
     * @param fastqTempDir fastq data directory, if input format is sra, this directory is a temporary directory
     * @param refGenomeFile reference genome file
     * @param gtfFile reference GTF file
     * @param gatkJar gatk executive jar
     * @param picardJar picard executive jar
     * @param samtools samtools executive file
     * @param outputDir output directory
     * @param prefix prefix of output result file
     * @param execThread number of working thread
     * @return IP alignment bam file
     */
    private static HashMap<String, String> alignment(String ipDataDir, String inputDataDir, String inputFormat,
                                                        String fastqTempDir,String refGenomeFile, String gtfFile,
                                                        String gatkJar, String picardJar, String samtools, String outputDir,
                                                        String prefix, int execThread) {
        String ipOutputDir = new File(outputDir, "IP").getAbsolutePath();
        mkDir(ipOutputDir);
//        String ipPrefix = prefix + "_ip";
        File[] fastqFiles = getFastqFiles(ipDataDir, inputFormat, fastqTempDir);

        if (fastqFiles == null) {
            logger.error("empty fastq data directory");
            System.exit(2);
        }
        String ipBamFile = readsMapping(fastqFiles, refGenomeFile, gtfFile, picardJar, gatkJar, samtools,
                                        ipOutputDir, prefix, execThread);
        if (inputFormat.equals("sra"))
            cleanUp(fastqTempDir);


        String inputOutputDir = new File(outputDir, "INPUT").getAbsolutePath();
        mkDir(inputOutputDir);
//        String inputPrefix = prefix + "_input";
        fastqFiles = getFastqFiles(inputDataDir, inputFormat, fastqTempDir);
        String inputBamFile = readsMapping(fastqFiles, refGenomeFile, gtfFile, picardJar, gatkJar, samtools,
                                           inputOutputDir, prefix, execThread);
        if (inputFormat.equals("sra"))
            cleanUp(fastqTempDir);

        HashMap<String, String> map = new HashMap<>();
        map.put("ip", ipBamFile);
        map.put("input", inputBamFile);

        return map;
    }

    /**
     * 将sra格式转化为fastq格式
     * @param sraInputDir the directory sra format files, it should contain two sub-directory "IP" and "INPUT"
     * @param fastqOutputDir the tmp output fastq file directory
     * @return true if trans successfully otherwise false
     */
    private static boolean sraToFastq(String sraInputDir, String fastqOutputDir) {
        boolean execSuccess;

        sra2fastq s2fq = new sra2fastq(sraInputDir, fastqOutputDir, logger);
        execSuccess = s2fq.transFormat();
        s2fq = null;

        return execSuccess;
    }

    /**
     * 获取到文件目录下的fastq文件
     * @param DataDir directory of sra or fastq data
     * @param inputFormat sra or fastq
     * @param fastqTempDir fastq data directory, if input format is sra, this directory is a temporary directory
     * @return fastq file list
     */
    private static File[] getFastqFiles(String DataDir, String inputFormat, String fastqTempDir) {
        File[] fastqFiles;
        if (inputFormat.equals("sra")) {
            mkDir(fastqTempDir);
            boolean sraTransRes = sraToFastq(DataDir, fastqTempDir);
            if (!sraTransRes) {
                logger.error("transform failed");
                System.exit(2);
            }
            fastqFiles = new File(fastqTempDir).listFiles();
        } else {
            fastqFiles = new File(DataDir).listFiles();
        }

        return fastqFiles;
    }

    /**
     * reads mapping
     * @param fastqFiles file list of fastq files
     * @param refGenomeFile reference genome file path
     * @param gtfFile reference GTF file path
     * @param picardJar executable picard jar
     * @param gatkJar executable gatk jar
     * @param samtools executable samtools jar
     * @param outputDir output directory
     * @param prefix output file prefix, the final alignment file names {prefix}_split.bam
     * @param execThread working thread
     * @return final alignment file path
     */
    private static String readsMapping(File[] fastqFiles, String refGenomeFile, String gtfFile, String picardJar,
                                       String gatkJar, String samtools, String outputDir, String prefix, int execThread) {
        String refGenomeDir = new File(refGenomeFile).getParent();
        int readLength = ReadsMapping.getFastqReadLength(fastqFiles[0].getAbsolutePath(), logger);
        ReadsMapping.readsMapping(refGenomeDir, refGenomeFile, gtfFile, fastqFiles, readLength, execThread, logger);

        // 返回最终生成的bam文件路径
        return ReadsMapping.filterMappedReads(refGenomeFile, picardJar, gatkJar, samtools, outputDir, prefix, logger);
    }

    /**
     * 如果目录不存在则创建
     * @param dirName directory name
     */
    private static void mkDir(String dirName) {
        File targetDir = new File(dirName);
        if (!targetDir.exists()) {
            boolean res = targetDir.mkdir();
        }
        targetDir = null;
    }

    /**
     * 清理临时目录中的fastq文件
     * @param dirName temporary directory name
     */
    private static void cleanUp(String dirName) {
        String cmd = "rm -rf " + dirName;
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitNum = p.waitFor();
            if (exitNum != 0) {
                logger.error("remove redundant fastq file failed");
            }
        } catch (IOException | InterruptedException ie) {
            logger.error(ie.getMessage());
        }
    }


    /**
     * 设置命令行参数
     * @param args arguments in command line
     * @param options Options instance
     * @return CommandLine instance
     * @throws ParseException
     */
    private static CommandLine setCommand(String[] args, Options options) throws ParseException{

        Option option = new Option("r", "ref_genome", true, "reference genome file absolute path");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("g", "gtf_file", true, "gene transfer format(GTF) file");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("ip", "ip_dir", true, "ip source data directory, ip data will be only used for reads mapping and peak calling");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("input", "input_dir", true, "input source data directory, input data will be used for reads mapping, peak calling and SNP calling");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("fmt", "input_format", true, "input file format, sra or fastq, default sra");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("p", "prefix", true, "output result file prefix, default 'AseM6a'");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("tmp", "temp_dir", true, "when input files in sra format, a temporary directory store fastq file generate by sra file, default PWD");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("o", "output", true, "absolute path of the output directory, default PWD");
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

        option = new Option("t", "threads", true, "number of working threads, default 2");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("h", "help", false, "display help text");
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
        return Logger.getLogger(RnaSeqSnpCalling.class);
    }
}

/**
 * 运行SNP calling的线程类
 */
class RunSnpCalling implements Runnable {
    private String refGenomeFile, inputBamFile, samtools, bcftools;
    private Logger logger;
    RunSnpCalling(String refGenomeFile, String inputBamFile, String samtools, String bcftools, Logger logger) {
        this.refGenomeFile = refGenomeFile;
        this.inputBamFile = inputBamFile;
        this.samtools = samtools;
        this.bcftools = bcftools;
        this.logger = logger;
    }
    public void run() {
        snpReadsCount();
    }

    /**
     * snp calling and filtration
     */
    private void snpReadsCount() {
        logger.debug("start SNP calling");
        String outputDir = new File(new File(inputBamFile).getParent()).getParent();
        SamtoolsPileupSNPCalling.snpCalling(refGenomeFile, outputDir, inputBamFile, samtools, bcftools, logger);
        logger.debug("SNP calling complete");
    }
}

/**
 * 运行peak calling的线程类
 */
class RunPeakCalling implements Runnable {
    private String ipBamFile, inputBamFile, gtfFile,  outputDir;
    private Logger logger;

    RunPeakCalling(String ipBamFile, String inputBamFile, String gtfFile, String outputDir, Logger logger) {
        this.inputBamFile = inputBamFile;
        this.ipBamFile = ipBamFile;
        this.gtfFile = gtfFile;
        this.outputDir = outputDir;
        this.logger = logger;
    }

    public void run() {
        getM6aPeak();
    }

    /**
     * get peak calling result
     */
    private void getM6aPeak() {
        String experimentName = "m6aPeak";
        peakCaller.peakCalling(gtfFile, outputDir, ipBamFile, inputBamFile, experimentName, logger);
    }
}
