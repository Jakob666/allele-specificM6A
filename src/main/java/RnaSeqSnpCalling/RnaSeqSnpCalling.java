package RnaSeqSnpCalling;

import CommonThreadClass.*;
import ReadsMapping.sra2fastq;
import org.apache.commons.cli.*;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;

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
        if (commandLine.hasOption('o'))
            outputDir = commandLine.getOptionValue('o');
        else {
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
        if (commandLine.hasOption("ip"))
            ipDataDir = commandLine.getOptionValue("ip");
        if (commandLine.hasOption("input"))
            inputDataDir = commandLine.getOptionValue("input");
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

        // 开启INPUT样本的分析
        String inputOutputDir = new File(outputDir, "INPUT").getAbsolutePath();
        mkDir(inputOutputDir);
        File[] fastqFiles = getFastqFiles(inputDataDir, inputFormat, fastqTempDir);
        RnaSeqReadsMapping rsrm = new RnaSeqReadsMapping(fastqFiles, genomeFile, gtfFile, picardJar, gatkJar, samtools,
                                                         inputOutputDir, prefix, execThread, logger);
        Thread inputThread = new Thread(rsrm);
        inputThread.start();
        try {
            inputThread.join();
        } catch (InterruptedException ie) {
            logger.error(ie.getMessage());
            System.exit(2);
        }
        if (inputFormat.equals("sra")) {
            cleanUp(fastqTempDir);
            mkDir(fastqTempDir);
        }

        // 开启IP样本的分析
        String ipOutputDir = new File(outputDir, "IP").getAbsolutePath();
        mkDir(inputOutputDir);
        fastqFiles = getFastqFiles(ipDataDir, inputFormat, fastqTempDir);
        rsrm = new RnaSeqReadsMapping(fastqFiles, genomeFile, gtfFile, picardJar, gatkJar, samtools, ipOutputDir,
                                      prefix, execThread, logger);
        Thread ipThread = new Thread(rsrm);

        // 创建线程用于SNP calling
        String inputBamFile = new File(inputOutputDir, prefix+"_alignment.bam").getAbsolutePath();
        RunSnpCalling runSnp = new RunSnpCalling(genomeFile, inputBamFile, samtools, bcftools, logger);
        Thread snpThread = new Thread(runSnp);

        ipThread.start();
        snpThread.start();
        // 主线程等待两个子线程运行完才结束
        try {
            ipThread.join();
            snpThread.join();
        } catch (InterruptedException ie) {
            logger.error(ie.getMessage());
            System.exit(2);
        }

        if (inputFormat.toLowerCase().equals("sra"))
            cleanUp(fastqTempDir);
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

