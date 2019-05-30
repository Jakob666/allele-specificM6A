package DnaSeqSnpCalling;

import RnaSeqSnpCalling.RnaSeqSnpCalling;
import org.apache.commons.cli.*;
import org.apache.log4j.Logger;
import ReadsMapping.sra2fastq;
import CommonThreadClass.DnaSeqReadsMapping;
import CommonThreadClass.RunSnpCalling;

import java.io.File;
import java.io.IOException;

public class DnaSeqSnpCalling {

    public static void main(String[] args) throws ParseException{
        Options options = new Options();
        CommandLine commandLine = setCommand(args, options);

        String ipDataDir, inputDataDir, genomeFile, fastqTempDir=System.getProperty("user.dir"), outputDir, prefix;
        String bwa = "bwa";
        String gatkJar = "gatk";
        String picardJar = "picard";
        String samtools = "samtools";
        String bcftools = "bcftools";
        String inputFormat = "sra";
        int execThread = 2;
        File[] ipFastqFiles, inputFastqFiles;

        genomeFile = commandLine.getOptionValue('r');
        String genomeFileDir = new File(genomeFile).getParent();

        // output result directory, default a new directory name "outputResult" in genome file directory
        if (commandLine.hasOption('o')) {
            outputDir = commandLine.getOptionValue('o');
        } else {
            outputDir = new File(genomeFileDir, "outputResult").getAbsolutePath();
            mkDir(outputDir);
        }

        Logger logger = initLog(outputDir);

        // source data directory
        ipDataDir = commandLine.getOptionValue("ip");
        inputDataDir = commandLine.getOptionValue("input");
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
        if (commandLine.hasOption("bwa")) {
            bwa = commandLine.getOptionValue("bwa");
        }
        if (commandLine.hasOption('t')) {
            execThread = Integer.parseInt(commandLine.getOptionValue('t'));
        }
        if (commandLine.hasOption("fmt"))
            inputFormat = commandLine.getOptionValue("fmt");

        // 创建INPUT样本输出路径，如果是SRA格式则需对其进行转换
        String inputOutputDir = new File(outputDir, "INPUT").getAbsolutePath();
        mkDir(inputOutputDir);
        if (inputFormat.toLowerCase().equals("sra")) {
            // 获取到转换为fastq格式的数据
            inputFastqFiles = sraTransformToFastq(inputDataDir, fastqTempDir, logger);
        } else {
            inputFastqFiles = new File(inputDataDir).listFiles();
        }
        // 对INPUT样本进行处理，得到比对的bam文件，首先需要保证INPUT执行完，这一过程主函数处于阻塞状态
        DnaSeqReadsMapping dsrm = new DnaSeqReadsMapping(inputFastqFiles, genomeFile, inputOutputDir, prefix, bwa,
                                                         samtools, picardJar, gatkJar, execThread, logger);
        Thread inputThread = new Thread(dsrm);
        inputThread.start();
        try {
            inputThread.join();
        } catch (InterruptedException ie) {
            logger.error(ie.getMessage());
            System.exit(2);
        }

        // 执行完成后清除多余的文件
        if (inputFormat.toLowerCase().equals("sra")) {
            cleanUp(fastqTempDir, logger);
            mkDir(fastqTempDir);
        }
        String mergedBamFile = new File(inputOutputDir, prefix + "_alignment.bam").getAbsolutePath();

        // 创建IP样本输出路径，如果是SRA格式则需对其进行转换
        String ipOutputDir = new File(outputDir, "IP").getAbsolutePath();
        mkDir(ipOutputDir);
        if (inputFormat.toLowerCase().equals("sra")) {
            // 获取到转换为fastq格式的数据
            ipFastqFiles = sraTransformToFastq(ipDataDir, fastqTempDir, logger);
        } else {
            ipFastqFiles = new File(ipDataDir).listFiles();
        }
        dsrm = new DnaSeqReadsMapping(ipFastqFiles, genomeFile, ipOutputDir, prefix, bwa, samtools, picardJar, gatkJar,
                                      execThread, logger);
        Thread ipThread = new Thread(dsrm);

        // 对INPUT样本进行SNP calling
        RunSnpCalling rsc = new RunSnpCalling(genomeFile, mergedBamFile, samtools, bcftools, logger);
        Thread snpCallingThread = new Thread(rsc);

        ipThread.start();
        snpCallingThread.start();
        try {
            ipThread.join();
            snpCallingThread.join();
        } catch (InterruptedException ie) {
            logger.error(ie.getMessage());
            System.exit(2);
        }

        // 执行完成后清除多余的文件
        if (inputFormat.toLowerCase().equals("sra")) {
            cleanUp(fastqTempDir, logger);
        }
    }

    private static CommandLine setCommand(String[] args, Options options) throws ParseException{

        Option option = new Option("r", "reference_genome", true, "reference genome file path");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("ip", "ip_files", true, "ip source data files, separate by ',' if exists multi-files");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("input", "input_files", true, "input source data directory, separate by ',' if exists multi-files");
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

        option = new Option("bwa", "bwatool", true, "BWA executive file path, default bwa");
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

        CommandLineParser parser = new DefaultParser();

        return parser.parse(options, args);
    }

    private static void mkDir(String dirName) {
        File targetDir = new File(dirName);
        if (!targetDir.exists()) {
            boolean res = targetDir.mkdir();
        }
        targetDir = null;
    }

    private static void cleanUp(String dirName, Logger logger) {
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

    private static Logger initLog(String logHome) {
        System.setProperty("log_home", logHome);
        return Logger.getLogger(RnaSeqSnpCalling.class);
    }

    private static File[] sraTransformToFastq(String sraFilesDir, String fastqTempDir, Logger logger) {
        sra2fastq sf = new sra2fastq(sraFilesDir, fastqTempDir, logger);
        boolean successTrans = sf.transFormat();
        if (!successTrans) {
            logger.error("sra to fastq failed");
            System.exit(2);
        }
        return new File(fastqTempDir).listFiles();
    }

}
