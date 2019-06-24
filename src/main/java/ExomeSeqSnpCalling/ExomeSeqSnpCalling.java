package ExomeSeqSnpCalling;

import CommonThreadClass.RunSnpCalling;
import CommonThreadClass.Sra2Fastq;
import org.apache.commons.cli.*;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

public class ExomeSeqSnpCalling {
    private static Logger logger;

    public static void main(String[] args) throws ParseException {
        Options options = new Options();
        CommandLine commandLine = setCommand(args, options);
        String helpMessage = "";
        if (commandLine.hasOption("h")) {
            System.out.println(helpMessage);
            System.exit(0);
        }
        String source_dir, genomeFile, fastqTempDir, outputDir, prefix, direct = "SE";
        String bwa = "bwa";
        String sraTool = "";
        String picardJar = "picard";
        String samtools = "samtools";
        String bcftools = "bcftools";
        String inputFormat = "sra";
        int execThread = 2;

        // get reference genome file
        genomeFile = commandLine.getOptionValue('r');
        // source data directory
        source_dir = commandLine.getOptionValue("s");
        // result output directory
        if (commandLine.hasOption('o')) {
            outputDir = commandLine.getOptionValue('o');
            mkDir(outputDir);
        } else {
            outputDir = new File(System.getProperty("user.dir"), "outputResult").getAbsolutePath();
            mkDir(outputDir);
        }
        logger = initLog(outputDir);
        File alignmentOutputDir = new File(outputDir, "alignmentResult");
        mkDir(alignmentOutputDir.getAbsolutePath());

        if (commandLine.hasOption('p')) {
            prefix = commandLine.getOptionValue('p');
        } else {
            prefix = "AseM6a";
        }
        if (commandLine.hasOption("b")) {
            bwa = commandLine.getOptionValue("b");
        }
        if (commandLine.hasOption("smt")) {
            samtools = commandLine.getOptionValue("smt");
        }
        if (commandLine.hasOption("bft")) {
            bcftools = commandLine.getOptionValue("bft");
        }
        if (commandLine.hasOption("picard")) {
            picardJar = commandLine.getOptionValue("picard");
        }
        if (commandLine.hasOption('t')) {
            execThread = Integer.parseInt(commandLine.getOptionValue('t'));
        }
        // get input file format(support sra, fastq)
        if (commandLine.hasOption("fmt")) {
            inputFormat = commandLine.getOptionValue("fmt").toLowerCase();
            if (!inputFormat.equals("sra") && !inputFormat.equals("fastq") && !inputFormat.equals("bam")) {
                logger.error("invalid input file format, the format should be sra, fastq or bam");
                System.exit(2);
            }
        }

        if (commandLine.hasOption("d")) {
            direct = commandLine.getOptionValue("d");
            if (direct.toLowerCase().equals("se")) {
                direct = "SE";
            } else if (direct.toLowerCase().equals("pe")) {
                direct = "PE";
            } else {
                logger.error("invalid input of sequencing direction, should be SE or PE");
                System.out.println(helpMessage);
                System.exit(2);
            }
        }

        // 如果输入数据是SRA格式，则需转为fastq格式，确定fastq文件的临时存放路径
        if (inputFormat.toLowerCase().equals("sra")) {
            if (commandLine.hasOption("tmp")) {
                fastqTempDir = commandLine.getOptionValue("tmp");
            } else {
                fastqTempDir = new File(outputDir, "temp").getAbsolutePath();
                mkDir(fastqTempDir);
            }
        } else if (inputFormat.toLowerCase().equals("fastq")) { // 如果输入数据本身就是fastq格式，则无需转换
            fastqTempDir = source_dir;
        } else {
            fastqTempDir = "";
        }

        String mergedBamFile;
        if (!inputFormat.toLowerCase().equals("bam")) {
            // 对测序文件进行比对
            logger.debug("start exome sequencing reads alignment");

            mergedBamFile = alignSequencingReads(fastqTempDir, bwa, samtools, picardJar, sraTool, source_dir, inputFormat,
                    alignmentOutputDir.getAbsolutePath(), genomeFile, execThread, direct, prefix);
            if (inputFormat.toLowerCase().equals("sra"))
                cleanUp(fastqTempDir);
            logger.debug("Exome sequencing reads alignment complete");
        } else {
            File mergedBam = new File(outputDir, prefix+"_alignment.bam");
            mergedBamFile = mergedBam.getAbsolutePath();
            File[] bamFile = new File(source_dir).listFiles();
            if (bamFile == null) {
                logger.error("empty bam File data");
                System.exit(2);
            }
            ArrayList<String> alignmentBamFile = new ArrayList<>();
            for (File f: bamFile) {
                if (f.getAbsolutePath().endsWith("bam")) {
                    alignmentBamFile.add(f.getAbsolutePath());
                }
            }
            if (alignmentBamFile.size() == 0) {
                logger.error("directory contains no bam File");
                System.exit(2);
            }
            logger.debug("merge alignment bam files");
            if (alignmentBamFile.size() == 1) {
                boolean res = new File(alignmentBamFile.get(0)).renameTo(mergedBam);
                if (!res) {
                    logger.error("rename file error: " + alignmentBamFile.get(0));
                }
            } else {
                ReadsAlignment readsAlignment = new ReadsAlignment(bwa, samtools, picardJar, genomeFile, execThread, logger);
                readsAlignment.mergeBamFiles(alignmentBamFile, mergedBam.getAbsolutePath());
            }
        }

        // 对INPUT样本进行SNP calling
        logger.debug("start exome SNP calling");
        RunSnpCalling rsc = new RunSnpCalling(genomeFile, mergedBamFile, samtools, bcftools, logger);
        Thread snpCallingThread = new Thread(rsc);
        snpCallingThread.start();

        try {
            snpCallingThread.join();
        } catch (InterruptedException ie) {
            logger.error(ie.getMessage());
            System.exit(2);
        }
        logger.debug("SNP calling complete");
    }

    /**
     * 对FastQ文件进行比对
     * @return merged BAM file
     */
    private static String alignSequencingReads(String fastqTempDir, String bwa, String samtools, String picardJar, String sraTool,
                                               String source_dir, String inputFormat, String outputDir, String genomeFile,
                                               int execThread, String direct, String prefix) {
        ReadsAlignment readsAlignment = new ReadsAlignment(bwa, samtools, picardJar, genomeFile, execThread, logger);

        // 先对参考基因建立索引。之后对每个FastQ文件进行比对，比对后对重复样本进行合并
        logger.debug("build index for reference genome: " + genomeFile);
        readsAlignment.genomeIndex();
        ArrayList<String> alignmentBamFile = new ArrayList<>();
        if (direct.toLowerCase().equals("se")) {    // single-end测序
            if (inputFormat.toLowerCase().equals("sra")) {
                Sra2Fastq sra2Fastq = Sra2Fastq.getInstance(sraTool, source_dir, fastqTempDir, direct, logger);
                sra2Fastq.transform();
            }
            File[] fastqFiles = new File(fastqTempDir).listFiles();
            if (fastqFiles == null) {
                logger.error("empty FastQ directory: " + fastqTempDir);
                System.exit(2);
            }
            for (File inputFile: fastqFiles) {
                String inputFileName = inputFile.getName();
                String alignmentSamFile = new File(outputDir, inputFileName.substring(0, inputFileName.lastIndexOf("."))+".sam").getAbsolutePath();
                String bamFile = readsAlignment.alignmentToGenome(inputFile.getAbsolutePath(), null, alignmentSamFile);
                alignmentBamFile.add(bamFile);
            }
        } else {    // pair-end测序
            if (inputFormat.toLowerCase().equals("sra")) {
                Sra2Fastq sra2Fastq = Sra2Fastq.getInstance(sraTool, source_dir, fastqTempDir, direct, logger);
                ArrayList<String> sraFilePaths = sra2Fastq.getSraFiles();
                for (String sraFile: sraFilePaths) {
                    sra2Fastq.switchToFastq(sraFile);
                    File[] fastqFiles = new File(fastqTempDir).listFiles();
                    String inputFile1 = fastqFiles[0].getAbsolutePath();
                    String inputFile2 = fastqFiles[1].getAbsolutePath();
                    String alignmentSamFile = inputFile1.substring(0, inputFile1.lastIndexOf("_")) + ".sam";
                    String bamFile = readsAlignment.alignmentToGenome(inputFile1, inputFile2, alignmentSamFile);
                    alignmentBamFile.add(bamFile);
                    new File(inputFile1).delete();
                    new File(inputFile2).delete();
                }
            } else {
                File[] fastqFiles = new File(fastqTempDir).listFiles();
                if (fastqFiles == null) {
                    logger.error("empty FastQ directory: " + fastqTempDir);
                    System.exit(2);
                }
                if (fastqFiles.length % 2 != 0) {
                    logger.error("pair-end sequence files must be odd number");
                    System.exit(2);
                }
                for (int i = 0; i < fastqFiles.length; i=i+2) {
                    String inputFile1 = fastqFiles[i].getAbsolutePath();
                    String inputFile2 = fastqFiles[i+1].getAbsolutePath();
                    String alignmentSamFile;
                    // 确定SAM文件名
                    if (inputFile1.contains("_")) {
                        alignmentSamFile = inputFile1.substring(0, inputFile1.lastIndexOf("_")) + ".sam";
                    } else {
                        alignmentSamFile = inputFile1.substring(0, inputFile1.lastIndexOf(".")) + ".sam";
                    }
                    String bamFile = readsAlignment.alignmentToGenome(inputFile1, inputFile2, alignmentSamFile);
                    alignmentBamFile.add(bamFile);
                }
            }
        }
        logger.debug("merge alignment bam files");
        File mergedBam = new File(outputDir, prefix+"_alignment.bam");
        if (alignmentBamFile.size() == 1) {
            boolean res = new File(alignmentBamFile.get(0)).renameTo(mergedBam);
            if (!res) {
                logger.error("rename file error: " + alignmentBamFile.get(0));
            }
        } else {
            readsAlignment.mergeBamFiles(alignmentBamFile, mergedBam.getAbsolutePath());
        }

        return mergedBam.getAbsolutePath();
    }

    /**
     * 如果目录不存在则创建
     * @param dirName directory name
     */
    private static void mkDir(String dirName) {
        File targetDir = new File(dirName);
        if (!targetDir.exists()) {
            boolean res = targetDir.mkdir();
            if (!res) {
                logger.error("can not make directory: " + dirName);
                System.exit(2);
            }
        }
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
     * 初始化log4j Logger 对象
     * @param logHome output directory of log file
     * @return Logger instance
     */
    private static Logger initLog(String logHome) {
        System.setProperty("log_home", logHome);
        return Logger.getLogger(ExomeSeqSnpCalling.class);
    }

    /**
     * 设置命令行参数
     * @param args arguments in command line
     * @param options Options instance
     * @return CommandLine instance
     * @throws ParseException throw out exception
     */
    private static CommandLine setCommand(String[] args, Options options) throws ParseException {
        Option option = new Option("r", "ref_genome", true, "reference genome file absolute path");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("s", "source_dir", true, "source exome-seq data directory, data will be used for reads mapping and SNP calling");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("fmt", "input_format", true, "input file format, SRA, FastQ and BAM, default SRA");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("d", "direct", true, "SE or PE, represent single-end and pair-end, default SE");
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

        option = new Option("b", "bwa", true, "BWA alignment tool executive file path, default bwa");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("smt", "samtool", true, "samtools executive file path, default samtools");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("bft", "bcftool", true, "bcftools executive file path, default bcftools");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("picard", "picard_tool", true, "PICARD executive file path, default picard.jar");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("t", "threads", true, "number of working threads, default 2");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("h", "help", false, "user guide");
        option.setRequired(false);
        options.addOption(option);

        CommandLineParser parser = new DefaultParser();

        return parser.parse(options, args);
    }
}
