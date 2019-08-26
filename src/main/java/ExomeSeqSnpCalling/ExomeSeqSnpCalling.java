package ExomeSeqSnpCalling;

import CommonThreadClass.RunSnpCalling;
import CommonThreadClass.Sra2Fastq;
import org.apache.commons.cli.*;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

public class ExomeSeqSnpCalling {
    private String sParamVal, genomeFile, fastqTempDir, outputDir, prefix, direct, bwa, sraTool, picard,
                   samtools, bcftools, inputFormat;
    private int execThread;
    private boolean compress;
    private static Logger logger;

    public static void main(String[] args) throws ParseException {
        Options options = new Options();
        CommandLine commandLine = setCommand(args, options);
        String helpMessage = "";
        if (commandLine.hasOption("h")) {
            System.out.println(helpMessage);
            System.exit(0);
        }
        String source_dir = null, genomeFile = null, fastqTempDir, outputDir, prefix, direct = "SE";
        String bwa = "bwa", sraTool = "", picardJar = "picard", samtools = "samtools", bcftools = "bcftools", inputFormat = "sra";
        int execThread = 2;
        boolean compressed = false;

        if (commandLine.hasOption("r"))
            genomeFile = commandLine.getOptionValue('r');
        if (commandLine.hasOption("s"))
            source_dir = commandLine.getOptionValue("s");
        if (genomeFile == null || source_dir == null) {
            System.out.println("reference genome file and source data dir is required");
            System.exit(2);
        }
        if (commandLine.hasOption('o')) {
            outputDir = commandLine.getOptionValue('o');
            mkDir(outputDir);
        } else {
            outputDir = new File(System.getProperty("user.dir"), "outputResult").getAbsolutePath();
            mkDir(outputDir);
        }
        logger = initLog(outputDir);


        if (commandLine.hasOption('p'))
            prefix = commandLine.getOptionValue('p');
        else
            prefix = "exome";
        if (commandLine.hasOption("b"))
            bwa = commandLine.getOptionValue("b");
        if (commandLine.hasOption("smt"))
            samtools = commandLine.getOptionValue("smt");
        if (commandLine.hasOption("bft"))
            bcftools = commandLine.getOptionValue("bft");
        if (commandLine.hasOption("picard"))
            picardJar = commandLine.getOptionValue("picard");
        if (commandLine.hasOption('t'))
            execThread = Integer.parseInt(commandLine.getOptionValue('t'));
        if (commandLine.hasOption("fmt")) {
            inputFormat = commandLine.getOptionValue("fmt").toLowerCase();
            if (!inputFormat.equals("sra") && !inputFormat.equals("fastq") && !inputFormat.equals("bam")) {
                logger.error("invalid input file format, the format should be sra, fastq or bam");
                System.exit(2);
            }
        }
        if (inputFormat.toLowerCase().equals("fastq") && commandLine.hasOption("c")) {
            String val = commandLine.getOptionValue("c");
            if (!val.toLowerCase().equals("true") && !val.toLowerCase().equals("false")) {
                logger.error("invalid input value of c parameter");
                System.exit(2);
            }
            compressed = Boolean.parseBoolean(val);
        }
        if (commandLine.hasOption("d")) {
            direct = commandLine.getOptionValue("d");
            if (direct.toLowerCase().equals("se"))
                direct = "se";
            else if (direct.toLowerCase().equals("pe"))
                direct = "pe";
            else {
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

        ExomeSeqSnpCalling essc = new ExomeSeqSnpCalling(source_dir, genomeFile, fastqTempDir, outputDir, prefix, direct,
                                                         bwa, sraTool, picardJar, samtools, bcftools, inputFormat,
                                                         execThread, compressed);
        essc.callSnp();
    }

    public ExomeSeqSnpCalling(String sParamVal, String genomeFile, String fastqTempDir, String outputDir, String prefix,
                              String direct, String bwa, String sraTool, String picard,String samtools, String bcftools,
                              String inputFormat, int execThread, boolean compressed) {
        this.sParamVal = sParamVal;
        this.genomeFile = genomeFile;
        this.fastqTempDir = fastqTempDir;
        this.outputDir = outputDir;
        this.prefix = prefix;
        this.direct = direct;
        this.bwa = bwa;
        this.sraTool = sraTool;
        this.picard = picard;
        this.samtools = samtools;
        this.bcftools = bcftools;
        this.inputFormat = inputFormat;
        this.execThread = execThread;
        this.compress = compressed;
    }

    public void callSnp() {
        File alignmentOutputDir = new File(outputDir, "alignmentResult");
        mkDir(alignmentOutputDir.getAbsolutePath());
        // 依据用户传入的 -s参数值得到分析的文件
        ArrayList<File> userFiles = getUserFiles(direct, sParamVal, inputFormat);

        String mergedBamFile;
        // 如果上传的数据不是BAM格式，则需先对其进行比对
        if (!inputFormat.toLowerCase().equals("bam")) {
            // 对测序文件进行比对
            logger.debug("start exome sequencing reads alignment");
            mergedBamFile = this.alignSequencingReads(userFiles, alignmentOutputDir.getAbsolutePath());
            if (inputFormat.toLowerCase().equals("sra"))
                cleanUp(fastqTempDir);
            logger.debug("Exome sequencing reads alignment complete");
        } else {
            File mergedBam = new File(outputDir, prefix+"_alignment.bam");
            mergedBamFile = mergedBam.getAbsolutePath();
            ArrayList<String> alignmentBamFile = new ArrayList<>();
            for (File f: userFiles) {
                if (f.getAbsolutePath().endsWith("bam")) {
                    alignmentBamFile.add(f.getAbsolutePath());
                }
            }
            if (alignmentBamFile.size() == 0) {
                logger.error("directory contains no bam File");
                System.exit(2);
            }
            logger.debug("merge alignment bam files");
            if (!this.inputFormat.toLowerCase().equals("bam") && alignmentBamFile.size() == 1) {
                boolean res = new File(alignmentBamFile.get(0)).renameTo(mergedBam);
                if (!res) {
                    logger.error("rename file error: " + alignmentBamFile.get(0));
                }
            } else if (this.inputFormat.toLowerCase().equals("bam") && alignmentBamFile.size() == 1) {
                String cmd = "cp " + alignmentBamFile.get(0) + " " + mergedBamFile;
                try {
                    Process p = Runtime.getRuntime().exec(cmd);
                    int res = p.waitFor();
                    if (res != 0) {
                        logger.error("cp BAM file to output directory failed");
                        System.exit(2);
                    }
                } catch (IOException | InterruptedException ie) {
                    logger.error(ie.getMessage());
                    System.exit(2);
                }
            } else {
                ReadsAlignment readsAlignment = new ReadsAlignment(bwa, samtools, picard, genomeFile, execThread, logger);
                readsAlignment.mergeBamFiles(alignmentBamFile, mergedBam.getAbsolutePath());
            }
        }

        // 对样本进行SNP calling
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
    private String alignSequencingReads(ArrayList<File> userFiles, String outputDir) {
        ReadsAlignment readsAlignment = new ReadsAlignment(bwa, samtools, picard, genomeFile, execThread, logger);

        // 先对参考基因建立索引。之后对每个FastQ文件进行比对，比对后对重复样本进行合并
        logger.debug("build index for reference genome: " + genomeFile);
        readsAlignment.genomeIndex();
        ArrayList<String> alignmentBamFiles = new ArrayList<>();
        ArrayList<String> fastqFileNames;
        if (direct.toLowerCase().equals("se")) {    // single-end测序
            if (inputFormat.toLowerCase().equals("sra")) {
                Sra2Fastq sra2Fastq = new Sra2Fastq(sraTool, fastqTempDir, direct, userFiles, logger);
                fastqFileNames = sra2Fastq.transformat();
            } else {
                fastqFileNames = new ArrayList<>();
                for (File f: userFiles)
                    fastqFileNames.add(f.getAbsolutePath());
            }
            File[] fastqFiles = new File[fastqFileNames.size()];
            for (int i = 0; i < fastqFileNames.size(); i++)
                fastqFiles[i] = new File(fastqFileNames.get(i));

            for (File inputFile: fastqFiles) {
                String inputFileName = inputFile.getName();
                String alignmentSamFile = new File(outputDir, inputFileName.substring(0, inputFileName.lastIndexOf("."))+".sam").getAbsolutePath();
                String bamFile = readsAlignment.alignmentToGenome(inputFile.getAbsolutePath(), null, alignmentSamFile);
                alignmentBamFiles.add(bamFile);
            }
        } else {    // pair-end测序
            if (inputFormat.toLowerCase().equals("sra")) {
                Sra2Fastq sra2Fastq = new Sra2Fastq(sraTool, fastqTempDir, direct, userFiles, logger);
                fastqFileNames = sra2Fastq.transformat();
            } else {
                fastqFileNames = new ArrayList<>();
                for (File f: userFiles)
                    fastqFileNames.add(f.getAbsolutePath());
            }
            File[] fastqFiles = new File[fastqFileNames.size()];
            for (int i = 0; i < fastqFileNames.size(); i++)
                fastqFiles[i] = new File(fastqFileNames.get(i));

            for (int i = 0; i < fastqFiles.length; i=i+2) {
                String inputFile1 = fastqFiles[i].getName();
                String inputFile2 = fastqFiles[i+1].getName();
                String alignmentSamFile;
                // 确定SAM文件名
                if (inputFile1.contains("_")) {
                    alignmentSamFile = new File(outputDir,inputFile1.substring(0, inputFile1.lastIndexOf("_")) + ".sam").getAbsolutePath();
                } else {
                    alignmentSamFile = new File(outputDir,inputFile1.substring(0, inputFile1.lastIndexOf(".")) + ".sam").getAbsolutePath();
                }
                String bamFile = readsAlignment.alignmentToGenome(fastqFiles[i].getAbsolutePath(), fastqFiles[i+1].getAbsolutePath(), alignmentSamFile);
                alignmentBamFiles.add(bamFile);
            }
        }
        logger.debug("merge alignment bam files");
        File mergedBam = new File(outputDir, prefix+"_alignment.bam");
        if (alignmentBamFiles.size() == 1) {
            boolean res = new File(alignmentBamFiles.get(0)).renameTo(mergedBam);
            if (!res) {
                logger.error("rename file error: " + alignmentBamFiles.get(0));
            }
        } else {
            readsAlignment.mergeBamFiles(alignmentBamFiles, mergedBam.getAbsolutePath());
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
            boolean res = targetDir.mkdirs();
            if (!res) {
                logger.error("can not make directory: " + dirName);
                System.exit(2);
            }
        }
    }

    /**
     * 依据用户指定的 -s参数值获取数据
     * @param direct 测序方向，se或pe
     * @param dir -s参数值
     * @param fileFormat 输入文件的格式 sra、fastq或bam
     * @return ArrayList<File>包含需要分析的文件
     */
    private ArrayList<File> getUserFiles(String direct, String dir, String fileFormat) {
        File files = new File(dir);
        ArrayList<File> userFiles = new ArrayList<>();
        // 如果传入的不是BAM格式
        if (!fileFormat.toLowerCase().equals("bam")) {
            // 对于单端测序而言
            if (direct.toLowerCase().equals("se")) {
                if (files.exists()) {
                    if (files.isDirectory()) { // 如果文件对象存在且是一个目录则遍历其中文件，如果是空目录则报错
                        File[] samples = files.listFiles();
                        if (samples == null) {
                            logger.error("empty IP data directory " + files);
                            System.exit(2);
                        }
                        if (fileFormat.toLowerCase().equals("fastq")) {
                            for (File f: samples) {
                                if (!this.compress && (f.getAbsolutePath().endsWith("fastq") || f.getAbsolutePath().endsWith("fq")))
                                    userFiles.add(f);
                                else if (this.compress && (f.getAbsolutePath().endsWith("fastq.gz") || f.getAbsolutePath().endsWith("fq.gz")))
                                    userFiles.add(f);
                            }
                            if (userFiles.size() == 0) {
                                logger.error("FastQ files not exist");
                                System.exit(2);
                            }
                        } else {
                            for (File f: samples) {
                                if (!f.getAbsolutePath().endsWith("sra"))
                                    continue;
                                userFiles.add(f);
                            }
                            if (userFiles.size() == 0) {
                                logger.error("SRA files not exist");
                                System.exit(2);
                            }
                        }
                    }
                    else if (files.isFile()) { // 如果文件对象存在且是一个文件，则将其作为要分析的数据
                        if (fileFormat.toLowerCase().equals("fastq")) {
                            if (! this.compress && (files.getAbsolutePath().endsWith("fq") || files.getAbsolutePath().endsWith("fastq")))
                                userFiles.add(files);
                            else if (this.compress && (files.getAbsolutePath().endsWith("fastq.gz") || files.getAbsolutePath().endsWith("fq.gz")))
                                userFiles.add(files);
                        } else if (fileFormat.toLowerCase().equals("sra")) {
                            if (files.getAbsolutePath().endsWith("sra"))
                                userFiles.add(files);
                        }
                        if (userFiles.size() == 0) {
                            logger.error("the format of file " + files + " is incompatible with file format parameter " + fileFormat);
                            System.exit(2);
                        }
                    }
                } // 如果文件对象不存在则用户传参形式为 file1,file2,...
                else {
                    String[] samples = dir.split(",");
                    if (fileFormat.toLowerCase().equals("fastq") && !this.compress) { // 如果是未压缩的FastQ格式
                        for (String sample: samples) {
                            if (!sample.toLowerCase().endsWith("fastq") && !sample.toLowerCase().endsWith("fq"))
                                continue;
                            File f = new File(sample);
                            if (f.exists())
                                userFiles.add(f);
                        }
                        if (userFiles.size() == 0) {
                            logger.error("FastQ files not exist " + String.join(", ", samples));
                            System.exit(2);
                        }
                    } else if (fileFormat.toLowerCase().equals("fastq") && this.compress) { // 如果是压缩的FastQ格式
                        for (String sample: samples) {
                            if (!sample.toLowerCase().endsWith("fastq.gz") && !sample.toLowerCase().endsWith("fq.gz"))
                                continue;
                            File f = new File(sample);
                            if (f.exists())
                                userFiles.add(f);
                        }
                        if (userFiles.size() == 0) {
                            logger.error("FastQ files not exist " + String.join(", ", samples));
                            System.exit(2);
                        }
                    } else {  // 如果是SRA格式
                        for (String sample: samples) {
                            if (!sample.toLowerCase().endsWith("sra"))
                                continue;
                            File f = new File(sample);
                            if (f.exists())
                                userFiles.add(f);
                        }
                        if (userFiles.size() == 0) {
                            logger.error("SRA files not exist " + String.join(", ", samples));
                            System.exit(2);
                        }
                    }
                }

            }  // 如果是双端测序
            else {
                String[] mates = dir.split(";");
                String[] mate1 = new String[mates.length], mate2 = new String[mates.length];
                for (int i = 0; i < mates.length; i++) {
                    String[] pair = mates[i].split(",");
                    if (pair.length != 2) {
                        logger.error("invalid input " + mates[i] + ", pair-end sequencing file need in pair");
                        System.exit(2);
                    }
                    mate1[i] = pair[0];
                    mate2[i] = pair[1];
                }
                if (fileFormat.toLowerCase().equals("fastq")) {
                    for (int i=0; i<mate1.length; i++) {
                        File mate1File = new File(mate1[i]);
                        File mate2File = new File(mate2[i]);
                        if (!this.compress) {
                            if (mate1File.getName().toLowerCase().endsWith("q") && mate1File.exists())
                                userFiles.add(mate1File);
                            if (mate2File.getName().toLowerCase().endsWith("q") && mate1File.exists())
                                userFiles.add(mate2File);
                        } else {
                            if (mate1File.getName().toLowerCase().endsWith("gz") && mate1File.exists())
                                userFiles.add(mate1File);
                            if (mate2File.getName().toLowerCase().endsWith("gz") && mate1File.exists())
                                userFiles.add(mate2File);
                        }
                    }
                    if (!(userFiles.size() % 2 == 0)) {
                        logger.error("files not in pairs");
                        System.exit(2);
                    }
                } else {
                    for (int i=0; i<mate1.length; i++) {
                        File mate1File = new File(mate1[i]);
                        File mate2File = new File(mate2[i]);
                        if (mate1File.getName().toLowerCase().endsWith("sra") && mate1File.exists())
                            userFiles.add(mate1File);
                        if (mate2File.getName().toLowerCase().endsWith("sra") && mate1File.exists())
                            userFiles.add(mate2File);
                    }
                    if (!(userFiles.size() % 2 == 0)) {
                        logger.error("files not in pairs");
                        System.exit(2);
                    }
                }
            }
        } else { // 如果是数据格式为BAM格式文件(无SE和PE)
            // 若传入的参数值是一个文件或目录
            if (files.exists()) {
                if (files.isFile()) {
                    if (files.getAbsolutePath().endsWith("bam")) {
                        userFiles.add(files);
                        return userFiles;
                    }
                    logger.error("file " + files + " is not a BAM format file");
                    System.exit(2);
                } else if (files.isDirectory()) {
                    File[] bamFiles = files.listFiles();
                    if (bamFiles == null) {
                        logger.error("empty BAM file directory " + files);
                        System.exit(2);
                    }
                    for (File f: bamFiles) {
                        if (f.getAbsolutePath().endsWith("bam"))
                            userFiles.add(f);
                    }
                }
            } else { // 若传入的不是文件或目录，则是使用逗号分隔的各个BAM文件
                String[] bamFiles = dir.split(",");
                for (String bam: bamFiles) {
                    if (!bam.toLowerCase().endsWith("bam"))
                        continue;
                    File f = new File(bam);
                    if (f.exists())
                        userFiles.add(f);
                }
                if (userFiles.size() == 0) {
                    logger.error("files not existed " + String.join(", ", bamFiles));
                    System.exit(2);
                }
            }
        }

        return userFiles;
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

        option = new Option("s", "source_data", true, "When single-end sequencing, parameter value in 2 format: \n\t1) source exome-seq data directory: /path/to/source_data_dir" +
                            "\n\t2) sequencing data separate by comma: file1,file2,...\nWhen pair-end sequencing, parameter value like: mate1_1,mate1_2;mate2_1,mate2_2;... set parameter -d or --direct to PE");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("fmt", "input_format", true, "input file format, SRA, FastQ and BAM, default SRA");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("c", "compress", true, "Whether the input FastQ files are compressed, default false");
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
