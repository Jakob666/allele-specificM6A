package RnaSeqSnpCalling;

import CommonThreadClass.*;
import CommonThreadClass.Sra2Fastq;
import org.apache.commons.cli.*;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

public class RnaSeqSnpCalling {

    private String ipParamVal, inputParamVal, genomeFile, gtfFile, fastqTempDir, outputDir, direct, gatk, picard, samtools, bcftools, inputFormat;
    private int execThread;
    private boolean zip;
    private Logger logger;

    public static void main(String[] args) throws ParseException {
        Options options = new Options();
        CommandLine commandLine = setCommand(args, options);

        // initialize the default params
        String ipDataDir=null, inputDataDir=null, genomeFile = null, gtfFile = null, fastqTempDir = null, outputDir, direct="se";
        String gatkJar = "gatk";
        String picardJar = "picard";
        String samtools = "samtools";
        String bcftools = "bcftools";
        String inputFormat = "sra";
        int execThread = 2;
        boolean zipFormat = false;

        if (commandLine.hasOption('h')) {
            new HelpFormatter().printHelp("flume-ng agent", options, true);
            return;
        }

        // get reference genome file and GTF file
        if (commandLine.hasOption("r"))
            genomeFile = commandLine.getOptionValue('r');
        if (commandLine.hasOption("g"))
            gtfFile = commandLine.getOptionValue('g');
        if (genomeFile == null) {
            System.out.println("require reference genome file");
            System.exit(2);
        } else if (new File(genomeFile).exists() && !new File(genomeFile).isFile()) {
            System.out.println(genomeFile + " should be a file, not a directory.");
            System.exit(2);
        } else if (!new File(genomeFile).exists()){
            System.out.println("file " + genomeFile + " doesn't exist");
            System.exit(2);
        }
        if (gtfFile == null) {
            System.out.println("require reference genome file");
            System.exit(2);
        } else if (new File(gtfFile).exists() && !new File(gtfFile).isFile()) {
            System.out.println(gtfFile + " should be a file, not a directory.");
            System.exit(2);
        } else if (!new File(gtfFile).exists()){
            System.out.println("file " + gtfFile + " doesn't exist");
            System.exit(2);
        }

        String genomeFileDir = new File(genomeFile).getParent();
        // output result directory, default a new directory name "outputResult" in genome file directory
        if (commandLine.hasOption('o')) {
            outputDir = commandLine.getOptionValue('o');
            mkDir(outputDir);
        }
        else {
            outputDir = new File(genomeFileDir, "outputResult").getAbsolutePath();
            mkDir(outputDir);
        }
        Logger logger = initLog(outputDir);
        if (commandLine.hasOption("d")) {
            direct = commandLine.getOptionValue("d");
            if (!direct.toLowerCase().equals("se") && !direct.toLowerCase().equals("pe")) {
                logger.error("invalid input sequencing direction, must be SE(single-end) or PE(pair-end)");
                System.exit(2);
            }
        }

        // get input file format(support sra, fastq, fasta)
        if (commandLine.hasOption("fmt")) {
            inputFormat = commandLine.getOptionValue("fmt").toLowerCase();
            if (!inputFormat.equals("sra") && !inputFormat.equals("fastq") && !inputFormat.equals("bam")) {
                logger.error("invalid input file format, the format should be sra, fastq or bam");
                System.exit(2);
            }
        }
        if (commandLine.hasOption("zip")) {
            String val = commandLine.getOptionValue("zip");
            if (!val.equals("true") && !val.equals("false")) {
                logger.error("invalid value. -zip parameter should be true or false");
                System.exit(2);
            }
            zipFormat = Boolean.parseBoolean(val);
        }
        if (inputFormat.equals("sra"))
            zipFormat = true;
        // source data directory
        if (commandLine.hasOption("ip"))
            ipDataDir = commandLine.getOptionValue("ip");
        if (commandLine.hasOption("input"))
            inputDataDir = commandLine.getOptionValue("input");
        // 如果没有给出IP或者INPUT样本则报错
        if (ipDataDir == null || inputDataDir == null) {
            logger.error("Specifies the value of at least one of the parameters -ip(--ip_dir) and -input(--input_dir)");
            System.exit(2);
        }

        if (inputFormat.equals("sra")) {
            if (commandLine.hasOption("tmp"))
                fastqTempDir = commandLine.getOptionValue("tmp");
            else {
                fastqTempDir = new File(genomeFileDir, "temp").getAbsolutePath();
                mkDir(fastqTempDir);
            }
        } else {
            fastqTempDir = ipDataDir;
        }
        if (commandLine.hasOption("smt"))
            samtools = commandLine.getOptionValue("smt");
        if (commandLine.hasOption("bft"))
            bcftools = commandLine.getOptionValue("bft");
        if (commandLine.hasOption("gatk"))
            gatkJar = commandLine.getOptionValue("gatk");
        if (commandLine.hasOption("picard"))
            picardJar = commandLine.getOptionValue("picard");
        if (commandLine.hasOption('t'))
            execThread = Integer.parseInt(commandLine.getOptionValue('t'));

        RnaSeqSnpCalling rssc = new RnaSeqSnpCalling(ipDataDir, inputDataDir, genomeFile, gtfFile, fastqTempDir, outputDir,
                direct, gatkJar, picardJar, samtools, bcftools, inputFormat, execThread, zipFormat, logger);
        rssc.callSnp();
    }

    public RnaSeqSnpCalling(String ipParamVal, String inputParamVal, String genomeFile, String gtfFile, String fastqTempDir,
                            String outputDir, String direct, String gatk, String picard, String samtools, String bcftools,
                            String inputFormat, int execThread, boolean zip, Logger logger) {
        this.ipParamVal = ipParamVal;
        this.inputParamVal = inputParamVal;
        this.genomeFile = genomeFile;
        this.gtfFile = gtfFile;
        this.fastqTempDir = fastqTempDir;
        this.outputDir = outputDir;
        this.direct = direct;
        this.gatk = gatk;
        this.picard = picard;
        this.samtools = samtools;
        this.bcftools = bcftools;
        this.inputFormat = inputFormat;
        this.execThread = execThread;
        this.zip = zip;
        this.logger = logger;
    }

    public void callSnp() {
        // 依据用户传入的 -ip和-input参数值得到此次分析的文件
        ArrayList<File> ipSamples = getUserFiles(this.direct, this.ipParamVal, this.inputFormat);
        ArrayList<File> inputSamples = getUserFiles(this.direct, this.inputParamVal, this.inputFormat);

        // 开启INPUT样本的分析
        String inputOutputDir = new File(outputDir, "INPUT").getAbsolutePath();
        mkDir(inputOutputDir);
        File[] fastqFiles;
        ArrayList<String> fastqFileNames;
        RnaSeqReadsMapping rsrm;
        // 对于非BAM格式文件需要先reads Mapping
        if (!inputFormat.toLowerCase().equals("bam")) {
            if (inputFormat.toLowerCase().equals("sra")) {
                fastqFileNames = transToFormat(inputSamples, fastqTempDir, direct);
            } else {
                fastqFileNames = new ArrayList<>();
                for (File f: inputSamples)
                    fastqFileNames.add(f.getAbsolutePath());
            }
            System.out.println(fastqFileNames);
            fastqFiles = new File[fastqFileNames.size()];
            for (int i=0; i< fastqFileNames.size(); i++) {
                File f = new File(fastqFileNames.get(i));
                fastqFiles[i] = f;
            }

            rsrm = new RnaSeqReadsMapping(fastqFiles, genomeFile, gtfFile, picard, gatk, samtools,
                    inputOutputDir, "ase", direct, execThread, zip, logger);
            Thread inputThread = new Thread(rsrm);
            inputThread.start();
            try {
                inputThread.join();
            } catch (InterruptedException ie) {
                logger.error(ie.getMessage());
                System.exit(2);
            }
            if (inputFormat.toLowerCase().equals("sra")) {
                cleanUp(fastqTempDir);
                mkDir(fastqTempDir);
            }
        } else { // BAM格式文件如果有多个，则需要进行合并
            File mergedBamFile = new File(inputOutputDir, "ase_alignment.bam");
            this.mergeBamFiles(inputSamples, mergedBamFile.getAbsolutePath());
        }

        // 开启IP样本的分析
        Thread ipThread = null;
        String ipOutputDir = new File(outputDir, "IP").getAbsolutePath();
        mkDir(ipOutputDir);
        if (!inputFormat.toLowerCase().equals("bam")) {
            if (inputFormat.toLowerCase().equals("sra")) {
                fastqFileNames = transToFormat(ipSamples, fastqTempDir, direct);
            } else {
                fastqFileNames = new ArrayList<>();
                for (File f: ipSamples)
                    fastqFileNames.add(f.getAbsolutePath());
            }
            fastqFiles = new File[fastqFileNames.size()];
            for (int i=0; i< fastqFileNames.size(); i++) {
                File f = new File(fastqFileNames.get(i));
                fastqFiles[i] = f;
            }

            rsrm = new RnaSeqReadsMapping(fastqFiles, genomeFile, gtfFile, picard, gatk, samtools, ipOutputDir,
                    "asm", direct, execThread, zip, logger);
            ipThread = new Thread(rsrm);

        } else {
            File mergedBamFile = new File(ipOutputDir, "asm_alignment.bam");
            this.mergeBamFiles(ipSamples, mergedBamFile.getAbsolutePath());
        }

        // 创建线程用于INPUT样本SNP calling
        String inputBamFile = new File(inputOutputDir, "ase_alignment.bam").getAbsolutePath();
        RunSnpCalling inputRunSnp = new RunSnpCalling(genomeFile, inputBamFile, samtools, bcftools, logger);
        Thread inputSnpThread = new Thread(inputRunSnp);

        // 创建线程用于IP样本SNP calling
        String ipBamFile = new File(ipOutputDir, "asm_alignment.bam").getAbsolutePath();
        RunSnpCalling ipRunSnp = new RunSnpCalling(genomeFile, ipBamFile, samtools, bcftools, logger);
        Thread ipSnpThread = new Thread(ipRunSnp);

        if (ipThread != null)
            ipThread.start();
        inputSnpThread.start();
        // 主线程等待子线程运行完才结束
        try {
            if (ipThread != null)
                ipThread.join();

            if (ipThread != null) {
                while (ipThread.isAlive())
                    Thread.sleep(60000);
            }
            ipSnpThread.start();
            inputSnpThread.join();
            ipSnpThread.join();
        } catch (InterruptedException ie) {
            logger.error(ie.getMessage());
            System.exit(2);
        }

        if (inputFormat.toLowerCase().equals("sra"))
            cleanUp(fastqTempDir);
    }

    /**
     * 依据用户指定的 -ip或-input参数值获取数据
     * @param direct 测序方向，se或pe
     * @param dir -ip或-input参数值
     * @param fileFormat 输入文件的格式 sra、fastq或bam
     * @return ArrayList<File>包含需要分析的文件
     */
    public ArrayList<File> getUserFiles(String direct, String dir, String fileFormat) {
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
                                if (!this.zip && (f.getAbsolutePath().endsWith("fastq") || f.getAbsolutePath().endsWith("fq"))) {
                                    userFiles.add(f);
                                } else if (this.zip && (f.getAbsolutePath().endsWith("fastq.gz") || f.getAbsolutePath().endsWith("fq.gz"))) {
                                    userFiles.add(f);
                                }
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
                            if (!this.zip && files.getAbsolutePath().endsWith("q"))
                                userFiles.add(files);
                            else if (this.zip && files.getAbsolutePath().endsWith("fastq.gz"))
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
                    if (fileFormat.toLowerCase().equals("fastq") && !this.zip) { // 如果是未压缩的FastQ格式
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
                    } else if (fileFormat.toLowerCase().equals("fastq") && this.zip) { // 如果是压缩的FastQ格式
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
                    }else {  // 如果是SRA格式
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
                        if (!this.zip) {
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
     * 将SRA文件转为FastQ格式
     * @param sraFiles sra文件列表
     * @param fastqTempDir 临时存放目录
     * @param direct 测序方向
     * @return FastQ文件列表
     */
    public ArrayList<String> transToFormat(ArrayList<File> sraFiles, String fastqTempDir, String direct) {
        Sra2Fastq sf = new Sra2Fastq("", fastqTempDir, direct, sraFiles, logger);

        return sf.transformat();
    }

    /**
     * 合并用户传入的BAM文件
     * @param bamFiles BAM文件列表
     * @param mergedBamFilePath 输出文件路径
     */
    public void mergeBamFiles (ArrayList<File> bamFiles, String mergedBamFilePath) {
        String cmd;
        if (this.inputFormat.toLowerCase().equals("bam")) {
            cmd = "cp " + bamFiles.get(0).getAbsolutePath() + " " + mergedBamFilePath;
            try {
                Process p = Runtime.getRuntime().exec(cmd);
                int res = p.waitFor();
                if (res != 0) {
                    logger.error("rename file error: " + bamFiles.get(0));
                    System.exit(2);
                }
            } catch (IOException | InterruptedException ie) {
                logger.error(ie.getMessage());
                System.exit(2);
            }
        } else {
            String[] filepaths = new String[bamFiles.size()];
            for (int i = 0; i < bamFiles.size(); i++) {
                filepaths[i] = bamFiles.get(i).getAbsolutePath();
            }
            cmd = String.join(" ", new String[]{samtools, "merge", mergedBamFilePath, String.join(" ", filepaths)});
            try {
                Process p = Runtime.getRuntime().exec(cmd);
                int res = p.waitFor();
                if (res != 0) {
                    logger.error("merge bam files failed");
                    System.exit(2);
                }
            } catch (IOException | InterruptedException ie) {
                logger.error(ie.getMessage());
                System.exit(2);
            }
        }
    }

    /**
     * 如果目录不存在则创建
     * @param dirName directory name
     */
    public static void mkDir(String dirName) {
        File targetDir = new File(dirName);
        if (!targetDir.exists()) {
            boolean res = targetDir.mkdirs();
        }
        targetDir = null;
    }

    /**
     * 清理临时目录中的fastq文件
     * @param dirName temporary directory name
     */
    public void cleanUp(String dirName) {
        if (dirName == null)
            return;
        String cmd = "rm -rf " + dirName;
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitNum = p.waitFor();
            if (exitNum != 0) {
                this.logger.error("remove redundant fastq file failed");
            }
        } catch (IOException | InterruptedException ie) {
            this.logger.error(ie.getMessage());
        }
    }

    /**
     * 设置命令行参数
     * @param args arguments in command line
     * @param options Options instance
     * @return CommandLine instance
     * @throws ParseException throws Exception
     */
    private static CommandLine setCommand(String[] args, Options options) throws ParseException{

        Option option = new Option("r", "ref_genome", true, "reference genome file absolute path");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("g", "gtf_file", true, "gene transfer format(GTF) file");
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

        option = new Option("fmt", "input_format", true, "input file format, SRA, FastQ or BAM, default SRA");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("zip", "fastq_zip", true, "fastQ data in GZIP compressed format. default false");
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

