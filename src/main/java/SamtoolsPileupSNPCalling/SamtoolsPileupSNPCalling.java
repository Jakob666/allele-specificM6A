package SamtoolsPileupSNPCalling;

import GatkSNPCalling.sra2fastq;
import org.apache.commons.cli.*;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;

public class SamtoolsPileupSNPCalling {

    private static Logger logger;

    public static void main(String[] args) throws ParseException {
        Options options = new Options();
        CommandLine commandLine = setCommand(args, options);

        String sourceDataDir = "";
        String fastqTempDir = null;
        String genomeFile = "";
        String gtfFileDir = null;
        String outputDir = "";
        String samtools = "samtools";
        String inputFormat = "sra";
        int execThread = 2;

        if (commandLine.hasOption('h')) {
            new HelpFormatter().printHelp("flume-ng agent", options, true);
            return;
        }
        if (commandLine.hasOption('r')) {
            genomeFile = commandLine.getOptionValue('r');
        }
        if (commandLine.hasOption('g')) {
            gtfFileDir = commandLine.getOptionValue('g');
        }
        if (commandLine.hasOption("fmt")) {
            inputFormat = commandLine.getOptionValue("fmt").toLowerCase();
            if (!inputFormat.equals("sra") && !inputFormat.equals("fastq")) {
                System.out.println("invalid input file format, the format should be sra or fastq");
                System.exit(2);
            }
        }
        if (commandLine.hasOption('s')) {
            sourceDataDir = commandLine.getOptionValue('s');
            if (inputFormat.equals("sra")) {
                if (commandLine.hasOption("tmp")) {
                    fastqTempDir = commandLine.getOptionValue("tmp");
                } else {
                    try {
                        File directory = new File("");
                        fastqTempDir = directory.getCanonicalPath();
                    } catch (IOException ie) {
                        ie.printStackTrace();
                        return;
                    }
                }
            } else {
                fastqTempDir = sourceDataDir;
            }
        }

        if (commandLine.hasOption('o')) {
            outputDir = commandLine.getOptionValue('o');
        } else {
            try {
                File directory = new File("");
                outputDir = directory.getCanonicalPath();
            } catch (IOException ie) {
                ie.printStackTrace();
                return;
            }
        }
        if (commandLine.hasOption("smtool")) {
            samtools = commandLine.getOptionValue("smtool");
        }
        if (commandLine.hasOption('t')) {
            execThread = Integer.parseInt(commandLine.getOptionValue('t'));
        }

        logger = initLog(outputDir);
        mkDir(outputDir);

        // make directories for fastq files and alignment result
        if (inputFormat.equals("sra")) {
            mkDir(fastqTempDir);
            boolean sraTransRes = sraToFastq(sourceDataDir, fastqTempDir);
            if (!sraTransRes) {
                logger.error("transform failed");
            }
        }

        snpCalling(genomeFile, new File(fastqTempDir), outputDir, gtfFileDir, samtools, execThread);

        try {
            Process p = Runtime.getRuntime().exec("rm -rf " + fastqTempDir);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                System.out.println("remove redundant fastq file failed");
                System.exit(2);
            }
        } catch (IOException | InterruptedException ie) {
            ie.printStackTrace();
        }
    }

    /**
     * samtools SNP calling procedure
     * @param genomeFilePath reference genome file path
     * @param fastqDir fastq file directory
     * @param outputDir output directory
     * @param gtfDir gtf files directory
     * @param samtools samtools executive file
     * @param execThread number of working threads
     */
    private static void snpCalling(String genomeFilePath, File fastqDir, String outputDir, String gtfDir, String samtools, int execThread) {
        File[] fastqFiles = fastqDir.listFiles();
        if (fastqFiles == null) {
            System.out.println("empty fastq file Dir");
            System.exit(2);
        }

        for (File fq : fastqFiles) {
            String fileName = fq.getName();
            String prefix = fileName.substring(0, fileName.lastIndexOf("."));
            ReadsMapping.alignment(genomeFilePath, fq.getAbsolutePath(), execThread, logger);
            String refGenomeDir = new File(genomeFilePath).getParent();
            String aligmentResultFile = new File(refGenomeDir, "Aligned.out.sam").getAbsolutePath();
            String dedupBamFile = SamtoolsProcessing.samFileProcess(aligmentResultFile, outputDir, prefix, samtools, logger);
            String readsCountFile = AseInference.inferenceASE(genomeFilePath, dedupBamFile, samtools, logger);
            SnpFilter sf = new SnpFilter(gtfDir, readsCountFile, logger);
            sf.filterVcf();
        }
    }

    /**
     * make up directory if it is not existed
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
     * trans sra format to fastq format for the files under
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
     * set received params for command line
     * @param args arguments in command line
     * @param options Options instance
     * @return CommandLine instance
     * @throws ParseException
     */
    private static CommandLine setCommand(String[] args, Options options) throws ParseException{

        Option option = new Option("r", "ref-genome", true, "reference genome file absolute path");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("s", "source-dir", true, "sequence source data directory");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("g", "gtf-dir", true, "gtf data directory");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("fmt", "input-format", true, "input file format, sra or fastq, default sra");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("tmp", "temp-dir", true, "when input files in sra format, a temporary directory store fastq file generate by sra file, default PWD");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("o", "output", true, "absolute path of the output directory, default PWD");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("smtool", "samtools", true, "samtools executive file path, default samtools");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("t", "threads", true, "number of working threads, default 2");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("h", "help", false, "display help text");
        options.addOption(option);

        CommandLineParser parser = new DefaultParser();
        CommandLine commandLine = parser.parse(options, args);

        return commandLine;
    }

    /**
     * initialize a log4j Logger instance for further logging
     * @param logHome output directory of log file
     * @return Logger instance
     */
    private static Logger initLog(String logHome) {
        System.setProperty("log_home", logHome);
        Logger logger = Logger.getLogger(SamtoolsPileupSNPCalling.class);

        return logger;
    }
}
