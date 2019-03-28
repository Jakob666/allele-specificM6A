package GatkSNPCalling;

import org.apache.commons.cli.*;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;

public class GatkSNPCalling {
    private static Logger logger;

    public static void main(String[] args) throws ParseException {
        Options options = new Options();
        CommandLine commandLine = setCommand(args, options);

        // initialize the default params
        String sourceDataDir = "";
        String fastqTempDir = null;
        String genomeFile = "";
        String outputDir = "";
        String gatkLocalJar = "gatk";
        String picardLocalJar = "./picard.jar";
        String inputFormat = "sra";
        int execThread = 2;

        // change params according to the received arguments
        if (commandLine.hasOption('h')) {
            new HelpFormatter().printHelp("flume-ng agent", options, true);
            return;
        }
        if (commandLine.hasOption('r')) {
            genomeFile = commandLine.getOptionValue('r');
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
        if (commandLine.hasOption("gatk")) {
            gatkLocalJar = commandLine.getOptionValue("gatk");
        }
        if (commandLine.hasOption("picard")) {
            picardLocalJar = commandLine.getOptionValue("picard");
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
                System.out.println("transform failed");
            }
        }

        SNPCalling.snpCalling(genomeFile, fastqTempDir, outputDir, picardLocalJar, gatkLocalJar, execThread, logger);

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

        option = new Option("fmt", "input-format", true, "input file format, sra or fastq, default sra");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("tmp", "temp-dir", true, "when input files in sra format, a temporary directory store fastq file generate by sra file, default PWD");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("o", "output", true, "absolute path of the output directory, default PWD");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("gatk", "gatk-tool", true, "GATK executive file path, default gatk");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("picard", "picard-tool", true, "PICARD executive file path, default ./picard.jar");
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
     * initialize a log4j Logger instance for further logging
     * @param logHome output directory of log file
     * @return Logger instance
     */
    private static Logger initLog(String logHome) {
        System.setProperty("log_home", logHome);
        Logger logger = Logger.getLogger(GatkSNPCalling.class);

        return logger;
    }
}
