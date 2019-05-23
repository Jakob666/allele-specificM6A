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
        String sourceDataDir, genomeFile, gtfFile, refVcfFile, fastqTempDir, outputDir, prefix;
        String gatkJar = "gatk";
        String picardJar = "picard";
        String samtools = "samtools";
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
        refVcfFile = commandLine.getOptionValue("rv");

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
            if (!inputFormat.equals("sra") && !inputFormat.equals("fastq")) {
                logger.error("invalid input file format, the format should be sra or fastq");
                System.exit(2);
            }
        }

        // source data directory
        sourceDataDir = commandLine.getOptionValue('s');
        if (inputFormat.equals("sra")) {
            if (commandLine.hasOption("tmp")) {
                fastqTempDir = commandLine.getOptionValue("tmp");
            } else {
                fastqTempDir = new File(genomeFileDir, "temp").getAbsolutePath();
                mkDir(fastqTempDir);
            }
        } else {
            fastqTempDir = sourceDataDir;
        }

        if (commandLine.hasOption('p')) {
            prefix = commandLine.getOptionValue('p');
        } else {
            prefix = "AseM6a";
        }
        if (commandLine.hasOption("smt")) {
            samtools = commandLine.getOptionValue("smt");
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

        // make directories for fastq files and alignment result
        if (inputFormat.equals("sra")) {
            mkDir(fastqTempDir);
            boolean sraTransRes = sraToFastq(sourceDataDir, fastqTempDir);
            if (!sraTransRes) {
                logger.error("transform failed");
                System.exit(2);
            }
        }

        SNPCalling.snpCalling(genomeFile, gtfFile, refVcfFile, prefix, fastqTempDir, outputDir, picardJar, gatkJar, samtools, execThread, logger);

        if (!fastqTempDir.equals(sourceDataDir)) {
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

        option = new Option("g", "gtf-file", true, "gene transfer format(GTF) file");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("rv", "ref_vcf", true, "reference VCF file, recommend dbSNP");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("s", "source-dir", true, "sequence source data directory");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("fmt", "input-format", true, "input file format, sra or fastq, default sra");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("p", "prefix", true, "output result file prefix, default 'AseM6a'");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("tmp", "temp-dir", true, "when input files in sra format, a temporary directory store fastq file generate by sra file, default PWD");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("o", "output", true, "absolute path of the output directory, default PWD");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("smt", "smtool", true, "samtools executive file path, default samtools");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("gatk", "gatk-tool", true, "GATK executive file path, default gatk");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("picard", "picard-tool", true, "PICARD executive file path, default picard.jar");
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
