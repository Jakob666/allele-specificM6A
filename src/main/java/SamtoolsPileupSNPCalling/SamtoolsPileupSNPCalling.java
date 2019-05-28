package SamtoolsPileupSNPCalling;

import org.apache.commons.cli.*;
import org.apache.log4j.Logger;

public class SamtoolsPileupSNPCalling {

    public static void main(String[] args) throws ParseException{
        Options options = new Options();
        CommandLine commandLine = setCommand(args, options);
        String genomeFile, bamFile, samtools = "samtools", bcftools = "bcftools", outputDir = System.getProperty("user.dir");
        Logger log;

        genomeFile = commandLine.getOptionValue("r");
        bamFile = commandLine.getOptionValue("i");
        if (commandLine.hasOption("smt"))
            samtools = commandLine.getOptionValue("smt");
        if (commandLine.hasOption("bft"))
            bcftools = commandLine.getOptionValue("bft");
        if (commandLine.hasOption("o"))
            outputDir = commandLine.getOptionValue("o");

        log = initLog(outputDir);
        String outputFile = snpCalling(genomeFile, outputDir, bamFile, samtools, bcftools, log);
        log.debug("SNP calling completed. Output file " + outputFile);
    }

    /**
     * samtools SNP calling procedure
     * @param genomeFilePath reference genome file path
     * @param outputDir output directory
     * @param samtools samtools executive file
     */
    public static String snpCalling(String genomeFilePath, String outputDir, String bamFile, String samtools,
                                   String bcftools, Logger logger) {

        BcftoolSNP bcftoolSNP = BcftoolSNP.createCaller(genomeFilePath, bamFile, outputDir, samtools, bcftools, logger);
        bcftoolSNP.mpileUp();
        bcftoolSNP.callSnp();

        return bcftoolSNP.filterSnp();
    }

    private static CommandLine setCommand(String[] args, Options options) throws ParseException {
        Option option = new Option("r", "reference_genome", true, "reference genome file path");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("i", "input_bam", true, "INPUT sample alignment BAM output file path");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("o", "output_dir", true, "result output directory, default present working directory");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("smt", "samtools", true, "samtools executive file path, default samtools");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("bft", "bcftools", true, "bcftools executive file path, default bcftools");
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
        return Logger.getLogger(SamtoolsPileupSNPCalling.class);
    }
}
