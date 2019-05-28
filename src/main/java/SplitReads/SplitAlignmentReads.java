package SplitReads;

import ReadsMapping.ReadsMapping;
import org.apache.commons.cli.*;
import org.apache.log4j.Logger;

import java.io.File;

public class SplitAlignmentReads {
    private String genomeFile, inputBamFile, outputBamFile, gatk;
    private Logger logger;

    public static void main(String[] args) throws ParseException {
        Options options = new Options();
        CommandLine commandLine = setCommand(args, options);
        String genomeFile, inputBamFile, outputBamFile, gatk = "gatk";

        genomeFile = commandLine.getOptionValue("r");
        inputBamFile = commandLine.getOptionValue("i");
        if (commandLine.hasOption("o"))
            outputBamFile = commandLine.getOptionValue("o");
        else
            outputBamFile = new File(new File(inputBamFile).getParent(), "split.bam").getAbsolutePath();
        if (commandLine.hasOption("g"))
            gatk = commandLine.getOptionValue("g");

        SplitAlignmentReads sar = new SplitAlignmentReads(genomeFile, inputBamFile, outputBamFile, gatk);
        sar.splitAlignmentReads();

    }

    public SplitAlignmentReads(String genomeFile, String inputFile, String outputFile, String gatk) {
        this.genomeFile = genomeFile;
        this.inputBamFile = inputFile;
        this.outputBamFile = outputFile;
        this.gatk = gatk;
        String outputDir = new File(outputBamFile).getParent();
        this.logger = initLog(outputDir);
        this.splitAlignmentReads();
    }

    private void splitAlignmentReads() {
        ReadsMapping.readsTrimReassign(gatk, genomeFile, inputBamFile, outputBamFile, logger);
        logger.debug("Complete split alignment reads. Output file " + this.outputBamFile);
    }

    private static CommandLine setCommand(String[] args, Options options) throws ParseException{
        Option option = new Option("r", "reference_genome", true, "reference genome file path");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("i", "input_bam", true, "input bam file path");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("o", "output_bam", true, "output bam file path, default readSplit.bam");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("g", "gatk", true, "GATK executive file path, default gatk");
        option.setRequired(false);
        options.addOption(option);

        CommandLineParser parser = new DefaultParser();

        return parser.parse(options, args);
    }

    private static Logger initLog(String logHome) {
        System.setProperty("log_home", logHome);
        return Logger.getLogger(SplitAlignmentReads.class);
    }
}
