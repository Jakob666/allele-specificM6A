package AseM6aPeakDetector;

import AseSeqSimulator.AseSeqSimulator;
import HierarchicalBayesianAnalysis.AseGeneDetection;
import HierarchicalBayesianAnalysis.AsmPeakDetection;
import org.apache.commons.cli.*;

public class RunExecution {
    public static void main(String[] args) {
        Options options = new Options();
        CommandLine commandLine = null;
        HelpFormatter help = new HelpFormatter();
        String header = "renlab.m6a_allele-1.0.jar provides the following tools: ";
        String footer = "";

        try {
            commandLine = setCommandLine(options, args);
        } catch (ParseException pe) {
            System.err.println(pe.getMessage());
            help.printHelp("java -jar renlab.m6a_allele-1.0.jar", header, options, footer, true);
            System.exit(2);
        }

        int runMode = 0;

        if (commandLine.hasOption("AseGeneDetection"))
            runMode = 1;
        else if (commandLine.hasOption("AsmPeakDetection"))
            runMode = 2;
        else if (commandLine.hasOption("AseSeqSimulator"))
            runMode = 3;
        else if (commandLine.hasOption("h")) {
            help.printHelp("java -jar renlab.m6a_allele-1.0.jar", header, options, footer, true);
            System.exit(0);
        } else {
            help.printHelp("java -jar renlab.m6a_allele-1.0.jar", header, options, footer, true);
            System.exit(0);
        }

        execute(runMode, args);
    }

    private static void execute(int runMode, String[] args) {
        if (runMode == 1)
            AseGeneDetection.main(args);
        else if (runMode == 2)
            AsmPeakDetection.main(args);
        else
            AseSeqSimulator.main(args);
    }

    private static CommandLine setCommandLine(Options options, String[] args) throws ParseException {
        Option option = new Option("AseGeneDetection", "detect allele-specific expression genes in test data");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("AsmPeakDetection", "detect allele-specific modification m6A signals in test data");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("AseSeqSimulator", "generate simulation data for test");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("h", "help", false, "show help message and exit program");
        option.setRequired(false);
        options.addOption(option);

        CommandLineParser parser = new DefaultParser();

        return parser.parse(options, args);
    }
}
