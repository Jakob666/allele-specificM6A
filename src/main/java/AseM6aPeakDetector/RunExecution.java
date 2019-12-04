package AseM6aPeakDetector;

import AseSeqSimulator.AseSeqSimulator;
import HierarchicalBayesianAnalysis.AseGeneDetection;
import HierarchicalBayesianAnalysis.AsmPeakDetection;
import org.apache.commons.cli.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class RunExecution {
    public static void main(String[] args) {
        Options options = new Options();
        CommandLine commandLine = null;
        HelpFormatter help = new HelpFormatter();
        String header = "renlabm6a_allele.jar provides the following tools: ";
        String footer = "";

        try {
            commandLine = setCommandLine(options, args);
        } catch (ParseException pe) {
            System.err.println(pe.getMessage());
            help.printHelp("java -jar renlabm6a_allele.jar", header, options, footer, true);
            System.exit(2);
        }

        if (!checkArguments(args)) {
            System.err.println("tools can not be used simultaneously");
            help.printHelp("java -jar renlabm6a_allele.jar", header, options, footer, true);
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
            help.printHelp("java -jar renlabm6a_allele.jar", header, options, footer, true);
            System.exit(0);
        } else {
            help.printHelp("java -jar renlabm6a_allele.jar", header, options, footer, true);
            System.exit(0);
        }

        execute(runMode, args);
    }

    private static boolean checkArguments(String[] args) {
        List list = Arrays.asList(args);
        int[] useTools = new int[3];
        useTools[0] = list.contains("-AseGeneDetection")? 1: 0;
        useTools[1] = list.contains("-AsmPeakDetection")? 1: 0;
        useTools[2] = list.contains("-AseSeqSimulator")? 1: 0;
        list = null;
        int tools = 0;
        for (int i: useTools) {
            tools += i;
        }

        return tools <= 1;
    }

    private static void execute(int runMode, String[] args) {
        String[] arr;
        if (runMode == 1) {
            arr = delRedundantOption(args, "-AseGeneDetection");
            AseGeneDetection.main(arr);
        } else if (runMode == 2) {
            arr = delRedundantOption(args, "-AsmPeakDetection");
            AsmPeakDetection.main(arr);
        } else {
            arr = delRedundantOption(args, "-AseSeqSimulator");
            AseSeqSimulator.main(arr);
        }
    }

    public static String[] delRedundantOption(String[] args, String targetOption) {
        ArrayList<String> list = new ArrayList<>(Arrays.asList(args));
        list.remove(targetOption);
        String[] arr = new String[list.size()];
        arr = list.toArray(arr);
        list = null;

        return arr;

    }

    private static CommandLine setCommandLine(Options options, String[] args) throws ParseException {
        Option option = new Option("AseGeneDetection", false,"detect allele-specific expression genes in test data");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("AsmPeakDetection", false, "detect allele-specific modification m6A signals in test data");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("AseSeqSimulator", false, "generate simulation data for test");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("h", "help", false, "show help message and exit program");
        option.setRequired(false);
        options.addOption(option);

        CommandLineParser parser = new DefaultParser();

        return parser.parse(options, args);
    }
}
