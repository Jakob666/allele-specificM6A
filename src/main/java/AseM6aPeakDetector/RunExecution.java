package AseM6aPeakDetector;

import AseSeqSimulator.AseSeqSimulator;
import DifferentiationAnalysis.SampleSpecificASE;
import DifferentiationAnalysis.SampleSpecificASM;
import HierarchicalBayesianAnalysis.AseGeneDetection;
import HierarchicalBayesianAnalysis.AsmPeakDetection;

import java.util.ArrayList;
import java.util.Arrays;

public class RunExecution {
    public static void main(String[] args) {
        String version = "version 1.0";
        String packageUsage = "renlabm6a_allele.jar provides the following tools: \n" +
                "usage: java -jar renlabm6a_allele.jar [AseGeneDetection] [AseSeqSimulator] [AsmPeakDetection] [-h]\n" +
                "renlabm6a_allele.jar provides the following tools:\n" +
                " AseGeneDetection   detect allele-specific expression (ASE) genes (one sample test)\n" +
                " AsmPeakDetection   detect allele-specific modification (ASM) m6A signals  (one sample test)\n" +
                " SampleSpecificASE  detect sample-specific ASE genes  (two sample test)\n" +
                " SampleSpecificASM  detect sample-specific ASM m6A signals  (two sample test)\n" +
                " AseSeqSimulator    generate simulation data for test\n" +
                " -h,--help          show help message and exit program\n" +
                " -v,--version       release version\n\n" + version;

        ArrayList<String> argsArr = new ArrayList<>(Arrays.asList(args));

        if (!checkArguments(argsArr)) {
            System.err.println("tools can not be used simultaneously");
            System.err.println(packageUsage);
            System.exit(2);
        }

        int runMode = 0;

        if (argsArr.contains("AseGeneDetection"))
            runMode = 1;
        else if (argsArr.contains("AsmPeakDetection"))
            runMode = 2;
        if (argsArr.contains("SampleSpecificASE"))
            runMode = 3;
        else if (argsArr.contains("SampleSpecificASM"))
            runMode = 4;
        else if (argsArr.contains("AseSeqSimulator"))
            runMode = 5;
        else if (argsArr.contains("-h") | argsArr.contains("--help")) {
            System.out.println(packageUsage);
            System.exit(0);
        } else if (argsArr.contains("-v") | argsArr.contains("--version")) {
            System.out.println(version);
            System.exit(0);
        } else {
            System.err.println(packageUsage);
            System.exit(0);
        }

        execute(runMode, argsArr);
    }

    private static boolean checkArguments(ArrayList<String> list) {
        int[] useTools = new int[5];
        useTools[0] = list.contains("AseGeneDetection")? 1: 0;
        useTools[1] = list.contains("AsmPeakDetection")? 1: 0;
        useTools[2] = list.contains("SampleSpecificASE")? 1: 0;
        useTools[3] = list.contains("SampleSpecificASM")? 1: 0;
        useTools[4] = list.contains("AseSeqSimulator")? 1: 0;
        int tools = Arrays.stream(useTools).reduce((x, y) -> x+y).getAsInt();

        return tools <= 1;
    }

    private static void execute(int runMode, ArrayList<String> args) {
        String[] arr;
        if (runMode == 1) {
            arr = delRedundantOption(args, "AseGeneDetection");
            AseGeneDetection.main(arr);
        } else if (runMode == 2) {
            arr = delRedundantOption(args, "AsmPeakDetection");
            AsmPeakDetection.main(arr);
        } else if (runMode == 3) {
            arr = delRedundantOption(args, "SampleSpecificASE");
            SampleSpecificASE.main(arr);
        } else if (runMode == 4) {
            arr = delRedundantOption(args, "SampleSpecificASM");
            SampleSpecificASM.main(arr);
        } else {
            arr = delRedundantOption(args, "AseSeqSimulator");
            AseSeqSimulator.main(arr);
        }
    }

    public static String[] delRedundantOption(ArrayList<String> list, String targetOption) {
        list.remove(targetOption);
        String[] arr = new String[list.size()];
        arr = list.toArray(arr);
        list = null;

        return arr;

    }
}
