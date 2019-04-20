package AseSeqSimulator;

import org.apache.commons.cli.*;

import java.io.*;

public class AseSeqSimulator {

    public static void main(String[] args) throws Exception {

        CommandLine commandLine = AseSeqSimulator.parseCommandLine(args);

        int librarySize = 100000, readLength = 50, fragmentMean = 250, fragmentStd = 25, geneNum = 0, multiple = 1, repeat = 1, mutateGeneNum;
        // default 20% gene has SNP site;
        double mutProportion = 0.2, pcrErrorProb = 0.005;
        String gtfFile, twoBitFile;
        File outputDir = new File("./AseSeqReads");

        gtfFile = commandLine.getOptionValue('g');
        twoBitFile = commandLine.getOptionValue('t');
        if (commandLine.hasOption('o'))
            outputDir = new File(commandLine.getOptionValue('o'));
        if (commandLine.hasOption("ls"))
            librarySize = Integer.parseInt(commandLine.getOptionValue("ls"));
        if (commandLine.hasOption("rl"))
            readLength = Integer.parseInt(commandLine.getOptionValue("rl"));
        if (commandLine.hasOption("fm"))
            fragmentMean = Integer.parseInt(commandLine.getOptionValue("fm"));
        if (commandLine.hasOption("ft"))
            fragmentStd = Integer.parseInt(commandLine.getOptionValue("ft"));
        if (commandLine.hasOption("mul"))
            multiple = Integer.parseInt(commandLine.getOptionValue("mul"));
        if (commandLine.hasOption("rep"))
            repeat = Integer.parseInt(commandLine.getOptionValue("rep"));
        if (commandLine.hasOption("gn")) {
            int number = Integer.parseInt(commandLine.getOptionValue("gn"));
            if (number > 0)
                geneNum = number;
            else {
                System.out.println("invalid input gene number, must large than 0");
                System.exit(2);
            }
        }
        if (commandLine.hasOption("mp")) {
            double proportion = Double.parseDouble(commandLine.getOptionValue("mp"));
            if (proportion > 0 && proportion < 1)
                mutProportion = proportion;
            else {
                System.out.println("invalid input mutate gene proportion, must in range 0 to 1");
                System.exit(2);
            }
        }
        if (commandLine.hasOption("pe")) {
            double pcrError = Double.parseDouble(commandLine.getOptionValue("pe"));
            if (pcrError > 0 && pcrError < 1)
                pcrErrorProb = pcrError;
            else {
                System.out.println("invalid input PCR sequencing error probability, must in range 0 to 1");
                System.exit(2);
            }
        }

        int gtfGeneNum = AseSeqSimulator.totalGeneNum(gtfFile);
        if (gtfGeneNum == 0) {
            System.out.println("invalid format of GTF file");
            System.exit(2);
        } else if (geneNum == 0) {
            geneNum = gtfGeneNum;
        } else if (geneNum > gtfGeneNum) {
            geneNum = gtfGeneNum;
        }

        mutateGeneNum = (int) (geneNum * mutProportion);

        ReadsGenerator readsGenerator = new ReadsGenerator(gtfFile, geneNum, twoBitFile);
        readsGenerator.simulateSequencing(outputDir.getAbsolutePath(), librarySize, readLength, fragmentMean, fragmentStd,
                                          mutateGeneNum, multiple, repeat, pcrErrorProb);

    }

    private static CommandLine parseCommandLine(String[] args) throws ParseException {
        Options options = new Options();

        Option option = new Option("g", "gtf", true, "gtf file path");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("t", "twobit", true, "UCSC 2bit file path");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("o", "outputDir", true, "output directory");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("ls", "library_size", true, "cDNA library size");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("rl", "read_length", true, "sequencing reads length");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("fm", "fragment_mean", true, "mean length of fragments");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("ft", "fragment_theta", true, "standard deviation of fragment length");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("gn", "gene_num", true, "numbers of gene select from total genes");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("mp", "mutate_proportion", true, "mutated gene proportion");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("pe", "pcr_error", true, "probability of PCR sequencing error");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("mul", "multiple_time", true, "multiple time");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("rep", "repeat_time", true, "experiment repeat time");
        option.setRequired(false);
        options.addOption(option);

        CommandLineParser parser = new DefaultParser();

        return parser.parse(options, args);
    }

    private static int totalGeneNum(String gtfFile) {
        BufferedReader bfr = null;
        int geneNum = 0;
        try {
            bfr = new BufferedReader(
                    new InputStreamReader(new FileInputStream(new File(gtfFile)))
            );
            String line = "";
            String[] info;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    info = line.split("\t");
                    if (info[2].equals("gene"))
                        geneNum ++;
                }
            }
            bfr.close();
        } catch (IOException ie) {
            ie.printStackTrace();
        } finally {
            if (bfr != null) {
                try {
                    bfr.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }

        return geneNum;
    }
}
