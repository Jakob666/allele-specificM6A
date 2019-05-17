package AseSeqSimulator;

import org.apache.commons.cli.*;

import java.io.*;

public class AseSeqSimulator {

    public static void main(String[] args) throws Exception {

        CommandLine commandLine = AseSeqSimulator.parseCommandLine(args);

        int librarySize = 1000000, readLength = 50, fragmentMean = 250, fragmentStd = 25, multiple = 1,
            repeat = 1, minimumMut = 5, maximumMut = 15, peakLength = 500, depth = 0;
        // default 20% gene has SNP site;
        double geneProp = 0.2, mutProportion = 0.4, pcrErrorProb = 0.005;
        String gtfFile, twoBitFile, vcfFile = null, geneExpFile = null;
        boolean overlap = false, singleEnd = true;
        File outputDir = new File("./AseSeqReads");

        gtfFile = commandLine.getOptionValue('g');
        twoBitFile = commandLine.getOptionValue('t');
        if (commandLine.hasOption('v'))
            vcfFile = commandLine.getOptionValue('v');
        if (commandLine.hasOption("exp"))
            geneExpFile = commandLine.getOptionValue("exp");
        if (commandLine.hasOption('o'))
            outputDir = new File(commandLine.getOptionValue('o'));
        if (commandLine.hasOption("ls"))
            librarySize = Integer.parseInt(commandLine.getOptionValue("ls"));
        if (commandLine.hasOption("pl"))
            peakLength = Integer.parseInt(commandLine.getOptionValue("pl"));
        if (commandLine.hasOption("dep")) {
            depth = Integer.parseInt(commandLine.getOptionValue("dep"));
            if (depth < 0) {
                System.out.println("depth must larger than 0");
                System.exit(-1);
            }
        }
        if (commandLine.hasOption("rl"))
            readLength = Integer.parseInt(commandLine.getOptionValue("rl"));
        if (commandLine.hasOption("mi"))
            minimumMut = Integer.parseInt(commandLine.getOptionValue("mi"));
        if (commandLine.hasOption("mx"))
            maximumMut = Integer.parseInt(commandLine.getOptionValue("mx"));
        if (commandLine.hasOption("fm"))
            fragmentMean = Integer.parseInt(commandLine.getOptionValue("fm"));
        if (commandLine.hasOption("ft"))
            fragmentStd = Integer.parseInt(commandLine.getOptionValue("ft"));
        if (commandLine.hasOption("mul"))
            multiple = Integer.parseInt(commandLine.getOptionValue("mul"));
        if (commandLine.hasOption("rep"))
            repeat = Integer.parseInt(commandLine.getOptionValue("rep"));
        if (commandLine.hasOption("gp")) {
            double number = Double.parseDouble(commandLine.getOptionValue("gp"));
            if (number > 0)
                geneProp = number;
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
        if (commandLine.hasOption("ol"))
            overlap = Boolean.getBoolean(commandLine.getOptionValue("ol"));
        if (commandLine.hasOption("se"))
            singleEnd = Boolean.getBoolean(commandLine.getOptionValue("se"));
        if (commandLine.hasOption("pe")) {
            double pcrError = Double.parseDouble(commandLine.getOptionValue("pe"));
            if (pcrError > 0 && pcrError < 1)
                pcrErrorProb = pcrError;
            else {
                System.out.println("invalid input PCR sequencing error probability, must in range 0 to 1");
                System.exit(2);
            }
        }

        // if depth is assigned, generate reads
        if (depth != 0) {
            geneExpFile = null;
        }

        ReadsGenerator readsGenerator = new ReadsGenerator(gtfFile, geneProp, twoBitFile);
        readsGenerator.simulateSequencing(outputDir.getAbsolutePath(), vcfFile, librarySize, depth, peakLength, readLength, minimumMut,
                                          maximumMut, fragmentMean, fragmentStd, mutProportion, multiple, repeat, pcrErrorProb,
                                          overlap, singleEnd, geneExpFile);

    }

    private static CommandLine parseCommandLine(String[] args) throws ParseException {
        Options options = new Options();

        Option option = new Option("g", "gtf", true, "gtf file path");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("t", "twobit", true, "UCSC 2bit file path");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("v", "vcf_file", true, "vcf file path used for generate SNP, default null");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("exp", "gene_express", true, "cellLineExpFile downloading file from http://medicalgenomics.org/rna_seq_atlas/download");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("o", "outputDir", true, "output directory, default ./AseSeqRead");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("ls", "library_size", true, "cDNA library size, default 1000000");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("dep", "depth", true, "sequencing depth, default 0.");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("pl", "peak_length", true, "m6A peak length, default 500");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("rl", "read_length", true, "sequencing reads length, default 50");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("mi", "minimum_mutation", true, "minimum mutation sites on fragment, default 5");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("mx", "maximum_mutation", true, "maximum mutation sites on fragment, default 15");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("fm", "fragment_mean", true, "mean length of fragments, default 250");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("ft", "fragment_theta", true, "standard deviation of fragment length, default 25");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("gp", "gene_proportion", true, "The proportion of genes selected in the total number of genes on a particular chromosome, default 0.2");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("mp", "mutate_proportion", true, "The proportion of mutated genes in the total number of selected genes on a particular chromosome, default 0.4");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("ol", "overlap", true, "if the random select gene overlapped with each other, default false");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("se", "single_end", true, "if true, single-end reads, otherwise pair-end. Default true");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("pe", "pcr_error", true, "probability of PCR sequencing error, default 5‰");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("mul", "multiple_time", true, "multiple time, default 1");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("rep", "repeat_time", true, "experiment repeat time, default 1");
        option.setRequired(false);
        options.addOption(option);

        CommandLineParser parser = new DefaultParser();

        return parser.parse(options, args);
    }
}
