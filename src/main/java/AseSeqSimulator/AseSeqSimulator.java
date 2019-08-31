package AseSeqSimulator;

import org.apache.commons.cli.*;

import java.io.*;

public class AseSeqSimulator {

    public static void main(String[] args) throws Exception {

        CommandLine commandLine = AseSeqSimulator.parseCommandLine(args);

        int librarySize = 10000000, readLength = 75, fragmentMean = 250, fragmentStd = 25, peakLength = 250,
            minimumMut = 0, maximumMut = 3, depth = 100, minReadsCoverage = 10, maxReadsCoverage = 70;
        // default 60% gene has SNP site;
        double geneProp = 0.05, mutProportion = 0.6, pcrErrorProb = 0, aseInfimum = 0.55, aseSupremum = 0.85;
        String gtfFile, twoBitFile, vcfFile = null, geneExpFile = null;
        boolean overlap = true, singleEnd = true;
        File outputDir = new File("./simulateData");

        gtfFile = commandLine.getOptionValue('g');
        twoBitFile = commandLine.getOptionValue('t');
        if (commandLine.hasOption('v'))
            vcfFile = commandLine.getOptionValue('v');
        if (commandLine.hasOption("exp"))
            geneExpFile = commandLine.getOptionValue("exp");
        if (commandLine.hasOption('o'))
            outputDir = new File(commandLine.getOptionValue('o'));
        if (!outputDir.exists())
            outputDir.mkdir();
        if (commandLine.hasOption("ls"))
            librarySize = Integer.parseInt(commandLine.getOptionValue("ls"));
        if (commandLine.hasOption("dep")) {
            depth = Integer.parseInt(commandLine.getOptionValue("dep"));
            if (depth < 0) {
                System.out.println("depth must larger than 0");
                System.exit(-1);
            }
        }
        if (commandLine.hasOption("rl"))
            readLength = Integer.parseInt(commandLine.getOptionValue("rl"));
        if (commandLine.hasOption("pl"))
            peakLength = Integer.parseInt(commandLine.getOptionValue("pl"));
        if (commandLine.hasOption("min_mut"))
            minimumMut = Integer.parseInt(commandLine.getOptionValue("min_mut"));
        if (commandLine.hasOption("max_mut"))
            maximumMut = Integer.parseInt(commandLine.getOptionValue("max_mut"));
        if (commandLine.hasOption("min_cover"))
            minReadsCoverage = Integer.parseInt(commandLine.getOptionValue("min_cover"));
        if (commandLine.hasOption("max_cover"))
            maxReadsCoverage = Integer.parseInt(commandLine.getOptionValue("max_cover"));
        if (commandLine.hasOption("al"))
            aseInfimum = Double.parseDouble(commandLine.getOptionValue("al"));
        if (commandLine.hasOption("ah"))
            aseSupremum = Double.parseDouble(commandLine.getOptionValue("ah"));
        if (aseInfimum > aseSupremum || aseInfimum >= 1 || aseSupremum >= 1) {
            System.out.println("invalid ASE ratio");
            System.exit(2);
        }
        if (commandLine.hasOption("fm"))
            fragmentMean = Integer.parseInt(commandLine.getOptionValue("fm"));
        if (commandLine.hasOption("ft"))
            fragmentStd = Integer.parseInt(commandLine.getOptionValue("ft"));
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

        ReadsGenerator readsGenerator = new ReadsGenerator(gtfFile, geneProp, aseInfimum, aseSupremum, twoBitFile);
        readsGenerator.simulateSequencing(outputDir.getAbsolutePath(), vcfFile, librarySize, depth, readLength, minimumMut,
                                          maximumMut, fragmentMean, fragmentStd, minReadsCoverage, maxReadsCoverage,
                                          mutProportion, peakLength, pcrErrorProb, overlap, singleEnd, geneExpFile);
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

        option = new Option("dep", "depth", true, "sequencing depth, default 100.");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("rl", "read_length", true, "sequencing reads length, default 75");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("pl", "peak_length", true, "m6A signal peak length, default 250");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("min_mut", "minimum_mutation", true, "minimum mutation sites number on fragment under a m6A signal, default 0");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("max_mut", "maximum_mutation", true, "maximum mutation sites number on fragment under a m6A signal, default 3");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("min_cover", "minimum_coverage", true, "minimum reads coverage when generate RNA-seq data, default 10");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("max_cover", "maximum_coverage", true, "maximum reads coverage when generate RNA-seq data, default 70");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("al", "ase_infimum", true, "ASE ratio infimum, default 0.55");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("ah", "ase_supremum", true, "ASE ratio supremum, default 0.85");
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

        CommandLineParser parser = new DefaultParser();

        return parser.parse(options, args);
    }
}
