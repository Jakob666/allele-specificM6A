package meripSeqPeakCalling;

import org.apache.commons.cli.*;
import org.apache.log4j.Logger;

import java.io.*;

public class peakCaller {

    /**
     * entrance of the program
     * @param args arguments
     * @throws ParseException parse exception
     */
    public static void main(String[] args) throws ParseException{
        Options options = new Options();
        CommandLine commandLine = setCommandline(args, options);
        String ip, input, gtfFile, outputDir = System.getProperty("user.dir"), experimentName = "m6aPeak";
        Logger log = null;

        ip = commandLine.getOptionValue("ip");
        input = commandLine.getOptionValue("input");
        gtfFile = commandLine.getOptionValue("g");
        if (commandLine.hasOption("o"))
            outputDir = commandLine.getOptionValue("o");
        if (commandLine.hasOption("exp"))
            experimentName = commandLine.getOptionValue("exp");

        log = initLog(outputDir);
        peakCalling(gtfFile, outputDir, ip, input, experimentName, log);
    }

    /**
     * peak calling with metpeak in R script
     * @param gtfFilePath GTF file path
     * @param outputDirPath name of output result directory
     */
    public static void peakCalling(String gtfFilePath, String outputDirPath, String IPFile, String INPUTFile, String experimentName , Logger log) {

        File rScript = new File(outputDirPath, "peakCalling.R");
        boolean gen = generateRScript(rScript, gtfFilePath, IPFile, INPUTFile, experimentName);
        if (!gen) {
            log.error(" R script can not generate or set to be executable");
            System.exit(2);
        }
        try {
            Process p = Runtime.getRuntime().exec("chmod a+x " + rScript.getAbsolutePath());
            int res = p.waitFor();
            if (res != 0) {
                log.error("chmod failed");
                System.exit(2);
            }
        } catch (IOException | InterruptedException ie) {
            log.error(ie.getMessage());
            System.exit(2);
        }
        String[] params = new String[]{"Rscript", rScript.getAbsolutePath()};
        String cmd = String.join(" ", params);
        log.debug(cmd);
        String finalOutput = new File(outputDirPath, experimentName).getAbsolutePath();
        log.debug("peak calling procedure, may take a little bit long time. Result output in: " + finalOutput);
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int res = p.waitFor();
            if (res != 0) {
                log.error("Run peak calling failed");
                System.exit(2);
            }
        } catch (IOException | InterruptedException ie) {
            log.error(ie.getMessage());
            System.exit(2);
        }
        log.debug("peak calling complete");
    }

    /**
     * generate R script for MeTPeak peak calling
     * @param rScript R script file path
     * @param gtfFilePath GTF file path
     * @param IPFile IP bam file
     * @param INPUTFile INPUT bam file
     * @param experimentName output name
     */
    private static boolean generateRScript(File rScript, String gtfFilePath, String IPFile, String INPUTFile, String experimentName) {
        boolean res;
        try {
            if (!rScript.exists()) {
                res = rScript.createNewFile();
                if (!res) {
                    return res;
                }
            }

            BufferedWriter bfw = new BufferedWriter(
                    new OutputStreamWriter(new FileOutputStream(rScript, false))
            );

            bfw.write("library(\"MeTPeak\")");
            bfw.newLine();
            bfw.write("metpeak(GENE_ANNO_GTF=\"" + gtfFilePath + "\",IP_BAM=\"" + IPFile + "\", INPUT_BAM=\"" + INPUTFile +
                           "\", OUTPUT_DIR=\"" + rScript.getParent() + "\", EXPERIMENT_NAME=\"" + experimentName + "\")");
            bfw.newLine();
            bfw.flush();
            bfw.close();
        } catch (IOException ie) {
            ie.getMessage();
            System.exit(2);
        }
        res = rScript.setExecutable(true);

        return res;
    }

    /**
     * parse command line parameters
     * @param args input parameters
     * @param options Options instance
     * @return CommandLineParser
     * @throws ParseException
     */
    private static CommandLine setCommandline(String[] args, Options options) throws ParseException {

        Option option = new Option("ip", "ip_bam", true, "IP sample alignment bam file path");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("input", "input_bam", true, "INPUT sample alignment bam file path");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("g", "gtf", true, "GTF file path");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("o", "output_dir", true, "output directory, default present working directory");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("exp", "experiment_name", true, "experiment name, default 'm6aPeak'");
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
        return Logger.getLogger(peakCaller.class);
    }
}
