package meripSeqPeakCalling;

import org.apache.commons.cli.*;
import org.apache.log4j.Logger;
import org.rosuda.REngine.REXP;
import org.rosuda.REngine.Rserve.RConnection;
import org.rosuda.REngine.Rserve.RserveException;

import java.io.*;

/**
 * 通过RServe使Java与R进行交互，使用MeTPeak包进行peak calling
 */
public class PeakCaller {
    private String gtfFile, ipBam, inputBam, outputDir, experimentName;
    private Logger logger;

    /**
     * Constructor
     * @param gtfFile GTF文件
     * @param ipBam IP样本比对结果
     * @param inputBam  INPUT样本比对结果
     * @param outputDir 输出文件目录
     * @param experimentName 实验名称
     * @param logger log4j对象
     */
    public PeakCaller(String gtfFile, String ipBam, String inputBam, String outputDir, String experimentName, Logger logger) {
        this.gtfFile = gtfFile;
        this.ipBam = ipBam;
        this.inputBam = inputBam;
        this.outputDir = outputDir;
        this.experimentName = experimentName;
        this.logger = logger;
    }

    /**
     * entrance of the program
     * @param args arguments
     * @throws ParseException parse exception
     */
    public static void main(String[] args) throws ParseException{
        Options options = new Options();
        CommandLine commandLine = setCommandline(args, options);
        String ip, input, gtfFile, outputDir = System.getProperty("user.dir"), experimentName = "m6aPeak";
        Logger log;

        ip = commandLine.getOptionValue("ip");
        input = commandLine.getOptionValue("input");
        gtfFile = commandLine.getOptionValue("g");
        if (commandLine.hasOption("o"))
            outputDir = commandLine.getOptionValue("o");
        if (commandLine.hasOption("exp"))
            experimentName = commandLine.getOptionValue("exp");

        log = initLog(outputDir);
        File output = new File(outputDir);
        if (!output.exists()) {
            boolean res = output.mkdirs();
            if (!res) {
                log.error("can not make directory " + outputDir);
                System.exit(2);
            }
        } else {
            if (output.isFile()) {
                log.error(output + " should be a directory not a file");
                System.exit(2);
            }
        }

        PeakCaller pc = new PeakCaller(gtfFile, ip, input, outputDir, experimentName, log);
        pc.peakCalling();
    }

    /**
     * peak calling with metpeak in R script
     */
    public void peakCalling() {

        File rScript = new File(this.outputDir, "peakCalling.R");
        boolean gen = generateRScript(rScript);
        if (!gen) {
            this.logger.error(" R script can not generate or set to be executable");
            System.exit(2);
        }
        try {
            Process p = Runtime.getRuntime().exec("chmod a+x " + rScript.getAbsolutePath());
            int res = p.waitFor();
            if (res != 0) {
                this.logger.error("chmod failed");
                System.exit(2);
            }
        } catch (IOException | InterruptedException ie) {
            this.logger.error(ie.getMessage());
            System.exit(2);
        }

        RConnection rc = null;
        REXP x = null;
        try {
            rc = new RConnection();
            rc.eval("source('" + rScript + "')");
            this.logger.debug("start peak calling...");
            String finalOutput = new File(this.outputDir, this.experimentName).getAbsolutePath();
            this.logger.debug("peak calling procedure, may take a little bit long time. Result output in: " + finalOutput);
            rc.eval("peakCall()");
            rc.close();
        } catch (RserveException e) {  // | REXPMismatchException
            this.logger.error(e.getMessage());
            System.exit(2);
        }
        this.logger.debug("peak calling complete");

        boolean res = rScript.delete();
        if (!res)
            this.logger.error("can not delete redundant file " + rScript);
    }

    /**
     * generate R script for MeTPeak peak calling
     * @param rScriptFilePath R script file path
     */
    private boolean generateRScript(File rScriptFilePath) {
        boolean res;
        try {
            if (!rScriptFilePath.exists()) {
                res = rScriptFilePath.createNewFile();
                if (!res) {
                    return res;
                }
            }

            BufferedWriter bfw = new BufferedWriter(
                    new OutputStreamWriter(new FileOutputStream(rScriptFilePath, false))
            );

            bfw.write("library(\"MeTPeak\")");
            bfw.newLine();
            bfw.write("peakCall <- function() {\n\tmetpeak(GENE_ANNO_GTF=\"" + this.gtfFile + "\",IP_BAM=\"" + this.ipBam + "\", INPUT_BAM=\"" + this.inputBam +
                           "\", OUTPUT_DIR=\"" + this.outputDir + "\", EXPERIMENT_NAME=\"" + this.experimentName + "\")\n}");
            bfw.newLine();
            bfw.flush();
            bfw.close();
        } catch (IOException ie) {
            ie.getMessage();
            System.exit(2);
        }
        res = rScriptFilePath.setExecutable(true);

        return res;
    }

    /**
     * parse command line parameters
     * @param args input parameters
     * @param options Options instance
     * @return CommandLineParser
     * @throws ParseException parseException
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
        return Logger.getLogger(PeakCaller.class);
    }
}
