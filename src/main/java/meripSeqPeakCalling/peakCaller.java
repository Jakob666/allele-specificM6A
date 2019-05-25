package meripSeqPeakCalling;

import org.apache.log4j.Logger;

import java.io.*;

public class peakCaller {

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
        String[] params = new String[]{"nohup Rscript", rScript.getAbsolutePath(),"&"};
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
}
