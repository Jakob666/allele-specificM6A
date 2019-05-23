package GatkSNPCalling;

import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;

public class sra2fastq {
    private String sraDataDir, fastqDataDir;
    private Logger logger;

    /**
     * Constructor
     * @param rawDataDir directory stores sra raw data
     * @param fastqDir dorectory stores fastq data
     */
    public sra2fastq(String rawDataDir, String fastqDir, Logger logger) {
        this.sraDataDir = rawDataDir;
        this.fastqDataDir = fastqDir;
        this.logger = logger;
    }

    /**
     * trans sra format to fastq format for files under a particular directory.
     * the IP and INPUT sra data store in sub-directory "IP" and "INPUT" respectively.
     * @throws RuntimeException
     */
    public boolean transFormat() {
        boolean execSuccess;

        // make up fastq data directory if it not existed
        File outputFastqDir = new File(this.fastqDataDir);
        if (outputFastqDir.exists() && outputFastqDir.isFile()) {
            this.logger.error("invalid fastq data directory " + this.fastqDataDir + ", need a directory not a file");
            System.exit(2);
        }
        if (!outputFastqDir.exists()) {
            this.makeDir(outputFastqDir);
        }

        File rawData = new File(this.sraDataDir);
        if (rawData.isFile()) {
            this.logger.error("invalid sra data directory " + this.sraDataDir + ", need a directory not a file");
            System.exit(2);
        }

        File outputDir = new File(this.fastqDataDir);
        if (!outputDir.exists()) {
            this.makeDir(outputDir);
        }

        execSuccess = sra2fastq(rawData, outputDir);
        return execSuccess;
    }

    /**
     * Transform SRA file in a directory to FASTQ format
     * @param sraFileDir directory of sra files
     * @param fastqFileDir directory of fastq files
     */
    private boolean sra2fastq(File sraFileDir, File fastqFileDir) {
        File[] sraFiles = sraFileDir.listFiles();
        if (sraFiles == null) {
            this.logger.error("empty raw data directory " + sraFileDir.getAbsolutePath());
            return false;
        }

        try {
            for (File f : sraFiles) {
                if (f != null && f.isFile() && f.getName().endsWith(".sra")) {
                    String sraFilePath = f.getAbsolutePath();
                    this.transform(sraFilePath, fastqFileDir.getAbsolutePath());
                }
            }
        } catch (RuntimeException re) {
            return false;
        }

        return true;
    }

    /**
     * Transform a single sra file to FASTQ file with command "fastq-dump SRRxxxx.sra
     * @param sraFilePath the absolute path for the sra file
     * @param fastqDir directory of fastq files
     */
    private void transform(String sraFilePath, String fastqDir) {
        String sraFileName = new File(sraFilePath).getName();
        String targetFastqFileName = sraFileName.substring(0, sraFileName.lastIndexOf(".")) + ".fastq";

        if (new File(fastqDir, targetFastqFileName).exists()) {
            logger.debug("already existed file " + targetFastqFileName);
            return;
        }

        String cmd = "fastq-dump " + sraFilePath + " -O " + fastqDir;
        this.logger.debug("transform sra file " + sraFilePath + ", output in " + fastqDir);
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                this.logger.error("fail to transform sra file " + sraFilePath);
            }
            this.logger.debug("complete");
        }catch (IOException | InterruptedException ie) {
            this.logger.error(ie.getMessage());
        }
    }

    private void makeDir(File mkDirectory) {
        boolean res = mkDirectory.mkdir();
        if (!res) {
            logger.error("fail to make directory " + mkDirectory.getAbsolutePath());
            System.exit(2);
        }
    }
}
