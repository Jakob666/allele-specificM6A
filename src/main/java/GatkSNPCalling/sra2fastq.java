package GatkSNPCalling;

import java.io.File;
import java.io.IOException;

public class sra2fastq {
    private String sraDataDir, fastqDataDir;

    /**
     * Constructor
     * @param rawDataDir directory stores sra raw data
     * @param fastqDir dorectory stores fastq data
     */
    sra2fastq(String rawDataDir, String fastqDir) {
        this.sraDataDir = rawDataDir;
        this.fastqDataDir = fastqDir;
    }

    /**
     * trans sra format to fastq format for files under a particular directory.
     * the IP and INPUT sra data store in sub-directory "IP" and "INPUT" respectively.
     * @throws RuntimeException
     */
    public boolean transFormat() throws RuntimeException {
        boolean execSuccess;

        File outputFastqDir = new File(this.fastqDataDir);
        if (!outputFastqDir.exists()) {
            makeDir(outputFastqDir);
        }

        File rawData = new File(this.sraDataDir);
        if (rawData.isFile()) {
            throw new RuntimeException("invalid sra data directory");
        }
        String cellLineName = rawData.getName();

        File outputDir = new File(this.fastqDataDir, cellLineName);
        File ipOutputDir = new File(outputDir.getAbsolutePath(), "IP");
        File inputOutputDir = new File(outputDir.getAbsolutePath(), "INPUT");

        if (!outputDir.exists()) {
            makeDir(outputDir);
            makeDir(ipOutputDir);
            makeDir(inputOutputDir);
        }

        execSuccess = sra2fastq(rawData, outputDir);
        return execSuccess;
    }

    /**
     * Transform SRA file in a directory to FASTQ format
     * @param sraFileDir directory of sra files
     * @param fastqFileDir directory of fastq files
     */
    public boolean sra2fastq(File sraFileDir, File fastqFileDir) {
        // SRA raw data stores in sub-directory "IP" and "INPUT"
        File ipDir = new File(sraFileDir.getAbsolutePath(), "IP");
        File inputDir = new File(sraFileDir.getAbsolutePath(), "INPUT");

        File[] ipFiles = ipDir.listFiles();
        File[] inputFiles = inputDir.listFiles();

        File inputFastqDir = new File(fastqFileDir.getAbsolutePath(), "INPUT");
        File ipFastqDir = new File(fastqFileDir.getAbsolutePath(), "IP");

        try {
            for (File f : ipFiles) {
                if (f != null && f.isFile() && f.getName().endsWith(".sra")) {
                    String sraFilePath = f.getAbsolutePath();
                    transform(sraFilePath, ipFastqDir.getAbsolutePath());
                }
            }

            for (File f : inputFiles) {
                if (f != null && f.isFile() && f.getName().endsWith(".sra")) {
                    String sraFilePath = f.getAbsolutePath();
                    transform(sraFilePath, inputFastqDir.getAbsolutePath());
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
    private static void transform(String sraFilePath, String fastqDir) {
        String cmd = "fastq-dump " + sraFilePath + " -O " + fastqDir;
        System.out.println(cmd);
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                throw new RuntimeException("fail to transform sra file " + sraFilePath);
            }
        }catch (IOException | InterruptedException ie) {
            ie.printStackTrace();
        }
    }

    private static void makeDir(File mkDirectory) {
        boolean res = mkDirectory.mkdir();
        if (!res) {
            System.out.println("fail to make directory " + mkDirectory.getAbsolutePath());
            System.exit(2);
        }
    }
}
