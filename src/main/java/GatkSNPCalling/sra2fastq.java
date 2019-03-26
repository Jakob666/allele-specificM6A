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
    public sra2fastq(String rawDataDir, String fastqDir) {
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

        File outputDir = new File(this.fastqDataDir);
        if (!outputDir.exists()) {
            makeDir(outputDir);
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
        File[] sraFiles = sraFileDir.listFiles();
        if (sraFiles == null)
            return false;
        try {
            for (File f : sraFiles) {
                if (f != null && f.isFile() && f.getName().endsWith(".sra")) {
                    String sraFilePath = f.getAbsolutePath();
                    transform(sraFilePath, fastqFileDir.getAbsolutePath());
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
