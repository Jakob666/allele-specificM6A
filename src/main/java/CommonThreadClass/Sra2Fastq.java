package CommonThreadClass;

import org.apache.log4j.Logger;

import java.io.*;
import java.util.ArrayList;

public class Sra2Fastq {
    private String sraTool, fastqDir, direct;
    private ArrayList<File> sraFiles;
    private Logger log;

    public Sra2Fastq(String sraTool, String fastqTempDir, String direct, ArrayList<File> sraFiles, Logger logger) {
        this.sraTool = sraTool;
        this.fastqDir = fastqTempDir;
        this.direct = direct;
        this.sraFiles = sraFiles;
        this.log = logger;
    }

    public ArrayList<String> transformat() {
        ArrayList<String> fastqFiles = new ArrayList<>();
        String sraFileName, sraFilePath, fastqFileName, fastqFilePath;
        for (File f: this.sraFiles) {
            sraFileName = f.getName();
            sraFilePath = f.getAbsolutePath();
            fastqFileName = sraFileName.substring(0, sraFileName.lastIndexOf(".sra")) + ".fastq.gz";
            fastqFilePath = new File(this.fastqDir, fastqFileName).getAbsolutePath();
            this.switchToFastq(sraFilePath);
            fastqFiles.add(fastqFilePath);
        }

        return fastqFiles;
    }

    public void switchToFastq(String sraFile) {
        String cmd;
        try {
            if (direct.toLowerCase().equals("se"))
                cmd = String.join(" ", new String[]{this.sraTool+"fastq-dump", "-gzip", sraFile, "-O", fastqDir});
            else
                cmd = String.join(" ", new String[]{this.sraTool+"fastq-dump", "-gzip", "--split-3", sraFile, "-O", fastqDir});

            log.debug(cmd);
            log.debug("transform SRA file: " + sraFile + " into Fastq format. Output directory: " + fastqDir);
            Process p = Runtime.getRuntime().exec(cmd);
            int res = p.waitFor();
            if (res != 0) {
                log.error("Fail to transform file: " + sraFile);
                System.exit(2);
            }
        } catch (IOException | InterruptedException ie) {
            log.error(ie.getMessage());
            System.exit(2);
        }
    }
}
