package CommonThreadClass;

import org.apache.log4j.Logger;

import java.io.*;
import java.util.ArrayList;

public class Sra2Fastq {
    private static String sraTool, sraDataDir, fastqDir, direct;
    private static Logger log;
    private static Sra2Fastq sf = null;
    private Sra2Fastq(String sraToolkit, String sourceDataDir, String fastqTempDir, String seqDirect, Logger logger) {
        sraTool = sraToolkit;
        sraDataDir = sourceDataDir;
        fastqDir = fastqTempDir;
        direct = seqDirect;
        log = logger;
    }

    public static Sra2Fastq getInstance(String sraToolkit, String sourceDataDir, String fastqTempDir, String direct, Logger logger) {
        if (sf == null) {
            sf = new Sra2Fastq(sraToolkit, sourceDataDir, fastqTempDir, direct, logger);
        }
        return sf;
    }

    /**
     * 获取目录下的SRA文件
     * @return 返回全部SRA文件路径
     */
    public ArrayList<String> getSraFiles() {
        File sraSourceDir = new File(sraDataDir);
        File[] files = sraSourceDir.listFiles();
        if (files == null) {
            log.error("empty SRA directory");
            System.exit(2);
        }
        ArrayList<String> sraFiles = new ArrayList<>();
        for (File file: files) {
            log.debug(file.getAbsolutePath());
            if (file.getName().endsWith("sra"))
                sraFiles.add(file.getAbsolutePath());
        }
        if (sraFiles.size() == 0) {
            log.error("directory contains no SRA files");
            System.exit(2);
        }
        return sraFiles;
    }

    /**
     * 将原始数据目录下的SRA文件转为Fastq的gz压缩格式
     */
    public void transform() {
        ArrayList<String> sraFiles = getSraFiles();
        String cmd = null;
        for (String sraFile: sraFiles) {
            this.switchToFastq(sraFile);
        }
    }

    public void switchToFastq(String sraFile) {
        String cmd;
        try {
            if (direct.toLowerCase().equals("se"))
                cmd = String.join(" ", new String[]{sraTool+"fastq-dump", "-gzip", sraFile, "-O", fastqDir});
            else
                cmd = String.join(" ", new String[]{sraTool+"fastq-dump", "-gzip", "--split-3", sraFile, "-O", fastqDir});

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
