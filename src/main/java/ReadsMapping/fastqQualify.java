package ReadsMapping;

import java.io.*;
import java.util.ArrayList;

public class fastqQualify {

    public static void main(String[] args) {
        runFastqc("/home/hbs/fastqc", "/home/hbs", 8);
    }

    /**
     * run FASTQC software with command "fastqc -t 8 -o path/fastqc sample1_R1.fq sample1_R2.fq"
     * @param fastqcDir Output result directory of FASTQC
     * @param rawFastqDir directory of FASTQC files
     * @param execThread number of working threads
     */
    public static void runFastqc(String fastqcDir, String rawFastqDir, int execThread) {
        File fastqcResultDir = new File(fastqcDir);
        File rawFastqDataDir = new File(rawFastqDir);

        if (!rawFastqDataDir.exists() | rawFastqDataDir.isFile()) {
            throw new RuntimeException("param must be the directory of fastq files");
        }
        File[] fastqFiles = rawFastqDataDir.listFiles();
        StringBuilder fastqList = new StringBuilder();
        try {
            for (File f : fastqFiles) {
                if (f.getName().endsWith(".fastq")) {
                    String fPath = " " + f.getAbsolutePath();
                    fastqList.append(fPath);
                }
            }
        } catch (NullPointerException ne) {
            System.out.println("empty fastq directory.");
        }

        String fq = fastqList.toString();

        if (!fastqcResultDir.exists()) {
            boolean res = fastqcResultDir.mkdir();
            if (!res) {
                throw new RuntimeException("can't make output directory");
            }
        }

        String cmd = "fastqc -t " + execThread + " -o " + fastqcDir + fq;

        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                // output the fail reason into log file
                BufferedReader record = new BufferedReader(new InputStreamReader(p.getInputStream()));
                ArrayList<String> recordLine = new ArrayList<String>();
                String line;
                while ((line = record.readLine()) != null) {
                    recordLine.add(line);
                }
                String logFileName = new File(fastqcDir, "fail.log").getAbsolutePath();
                FileWriter fw = new FileWriter(logFileName);
                for (String l : recordLine) {
                    fw.write(l);
                }
                fw.close();
                throw new RuntimeException("fastqc process failed. Check reason in " + logFileName);
            }
        } catch (IOException | InterruptedException ie) {
            ie.printStackTrace();
        }
    }
}
