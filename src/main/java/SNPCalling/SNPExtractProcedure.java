package SNPCalling;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

public class SNPExtractProcedure {

    public static void main(String[] args) {
        String sraDataDir = args[0];
        String fastqTempDir = args[1];
        String alignmentOutputDir = args[2];
        String genomeFile = args[3];
        String picardLocalJar = args[4];
        String gatkLocalJar = args[5];
        int execThread = Integer.parseInt(args[6]);

        File sraDir = new File(sraDataDir);
        String cellLineName = sraDir.getName();

        // make directories for fastq files and alignment result
        mkDir(fastqTempDir);
        mkDir(alignmentOutputDir);

        boolean sraTransRes = sraToFastq(sraDataDir, fastqTempDir);
        if (!sraTransRes) {
            System.out.println("transform failed");
        }

        File fastqDir = new File(fastqTempDir, cellLineName);
        File fastqIPDir = new File(fastqDir.getAbsolutePath(), "IP");
        File fastqINPUTDir = new File(fastqDir.getAbsolutePath(), "INPUT");

        File alignDir = new File(alignmentOutputDir, cellLineName);
        mkDir(alignDir.getAbsolutePath());
        File alignIPDir = new File(alignDir.getAbsolutePath(), "IP");
        File alignINPUTDir = new File(alignDir.getAbsolutePath(), "INPUT");

        SNPCaller.snpCalling(genomeFile, fastqIPDir.getAbsolutePath(), alignIPDir.getAbsolutePath(), picardLocalJar, gatkLocalJar, execThread);
        SNPCaller.snpCalling(genomeFile, fastqINPUTDir.getAbsolutePath(), alignINPUTDir.getAbsolutePath(), picardLocalJar, gatkLocalJar, execThread);

        try {
            Process p = Runtime.getRuntime().exec("rm -rf " + fastqTempDir);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                System.out.println("remove redundant fastq file failed");
                System.exit(2);
            }
        } catch (IOException | InterruptedException ie) {
            ie.printStackTrace();
        }
    }

    /**
     * trans sra format to fastq format for the files under
     * @param sraInputDir the directory sra format files, it should contain two sub-directory "IP" and "INPUT"
     * @param fastqOutputDir the tmp output fastq file directory
     * @return true if trans successfully otherwise false
     */
    private static boolean sraToFastq(String sraInputDir, String fastqOutputDir) {
        boolean execSuccess;

        sra2fastq s2fq = new sra2fastq(sraInputDir, fastqOutputDir);
        execSuccess = s2fq.transFormat();
        s2fq = null;

        return execSuccess;
    }

    private static void mkDir(String dirName) {
        File targetDir = new File(dirName);
        if (!targetDir.exists()) {
            boolean res = targetDir.mkdir();
        }
        targetDir = null;
    }
}
