package SamtoolsPileupSNPCalling;

import GatkSNPCalling.sra2fastq;

import java.io.File;

public class SamtoolsPileupSNPCalling {
    public static void main(String[] args) {

        String sraDataDir = args[0];
        String fastqTempDir = args[1];
        String alignmentOutputDir = args[2];
        String genomeFile = args[3];
        String gtfFileDir = args[4];
        String samtools = args[5];
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

        snpCalling(genomeFile, fastqIPDir, alignIPDir.getAbsolutePath(), gtfFileDir, samtools, execThread);
        snpCalling(genomeFile, fastqINPUTDir, alignINPUTDir.getAbsolutePath(), gtfFileDir, samtools, execThread);

    }

    private static void snpCalling(String genomeFilePath, File fastqDir, String outputDir, String gtfDir, String samtools, int execThread) {
        File[] fastqFiles = fastqDir.listFiles();
        if (fastqFiles == null) {
            System.out.println("empty fastq file Dir");
            System.exit(2);
        }

        for (File fq : fastqFiles) {
            String fileName = fq.getName();
            String prefix = fileName.substring(0, fileName.lastIndexOf("."));
            ReadsMapping.alignment(genomeFilePath, fq.getAbsolutePath(), execThread);
            String refGenomeDir = new File(genomeFilePath).getParent();
            String aligmentResultFile = new File(refGenomeDir, "Aligned.out.sam").getAbsolutePath();
            String sortedBamFile = SamtoolsProcessing.samFileProcess(aligmentResultFile, outputDir, prefix, samtools);
            String pileupFile = AseInference.inferenceASE(genomeFilePath, sortedBamFile, samtools);
            SnpFilter sf = new SnpFilter(gtfDir, pileupFile);
            sf.filterVcf();
        }
    }

    private static void mkDir(String dirName) {
        File targetDir = new File(dirName);
        if (!targetDir.exists()) {
            boolean res = targetDir.mkdir();
        }
        targetDir = null;
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
}
