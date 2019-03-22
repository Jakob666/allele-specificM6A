package SamtoolsPileupSNPCalling;

import java.io.File;
import java.io.IOException;

public class SamtoolsProcessing {

    public static void samFileProcess(String alignmentResultFile, String outputDir, String sraNum, String samtools) {
        SamtoolsProcessing.sam2bam(alignmentResultFile, outputDir, sraNum, samtools);

        String bamFilePath = new File(outputDir, sraNum + "_alignment.bam").getAbsolutePath();
        SamtoolsProcessing.sorted(samtools, bamFilePath, sraNum);
    }

    /**
     * transform file format from sam to bam
     * @param alignmentResultFile alignment result file output by STAR
     * @param outputDir bam file directory
     * @param sraNum sra number
     * @param samtools samtools executive file path
     */
    private static void sam2bam(String alignmentResultFile, String outputDir, String sraNum, String samtools) {
        String bamFilePath = new File(outputDir, sraNum + "_alignment.bam").getAbsolutePath();
        String cmd = samtools + " view -bS -o " + bamFilePath + " " + alignmentResultFile;

        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                throw new RuntimeException("trans sam to bam failed");
            }
        } catch (IOException | InterruptedException ie) {
            ie.printStackTrace();
        }
    }

    /**
     * sort the alignment bam file with samtools
     * @param samtools samtools executive file path
     * @param bamFilePath alignment result bam file path
     * @param sraNum sra number
     */
    private static void sorted(String samtools, String bamFilePath, String sraNum) {
        String cmd = samtools + " sort " + bamFilePath + " -o " + sraNum+"_alignment.sort";
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                throw new RuntimeException("samtools sort failed");
            }
        }catch (IOException | InterruptedException ie) {
            ie.printStackTrace();
        }
    }
}
