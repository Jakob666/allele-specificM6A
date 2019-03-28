package SamtoolsPileupSNPCalling;

import java.io.File;
import java.io.IOException;

public class SamtoolsProcessing {

    public static String samFileProcess(String alignmentResultFile, String outputDir, String prefix, String samtools) {
        String bamFile = SamtoolsProcessing.sam2bam(alignmentResultFile, outputDir, prefix ,samtools);
        String sortedBamFile = SamtoolsProcessing.sorted(samtools, bamFile);
        String dedupBamFile = SamtoolsProcessing.deduplicate(samtools, sortedBamFile);

        return dedupBamFile;
    }

    /**
     * transform file format from sam to bam
     * @param alignmentResultFile alignment result file output by STAR
     * @param samtools samtools executive file path
     */
    private static String sam2bam(String alignmentResultFile, String outputDir, String prefix, String samtools) {
        String bamFilePath = new File(outputDir, prefix + "_alignment.bam").getAbsolutePath();
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

        return bamFilePath;
    }

    /**
     * sort the alignment bam file with samtools
     * @param samtools samtools executive file path
     */
    private static String sorted(String samtools, String bamFile) {
        String sortedBamFile = new File(bamFile.substring(0, bamFile.lastIndexOf("_"))+"_sort.bam").getAbsolutePath();

        String cmd = samtools + " sort " + bamFile + " -o " + sortedBamFile;
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                throw new RuntimeException("samtools sort failed");
            }
        }catch (IOException | InterruptedException ie) {
            ie.printStackTrace();
        }

        return sortedBamFile;
    }

    /**
     * mark duplicates reads in bam file with samtools
     * @param samtools samtools executive file path
     */
    private static String deduplicate(String samtools, String sortedFile) {
        String dedupFile = new File(sortedFile.substring(0, sortedFile.lastIndexOf("_"))+"_dedup.bam").getAbsolutePath();

        String cmd = samtools + " markdup " + sortedFile + " " + dedupFile;
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                throw new RuntimeException("samtools deduplicate failed");
            }
        }catch (IOException | InterruptedException ie) {
            ie.printStackTrace();
        }

        return dedupFile;
    }
}
