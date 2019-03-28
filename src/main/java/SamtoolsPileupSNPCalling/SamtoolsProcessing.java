package SamtoolsPileupSNPCalling;

import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;

public class SamtoolsProcessing {

    public static String samFileProcess(String alignmentResultFile, String outputDir, String prefix, String samtools, Logger log) {
        log.debug("pre-processing alignment result");
        String bamFile = SamtoolsProcessing.sam2bam(alignmentResultFile, outputDir, prefix ,samtools, log);
        String sortedBamFile = SamtoolsProcessing.sorted(samtools, bamFile, log);
        String dedupBamFile = SamtoolsProcessing.deduplicate(samtools, sortedBamFile, log);

        return dedupBamFile;
    }

    /**
     * transform file format from sam to bam
     * @param alignmentResultFile alignment result file output by STAR
     * @param samtools samtools executive file path
     */
    private static String sam2bam(String alignmentResultFile, String outputDir, String prefix, String samtools, Logger log) {
        String bamFilePath = new File(outputDir, prefix + "_alignment.bam").getAbsolutePath();
        String cmd = samtools + " view -bS -o " + bamFilePath + " " + alignmentResultFile;
        log.debug("transform sam file into bam format, output: " + bamFilePath);

        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                log.error("trans sam to bam failed");
                System.exit(2);
            }
        } catch (IOException | InterruptedException ie) {
            log.error("trans sam to bam failed\n" + ie.getMessage());
            System.exit(2);
        }

        return bamFilePath;
    }

    /**
     * sort the alignment bam file with samtools
     * @param samtools samtools executive file path
     */
    private static String sorted(String samtools, String bamFile, Logger log) {
        String sortedBamFile = new File(bamFile.substring(0, bamFile.lastIndexOf("_"))+"_sort.bam").getAbsolutePath();

        String cmd = samtools + " sort " + bamFile + " -o " + sortedBamFile;
        log.debug("sort bam file, output: " + sortedBamFile);
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                log.error("samtools sort failed");
                System.exit(2);
            }
            p = Runtime.getRuntime().exec("rm -f " + bamFile);
            exitVal = p.waitFor();
            if (exitVal != 0) {
                log.error("can not delete bam file: " + bamFile);
            }
        }catch (IOException | InterruptedException ie) {
            log.error("samtools sort failed\n" + ie.getMessage());
            System.exit(2);
        }

        return sortedBamFile;
    }

    /**
     * mark duplicates reads in bam file with samtools
     * @param samtools samtools executive file path
     */
    private static String deduplicate(String samtools, String sortedFile, Logger log) {
        String dedupFile = new File(sortedFile.substring(0, sortedFile.lastIndexOf("_"))+"_dedup.bam").getAbsolutePath();

        String cmd = samtools + " markdup " + sortedFile + " " + dedupFile;
        log.debug("mark duplicated reads in bam file, output: " + dedupFile);
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                log.error("samtools deduplicate failed");
                System.exit(2);
            }
            p = Runtime.getRuntime().exec("rm -f " + sortedFile);
            exitVal = p.waitFor();
            if (exitVal != 0) {
                log.error("can not delete bam file: " + sortedFile);
            }
        }catch (IOException | InterruptedException ie) {
            log.error("samtools deduplicate failed\n" + ie.getMessage());
            System.exit(2);
        }

        return dedupFile;
    }
}
