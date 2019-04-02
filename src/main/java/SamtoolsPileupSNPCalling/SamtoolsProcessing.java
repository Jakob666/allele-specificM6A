package SamtoolsPileupSNPCalling;

import GatkSNPCalling.SNPCalling;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;

public class SamtoolsProcessing {

    /**
     * sort, mark duplicates and split reads for the raw output sam file
     * @param refGenomeFile reference genome file
     * @param alignmentResultFile alignment sam file
     * @param outputDir output directory
     * @param prefix output file prefix
     * @param samtools executable samtools file
     * @param picard executable picard jar
     * @param gatk executable gatk jar
     * @param log Logger instance
     * @return pre-process procedure output file
     */
    public static String samFileProcess(String refGenomeFile, String alignmentResultFile, String outputDir, String prefix,
                                        String samtools, String picard, String gatk, Logger log) {
        log.debug("pre-processing alignment result");
        String bamFile = SamtoolsProcessing.sam2bam(alignmentResultFile, outputDir, prefix ,samtools, log);
        String sortedBamFile = SamtoolsProcessing.sorted(samtools, bamFile, log);
        String dedupBamFile = SamtoolsProcessing.deduplicate(picard, sortedBamFile, log);
        String splitFile = splitNCigar(samtools, picard, gatk, refGenomeFile, dedupBamFile, log);

        return splitFile;
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
     * mark duplicates reads in bam file with picard
     * @param picard picard executive file path
     */
    private static String deduplicate(String picard, String sortedFile, Logger log) {
        String dedupFile = new File(sortedFile.substring(0, sortedFile.lastIndexOf("_"))+"_dedup.bam").getAbsolutePath();
        SNPCalling.dropDuplicateReads(picard, sortedFile, dedupFile, log);

        return dedupFile;
    }

    /**
     * split N cigar with gatk tools to reduce false positive rate
     * @param samtools executable samtools file
     * @param picard executable picard jar
     * @param gatk executable gatk jar
     * @param refGenomeFile reference genome file
     * @param dedupFile deduplicated file generate by picard
     * @param log Logger instance
     * @return split N cigar file output by gatk
     */
    private static String splitNCigar(String samtools, String picard, String gatk, String refGenomeFile, String dedupFile, Logger log) {
        String splitFile = new File(dedupFile.substring(0, dedupFile.lastIndexOf("_"))+"_split.bam").getAbsolutePath();
        SNPCalling.refGenomeDict(picard, refGenomeFile, log);
        SNPCalling.createFastaiFile(samtools, refGenomeFile, log);
        SNPCalling.readsTrimReassign(gatk, refGenomeFile, dedupFile, splitFile, log);

        return splitFile;
    }
}
