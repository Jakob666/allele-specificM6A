package DnaSeqSnpCalling;

import java.io.File;
import java.io.IOException;

public class ReadsMapping {

    private String refGenomeFile, fastqFile, refGenomeIndexFile, alignmentFile, bamFilePath, samtools, picard;
    private int execthread;

    /**
     * Constructor
     * @param refGenomeFile reference genome file path
     * @param fastqFile fastq sequencing result
     * @param samtools samtool executive file
     * @param picard picard jar package
     * @param execthread working threads
     */
    public ReadsMapping(String refGenomeFile, String fastqFile, String samtools, String picard, int execthread) {
        this.refGenomeFile = refGenomeFile;
        this.fastqFile = fastqFile;
        this.execthread = execthread;
        String referenceDir = new File(refGenomeFile).getParent();
        this.refGenomeIndexFile = new File(referenceDir, "refGenomeidx").getAbsolutePath();
        this.alignmentFile = new File(referenceDir, "alignment.sam").getAbsolutePath();
        this.samtools = samtools;
        this.picard = picard;
    }

    /**
     * DNA-Seq reads mapping workflow
     */
    public void dnaReadsMapping() {
        genomeIndex();
        alignment();
        sam2bam();
        sortBam();
        rmDuplicate();
    }

    /**
     * generate genome index by bowtie
     */
    private void genomeIndex() {
        String cmd = "bowtie2-build " + this.refGenomeFile + " " + this.refGenomeIndexFile + " --threads " + this.execthread;
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                throw new RuntimeException("reference genome index failed");
            }
        } catch (IOException | InterruptedException ie) {
            ie.printStackTrace();
        }
    }

    /**
     * reads align to reference genome
     */
    private void alignment() {
        // bowtie2 -x hg38idx -U /data/hbs/SRR847362.fastq -S alignment.sam
        String cmd = "bowtie2 -x " + this.refGenomeIndexFile + " -U " + this.fastqFile + " -S " + this.alignmentFile;
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                throw new RuntimeException("alignment failed");
            }
        } catch (IOException | InterruptedException ie) {
            ie.printStackTrace();
        }
    }

    /**
     * trans alignment output result in sam format into bam format
     */
    private void sam2bam() {

        String outputDir = new File(this.alignmentFile).getParent();
        this.bamFilePath = new File(outputDir, "alignment.bam").getAbsolutePath();

        String cmd = this.samtools + " view -bS -o " + this.bamFilePath + " --threads " + this.execthread + " " + this.alignmentFile;
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                throw new RuntimeException("sam file transform failed");
            }
        } catch (IOException | InterruptedException ie) {
            ie.printStackTrace();
            System.exit(2);
        }
    }

    /**
     * sort bam file
     */
    private void sortBam() {

        String sortedFileName = this.bamFilePath.substring(0, this.bamFilePath.lastIndexOf(".")) + "_sorted.bam";
        String cmd = this.samtools + " sort " + this.bamFilePath + " --threads " + this.execthread + " -o " + sortedFileName;
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                throw new RuntimeException("sort bam file failed");
            }
        } catch (IOException | InterruptedException ie) {
            ie.printStackTrace();
            System.exit(2);
        }
    }

    /**
     * mark duplicate reads with picard tool which is said better than samtools
     */
    private void rmDuplicate() {

        String sortedFileName = this.bamFilePath.substring(0, this.bamFilePath.lastIndexOf(".")) + "_sorted.bam";
        String dedupFileName = this.bamFilePath.substring(0, this.bamFilePath.lastIndexOf(".")) + "_dedup.bam";
        String cmd = "java -jar " + this.picard + " MarkDuplicates I=" + sortedFileName + " O=" + dedupFileName +
                " CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics";
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                throw new RuntimeException("sort bam file failed");
            }
        } catch (IOException | InterruptedException ie) {
            ie.printStackTrace();
            System.exit(2);
        }
    }

}
