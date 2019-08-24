package CommonThreadClass;

import ReadsMapping.ReadsMapping;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

/**
 * 将DNA-seq的reads比对到参考基因组
 */
public class DnaSeqReadsMapping implements Runnable {
    private File[] fastqFiles;
    private String refGenomeFile, outputDir, prefix, bwa, samtools, picard, gatk;
    private int execthread;
    private Logger logger;

    public DnaSeqReadsMapping(File[] fastqFiles, String refGenomeFile, String outputDir, String prefix, String bwa,
                              String samtools, String picard, String gatk, int execthread, Logger logger) {
        this.fastqFiles = fastqFiles;
        this.refGenomeFile = refGenomeFile;
        this.outputDir = outputDir;
        this.prefix = prefix;
        this.bwa = bwa;
        this.samtools = samtools;
        this.picard = picard;
        this.gatk = gatk;
        this.execthread = execthread;
        this.logger = logger;
    }

    public void run() {
        this.mappingDnaSeqReads();
    }

    private void mappingDnaSeqReads() {
        this.genomeIndex();
        this.alignment();
    }

    /**
     * 为参考基因组建立索引
     */
    private void genomeIndex() {
        String cmd = String.join(" ", new String[]{this.bwa, "index", "-a bwtsw",this.refGenomeFile});
        this.logger.debug("index reference genome " + this.refGenomeFile);
        this.logger.debug(cmd);
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                this.logger.error("reference genome index failed");
                System.exit(2);
            }
        } catch (IOException | InterruptedException ie) {
            this.logger.error(ie.getMessage());
            System.exit(2);
        }
    }

    /**
     * reads align to reference genome
     */
    private void alignment() {
        String alignmentSamFile;

        // bwa 对每个fastq文件进行reads mapping
        for (File fastqFile: this.fastqFiles) {
            String fastqName = fastqFile.getAbsolutePath();
            alignmentSamFile = new File(this.outputDir,fastqName.substring(0, fastqName.lastIndexOf(".")) + ".sam").getAbsolutePath();
            String cmd = String.join(" ", new String[]{this.bwa, "mem -x ont2d -t", Integer.toString(this.execthread),
                                     "-M", this.refGenomeFile, fastqFile.getAbsolutePath(), ">", alignmentSamFile});

            this.logger.debug(fastqFile.getAbsolutePath() + " start reads mapping");
            this.logger.debug(cmd);
            try {
                Process p = Runtime.getRuntime().exec(cmd);
                int exitVal = p.waitFor();
                if (exitVal != 0) {
                    this.logger.error("alignment failed");
                }
            } catch (IOException | InterruptedException ie) {
                this.logger.error(ie.getMessage());
                System.exit(2);
            }
            this.logger.debug("alignment complete, output file " + alignmentSamFile);

            String head = alignmentSamFile.substring(0, alignmentSamFile.lastIndexOf("."));
            File sortedBamFile = new File(head + "_sorted.bam");
            File sortedBaiFile = new File(head + "_sorted.bai");
            ReadsMapping.readsGroup(this.picard, alignmentSamFile, sortedBamFile.getAbsolutePath(), this.logger);
            File deduplicatedBamFile = new File(head + "_deduplicated.bam");
            File deduplicatedBaiFile = new File(head + "_deduplicated.bai");
            ReadsMapping.dropDuplicateReads(this.picard, sortedBamFile.getAbsolutePath(), deduplicatedBamFile.getAbsolutePath(), this.logger);
            this.deleteFile(sortedBamFile);
            this.deleteFile(sortedBaiFile);
            File finalBamFile = new File(head + "_alignment.bam");
            File finalBaiFile = new File(head + "_alignment.bai");

            // 文件大于2G目前跳过split过程，因为时间太长
            if (deduplicatedBamFile.length() > new Long("2147483648")) {
                String cmd1 = "mv " + deduplicatedBamFile.getAbsolutePath() + " " + finalBamFile.getAbsolutePath();
                String cmd2 = "mv " + deduplicatedBaiFile.getAbsolutePath() + " " + finalBaiFile.getAbsolutePath();
                try {
                    Process p = Runtime.getRuntime().exec(cmd1);
                    int exitVal = p.waitFor();
                    if (exitVal != 0) {
                        this.logger.error("file rename failed");
                        System.exit(2);
                    }
                    p = Runtime.getRuntime().exec(cmd2);
                    exitVal = p.waitFor();
                    if (exitVal != 0) {
                        this.logger.error("file rename failed");
                        System.exit(2);
                    }
                } catch (IOException | InterruptedException ie) {
                    this.logger.error(ie.getMessage());
                    System.exit(2);
                }
            } else {
                ReadsMapping.refGenomeDict(this.picard, this.refGenomeFile, this.logger);
                ReadsMapping.createFastaiFile(samtools, this.refGenomeFile, this.logger);
                ReadsMapping.readsTrimReassign(this.gatk, this.refGenomeFile, deduplicatedBamFile.getAbsolutePath(), finalBamFile.getAbsolutePath(), this.logger);
                deleteFile(deduplicatedBamFile);
                deleteFile(deduplicatedBaiFile);
            }
        }
        // 将fastq文件的比对结果bam文件合并
        String mergedBamFile = new File(this.outputDir, this.prefix + "_alignment.bam").getAbsolutePath();
        File[] files = new File(this.outputDir).listFiles();
        if (files == null) {
            this.logger.error("empty alignment result directory");
            System.exit(2);
        }
        ArrayList<String> bamFiles = new ArrayList<>();
        for (File file: files) {
            if (file.getAbsolutePath().endsWith("bam"))
                bamFiles.add(file.getAbsolutePath());
        }
        this.mergeBamFiles(bamFiles, mergedBamFile);
    }

    /**
     * merge input and ip sample alignment bam files
     * @param bamFiles bam files list
     */
    private void mergeBamFiles(ArrayList<String> bamFiles, String mergedBamFile) {
        // java -jar picard.jar MergeSamFiles \
        //      I=input_1.bam \
        //      I=input_2.bam \
        //      O=output_merged_files.bam
        String cmd;
        if (bamFiles.size() == 1) {
            try {
                Process p = Runtime.getRuntime().exec("mv " + bamFiles.get(0) + " " + mergedBamFile);
                int exitVal = p.waitFor();
                if (exitVal != 0) {
                    this.logger.error("move failed");
                    System.exit(2);
                }
            } catch (IOException | InterruptedException ie) {
                this.logger.error(ie.getMessage());
                System.exit(2);
            }
        } else {
            StringBuilder sb = new StringBuilder();
            sb.append(this.picard);
            sb.append(" MergeSamFiles");
            for (String bamFile: bamFiles) {
                String input = " I=" + bamFile;
                sb.append(input);
            }
            String output = " O="+mergedBamFile;
            sb.append(output);
            cmd = sb.toString();
            this.logger.debug("merge bam files\n" + cmd);
            try {
                Process p = Runtime.getRuntime().exec(cmd);
                int exitVal = p.waitFor();
                if (exitVal != 0) {
                    this.logger.error("merge bam files failed");
                    System.exit(2);
                }
            } catch (IOException | InterruptedException ie) {
                this.logger.error(ie.getMessage());
            }
        }

        for (String bamFile: bamFiles) {
            cmd = "rm -f " + bamFile;
            try {
                Process p = Runtime.getRuntime().exec(cmd);
                int exitVal = p.waitFor();
                if (exitVal != 0) {
                    this.logger.error("can not remove redundant BAM file");
                    System.exit(2);
                }
            } catch (IOException | InterruptedException ie) {
                this.logger.error(ie.getMessage());
                System.exit(2);
            }
        }
    }

    private void deleteFile(File targetFile) {
        if (targetFile.exists()) {
            boolean res = targetFile.delete();
            if (!res)
                this.logger.error("can not remove redundant SAM file " + targetFile);
        }
    }
}
