package CommonThreadClass;

import org.apache.log4j.Logger;
import ReadsMapping.ReadsMapping;

import java.io.File;

/**
 * 将RNA-seq的reads比对到参考基因组
 */
public class RnaSeqReadsMapping implements Runnable {
    private File[] fastqFiles;
    private String refGenomeFile, gtfFile, picardJar, gatkJar, samtools, outputDir, prefix;
    private int execThread;
    private Logger logger;

    public void run() {
        String alignmentBamFile = mappingRnaSeqReads();
        this.logger.debug("alignment result output in " + alignmentBamFile);
    }

    public RnaSeqReadsMapping(File[] fastqFiles, String refGenomeFile, String gtfFile, String picardJar, String gatkJar,
                              String samtools, String outputDir, String prefix, int execThread, Logger logger) {
        this.fastqFiles = fastqFiles;
        this.refGenomeFile = refGenomeFile;
        this.gtfFile = gtfFile;
        this.picardJar = picardJar;
        this.gatkJar = gatkJar;
        this.samtools = samtools;
        this.outputDir = outputDir;
        this.prefix = prefix;
        this.execThread = execThread;
        this.logger = logger;
    }

    private String mappingRnaSeqReads() {
        String refGenomeDir = new File(this.refGenomeFile).getParent();
        int readLength = ReadsMapping.getFastqReadLength(this.fastqFiles[0].getAbsolutePath(), this.logger);
        ReadsMapping.readsMapping(refGenomeDir, this.refGenomeFile, this.gtfFile, this.fastqFiles, readLength,
                                  this.execThread, this.logger);
        // 返回最终生成的bam文件路径
        return ReadsMapping.filterMappedReads(this.refGenomeFile, this.picardJar, this.gatkJar, this.samtools,
                                              this.outputDir, this.prefix, this.logger);
    }
}
