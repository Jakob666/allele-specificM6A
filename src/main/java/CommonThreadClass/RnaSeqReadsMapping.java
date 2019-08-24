package CommonThreadClass;

import org.apache.log4j.Logger;
import ReadsMapping.ReadsMapping;

import java.io.File;

/**
 * 将RNA-seq的reads比对到参考基因组
 */
public class RnaSeqReadsMapping implements Runnable {
    private File[] fastqFiles;
    private String refGenomeFile, gtfFile, picardJar, gatkJar, samtools, outputDir, prefix, direct;
    private int execThread;
    private boolean zip;
    private Logger logger;

    public void run() {
        String alignmentBamFile = this.mappingRnaSeqReads();
        this.logger.debug("alignment result output in " + alignmentBamFile);
    }

    public RnaSeqReadsMapping(File[] fastqFiles, String refGenomeFile, String gtfFile, String picardJar, String gatkJar,
                              String samtools, String outputDir, String prefix, String direct, int execThread, boolean zip,
                              Logger logger) {
        this.fastqFiles = fastqFiles;
        this.refGenomeFile = refGenomeFile;
        this.gtfFile = gtfFile;
        this.picardJar = picardJar;
        this.gatkJar = gatkJar;
        this.samtools = samtools;
        this.outputDir = outputDir;
        this.prefix = prefix;
        this.direct = direct;
        this.execThread = execThread;
        this.zip = zip;
        this.logger = logger;
    }

    private String mappingRnaSeqReads() {
        String refGenomeDir = new File(this.refGenomeFile).getParent();
        int readLength = ReadsMapping.getFastqReadLength(this.fastqFiles[0].getAbsolutePath(), this.zip, this.logger);
        String outFileNamePrefix = new File(this.outputDir, this.prefix).getAbsolutePath();
        ReadsMapping.readsMapping(refGenomeDir, this.refGenomeFile, this.gtfFile, this.fastqFiles, outFileNamePrefix,
                                  this.direct, readLength, this.execThread, this.zip, this.logger);
        // 返回最终生成的bam文件路径
        return ReadsMapping.filterMappedReads(this.refGenomeFile, outFileNamePrefix, this.picardJar, this.gatkJar,
                                              this.samtools, this.logger);
    }
}
