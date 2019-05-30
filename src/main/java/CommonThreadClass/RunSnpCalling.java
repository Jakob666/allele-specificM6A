package CommonThreadClass;

import org.apache.log4j.Logger;
import SamtoolsPileupSNPCalling.SamtoolsPileupSNPCalling;

import java.io.File;

public class RunSnpCalling implements Runnable {
    private String refGenomeFile, inputBamFile, samtools, bcftools;
    private Logger logger;
    public RunSnpCalling(String refGenomeFile, String inputBamFile, String samtools, String bcftools, Logger logger) {
        this.refGenomeFile = refGenomeFile;
        this.inputBamFile = inputBamFile;
        this.samtools = samtools;
        this.bcftools = bcftools;
        this.logger = logger;
    }
    public void run() {
        snpReadsCount();
    }

    /**
     * snp calling and filtration
     */
    private void snpReadsCount() {
        logger.debug("start SNP calling");
        String outputDir = new File(new File(inputBamFile).getParent()).getParent();
        SamtoolsPileupSNPCalling.snpCalling(refGenomeFile, outputDir, inputBamFile, samtools, bcftools, logger);
        logger.debug("SNP calling complete");
    }
}
