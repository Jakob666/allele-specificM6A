package SamtoolsPileupSNPCalling;

import ReadsMapping.ReadsMapping ;
import org.apache.log4j.Logger;

import java.io.File;


public class ReadsMapping {

    /**
     * use STAR software mapping reads to the reference genome and output alignment SAM file named "Aligned.out.sam"
     * @param refGenomeFile absolute path of reference genome file
     * @param fastqFiles list of fastq file
     * @param execThread working thread number
     */
    public static void alignment(String refGenomeFile, String gtfFile, File[] fastqFiles, int execThread, Logger logger) {
        String refGenomeDir = new File(refGenomeFile).getParent();
        int readLength = ReadsMapping.getFastqReadLength(fastqFiles[0].getAbsolutePath(), logger);
        ReadsMapping.readsMapping(refGenomeDir, refGenomeFile, gtfFile, fastqFiles, readLength, execThread, logger);
    }
}
