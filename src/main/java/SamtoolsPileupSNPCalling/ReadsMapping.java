package SamtoolsPileupSNPCalling;

import GatkSNPCalling.SNPCalling;

import java.io.File;

public class ReadsMapping {

    /**
     * use STAR software mapping reads to the reference genome and output alignment SAM file named "Aligned.out.sam"
     * @param refGenomeFile absolute path of reference genome file
     * @param fastqFile absolute path of fastq file
     * @param execThread working thread number
     */
    public static void alignment(String refGenomeFile, String fastqFile, int execThread) {
        String refGenomeDir = new File(refGenomeFile).getParent();
        int readLength = SNPCalling.getFastqReadLength(fastqFile);
        SNPCalling.readsMapping(refGenomeDir, refGenomeFile, fastqFile, readLength, execThread);
    }
}
