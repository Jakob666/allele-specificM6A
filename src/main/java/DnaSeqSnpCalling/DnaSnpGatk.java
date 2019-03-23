package DnaSeqSnpCalling;

import GatkSNPCalling.GatkSNPCalling;

import java.io.File;

public class DnaSnpGatk {

    private String refGenomeFile, bamFile, outputDir, gatk;

    public DnaSnpGatk(String refGenomeFile, String bamFile, String samtool, String bcftool, String gatk, int execThread) {
        this.refGenomeFile = refGenomeFile;
        this.bamFile = bamFile;
        this.gatk = gatk;
        this.outputDir = new File(bamFile).getParent();
        DnaSnpSamtools dss = new DnaSnpSamtools(refGenomeFile, bamFile, samtool, bcftool, execThread);
    }

    private void snpCalling() {
        String prefix = this.bamFile.substring(0, this.bamFile.lastIndexOf("."));
        GatkSNPCalling.variantCalling(this.gatk, this.refGenomeFile, this.outputDir, prefix);
        GatkSNPCalling.variantFilter(this.gatk, this.refGenomeFile, this.outputDir, prefix);
    }
}
