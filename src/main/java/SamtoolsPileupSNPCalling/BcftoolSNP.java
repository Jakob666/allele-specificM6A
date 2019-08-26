package SamtoolsPileupSNPCalling;

import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;

public class BcftoolSNP {

    private static BcftoolSNP bcftoolSNP = null;
    private String referenceGenome, bamFile, outputDir, samtools, bcftools, prefix, bcfFile, vcfGzFile;
    private Logger logger;

    private BcftoolSNP(String refGenomeFile, String bamFile, String outputDir, String samtools, String bcftools, Logger logger) {
        this.referenceGenome = refGenomeFile;
        this.bamFile = bamFile;
        this.outputDir = outputDir;
        this.samtools = samtools;
        this.bcftools = bcftools;
        String bamFileName = new File(this.bamFile).getName();
        this.prefix = bamFileName.substring(0, bamFileName.lastIndexOf("_"));
        this.logger = logger;
    }

    public static BcftoolSNP createCaller(String refGenomeFile, String bamFile, String outputDir,
                                          String samtools, String bcftools, Logger logger) {
        if (BcftoolSNP.bcftoolSNP != null)
            BcftoolSNP.bcftoolSNP = null;
        BcftoolSNP.bcftoolSNP = new BcftoolSNP(refGenomeFile, bamFile, outputDir, samtools, bcftools, logger);

        return BcftoolSNP.bcftoolSNP;
    }

    public void mpileUp() {
        this.bcfFile = new File(outputDir, this.prefix + "_mpileup.bcf").getAbsolutePath();
        String cmd = this.samtools + " mpileup -go " + this.bcfFile + " -f " + this.referenceGenome + " -t DP -t SP " + this.bamFile;
        this.logger.debug(cmd);
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                this.logger.error("pileup failed");
                System.exit(2);
            }
        } catch (IOException | InterruptedException ie) {
            this.logger.error("pileup failed\n" + ie.getMessage());
            System.exit(2);
        }
        this.logger.debug("pile up reads complete");
    }

    public void callSnp() {
        this.vcfGzFile = new File(outputDir, this.prefix + "_raw.vcf.gz").getAbsolutePath();
        String cmd = this.bcftools + " call -vmO z -o " + this.vcfGzFile + " " + this.bcfFile;
        this.logger.debug(cmd);
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                this.logger.error("call SNP failed");
                System.exit(2);
            }
        } catch (IOException | InterruptedException ie) {
            this.logger.error("call SNP failed\n" + ie.getMessage());
            System.exit(2);
        }
        cmd = "rm -f " + this.bcfFile;
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0)
                this.logger.error("can not remove file: " + this.bcfFile);
        } catch (IOException | InterruptedException ie) {
            this.logger.error(ie.getMessage());
        }
        this.logger.debug("SNP calling complete");
    }

    public String filterSnp() {
        String outputVcfFile = new File(outputDir, this.prefix + "_filtered.vcf").getAbsolutePath();
        String cmd = String.join(" ", new String[]{this.bcftools, "filter -O v -o ", outputVcfFile, "--SnpGap 30", this.vcfGzFile});
        this.logger.debug(cmd);
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                this.logger.error("filter SNP failed");
                System.exit(2);
            }
        } catch (IOException | InterruptedException ie) {
            this.logger.error("filter SNP failed\n" + ie.getMessage());
            System.exit(2);
        }
        cmd = "rm -f " + this.vcfGzFile;
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0)
                this.logger.error("can not remove file: " + this.vcfGzFile);
        } catch (IOException | InterruptedException ie) {
            this.logger.error(ie.getMessage());
        }
        this.logger.debug("filter SNP complete");
        this.logger.debug("output VCF file " + outputVcfFile);

        return outputVcfFile;
    }
}
