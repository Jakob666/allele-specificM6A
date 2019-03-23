package DnaSeqSnpCalling;

import java.io.File;
import java.io.IOException;

public class DnaSnpSamtools {

    private String refGenomeFile, bamFile, samtool, bcftool;
    private int execThread;

    public DnaSnpSamtools(String refGenomeFile, String bamFile, String samtool, String bcftool, int execThread) {
        this.refGenomeFile = refGenomeFile;
        this.bamFile = bamFile;
        this.samtool = samtool;
        this.bcftool = bcftool;
        this.execThread = execThread;
    }

    /**
     * SNP calling with samtool using DNA alignment sam file
     */
    public void snpCalling() {
        formBcfFile();
        bcf2Vcf();
        snpFilter();
        finalVcf();
    }

    /**
     * generate bcf file from alignment result
     */
    private void formBcfFile() {

        int lastDot = this.bamFile.lastIndexOf(".");
        String bcfFile = new File(this.bamFile.substring(0, lastDot) + "_output.bcf").getAbsolutePath();

        String cmd = this.samtool + " mpileup -go " + bcfFile + " -f " + this.refGenomeFile + " " + this.bamFile;
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                throw new RuntimeException("generate bcf file failed");
            }
        } catch (IOException | InterruptedException ie) {
            ie.printStackTrace();
            System.exit(2);
        }
    }

    /**
     * generate raw vcf file from bcf file
     */
    private void bcf2Vcf() {

        int lastDot = this.bamFile.lastIndexOf(".");
        String bcfFile = new File(this.bamFile.substring(0, lastDot) + "_output.bcf").getAbsolutePath();
        String vcfFile = new File(this.bamFile.substring(0, bcfFile.lastIndexOf(".")) + "_raw.vcf.gz").getAbsolutePath();
        String cmd = this.bcftool + " call -vmO z --threads " + this.execThread + " -o " + vcfFile + " " + bcfFile;

        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                throw new RuntimeException("generate vcf from bcf file failed");
            }
        } catch (IOException | InterruptedException ie) {
            ie.printStackTrace();
            System.exit(2);
        }
    }

    /**
     * filter SNP with several parameters, set low quality threshold and minimize distance between two SNP sites
     */
    private void snpFilter() {

        int lastDot = this.bamFile.lastIndexOf(".");
        String bcfFile = new File(this.bamFile.substring(0, lastDot) + "_output.bcf").getAbsolutePath();
        String rawVcf = new File(this.bamFile.substring(0, bcfFile.lastIndexOf(".")) + "_raw.vcf.gz").getAbsolutePath();
        String filteredVcf = new File(this.bamFile.substring(0, bcfFile.lastIndexOf(".")) + "_filtered.vcf").getAbsolutePath();

        String cmd = this.bcftool + " filter -O v -o " + filteredVcf + " -s LOWQUAL -e 'QUAL<10' --SnpGap 20 --set-GTs . " + rawVcf;

        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                throw new RuntimeException("filter raw vcf failed");
            }
        } catch (IOException | InterruptedException ie) {
            ie.printStackTrace();
            System.exit(2);
        }
    }

    /**
     * select all SNP with high quality and put into the final output file
     */
    private void finalVcf() {

        int lastDot = this.bamFile.lastIndexOf(".");
        String bcfFile = new File(this.bamFile.substring(0, lastDot) + "_output.bcf").getAbsolutePath();
        String filteredVcf = new File(this.bamFile.substring(0, bcfFile.lastIndexOf(".")) + "_filtered.vcf").getAbsolutePath();
        String highQualVcf = new File(this.bamFile.substring(0, bcfFile.lastIndexOf(".")) + "_HQ.vcf").getAbsolutePath();


        String cmd = this.bcftool + " view -v snps " + filteredVcf + " > " + highQualVcf;
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                throw new RuntimeException("filter raw vcf failed");
            }
        } catch (IOException | InterruptedException ie) {
            ie.printStackTrace();
            System.exit(2);
        }
    }

}
