package heterozygoteSiteAnalysis;

import java.io.*;

public class MergeVCFFile {
    private File VCFDir;
    private File outputFile;

    /**
     * Constructor
     * @param VCFDirectory directory which stores vcf files
     * @param outputVCFFile output vcf file
     */
    public MergeVCFFile(String VCFDirectory, String outputVCFFile) {
        this.VCFDir = new File(VCFDirectory);
        if (this.VCFDir.isFile()) {
            throw new RuntimeException("invalid input vcf directory");
        }
        this.outputFile = new File(outputVCFFile);
    }

    /**
     * merge vcf data to a common file
     * @throws IOException
     */
    public void mergeVCF() throws IOException{
        File[] dirFiles = this.VCFDir.listFiles();
        for (File f : dirFiles) {
            if (f.getName().endsWith("vcf")) {
                this.write2output(f);
            }else
                continue;
        }
    }

    /**
     * write vcf from original file to common output
     * @param vcfFile vcf file File object
     * @throws IOException
     */
    private void write2output(File vcfFile) throws IOException {

        String line;
        try {
            FileReader fr = new FileReader(vcfFile);
            BufferedReader bfr = new BufferedReader(fr);

            FileWriter fw = new FileWriter(this.outputFile, true);
            BufferedWriter bfw = new BufferedWriter(fw);
            while ((line = bfr.readLine()) != null) {
                if (!line.startsWith("#")) {
                    bfw.write(line);
                    bfw.newLine();
                }
            }

            bfw.flush();
            bfw.close();
            bfr.close();

        } catch (FileNotFoundException fe) {
            fe.printStackTrace();
        }
    }

}
