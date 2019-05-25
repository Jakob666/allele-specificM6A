package SamtoolsPileupSNPCalling;


import org.apache.log4j.Logger;

import java.io.*;
import java.util.HashMap;

/**
 * SamTools SNP calling deprecated
 */
public class AseInference {

    /**
     * get raw SNP reads count with samtools mpileup function
     * @param refGenomeFilePath reference genome file
     * @param sortedBamFile alignment result bam file
     * @param samtools samtools executive file
     * @return abundant file name which records raw SNP sites and reads count
     */
    public static String inferenceASE(String refGenomeFilePath, String sortedBamFile, String samtools, Logger log) {
        log.debug("analysis ASE reads abundant");
        pileup(refGenomeFilePath, sortedBamFile, samtools, log);
        String pileFile = new File(sortedBamFile.substring(0, sortedBamFile.lastIndexOf("_"))+"_pileup.txt").getAbsolutePath();
        String abundantFile = alleleAbundant(pileFile, log);

        return abundantFile;
    }

    /**
     * ASE inference can be carried out using the mpileup function available with SAMtools. The mpileup function requires a
     * sorted BAM file (created in the previous step) and a reference genome in FASTA format.
     * @param refGenomeFilePath reference fasta file absolute path
     * @param sortedBamFile sorted bam file absolute path
     * @param samtools executive samtools file
     */
    private static void pileup(String refGenomeFilePath, String sortedBamFile, String samtools, Logger log) {
        String outputText = new File(sortedBamFile.substring(0, sortedBamFile.lastIndexOf("_"))+"_pileup.txt").getAbsolutePath();
        // samtools mpileup -o output.txt -f /data1/hubs/reference_genome/hg38.fa /data1/hubs/samtoolsTest/alignment_sorted.bam
        String cmd = samtools + " mpileup -o " + outputText + " -f " + refGenomeFilePath + " " + sortedBamFile;
        log.debug("pile up reads, output: " + outputText);

        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                log.debug("pileup failed");
                System.exit(2);
            }
        } catch (IOException | InterruptedException ie) {
            log.debug("pileup failed\n" + ie.getMessage());
            System.exit(2);
        }
    }


    /**
     * Each symbol in the fifth column represents the nucleotide from one read. The fifth column is encoded relative to
     * the reference nucleotide
     * @param pileupFilePath pileupFile path
     */
    private static String alleleAbundant(String pileupFilePath, Logger log) {

        String prefix = pileupFilePath.substring(0, pileupFilePath.lastIndexOf("_"));
        File abundantFile = new File(prefix + "_abundant.txt");
        log.debug("calculate ASE reads abundant, output: " + abundantFile);

        // chr6 410512 T 25 .,,,,,,.,.,.,,,..,......^S. ""!#&%%%%%"&%!"!$%#%%!!"
        try {
            BufferedReader pileupFile = new BufferedReader(
                    new InputStreamReader(new FileInputStream(new File(pileupFilePath)))
            );

            BufferedWriter outputAbundantFile = new BufferedWriter(
                    new OutputStreamWriter(new FileOutputStream(abundantFile))
            );

            String line = "";
            HashMap<String, Integer> col5Info;
            String[] usefulInfo;
            String outputLine;
            int total, ambiguous, maxRead;
            int[] nuclearRead;
            while (line != null) {
                line = pileupFile.readLine();
                if (line != null) {
                    total = 0;
                    ambiguous = 0;
                    String[] colInfo = line.split("\t");
                    String refNuclear = colInfo[2];
                    col5Info = decodeCol5(colInfo[4]);

                    int match = col5Info.getOrDefault("matched", 0);
                    ambiguous = col5Info.getOrDefault("ambiguous", 0);
                    total = col5Info.getOrDefault("A", 0) + col5Info.getOrDefault("T", 0) + col5Info.getOrDefault("C", 0) + col5Info.getOrDefault("G", 0) + match + ambiguous;
                    nuclearRead = new int[]{col5Info.getOrDefault("A", 0), col5Info.getOrDefault("T", 0), col5Info.getOrDefault("C", 0), col5Info.getOrDefault("G", 0)};
                    maxRead = getMaxRead(nuclearRead);

                     col5Info.put(refNuclear, col5Info.getOrDefault(refNuclear, 0) + match);

                    if (0.15 * total <= maxRead && 0.85 * total >= maxRead && total >= 20) {
                        usefulInfo = new String[]{colInfo[0], colInfo[1], colInfo[2], colInfo[3], col5Info.getOrDefault("A", 0).toString(), col5Info.getOrDefault("C", 0).toString(),
                                col5Info.getOrDefault("T", 0).toString(), col5Info.getOrDefault("G", 0).toString(), Integer.toString(ambiguous)};

                        StringBuffer sb = new StringBuffer();

                        for (String str : usefulInfo) {
                            sb.append(str).append("\t");
                        }
                        outputLine = sb.toString();

                        outputAbundantFile.write(outputLine);
                        outputAbundantFile.newLine();
                        sb = null;
                        outputLine = null;
                    }
                }
            }
            pileupFile.close();
            outputAbundantFile.flush();
            outputAbundantFile.close();

            Process p = Runtime.getRuntime().exec("rm -f " + pileupFilePath);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                log.error("rm pileup.txt failed");
            }

        } catch (FileNotFoundException fne) {
            log.debug("reads abundant failed\n" + fne.getMessage());
            System.exit(2);
        } catch (IOException | InterruptedException ie) {
            log.debug("reads abundant failed\n" + ie.getMessage());
            System.exit(3);
        }

        return abundantFile.getAbsolutePath();
    }


    /**
     * given the string of the fifth column returns match, mismatch, ambiguous
     * @param col5Text text of the fifth column
     */
    private static HashMap<String, Integer> decodeCol5(String col5Text) {
        char[] columnText = col5Text.toCharArray();
        if (columnText.length == 0)
            throw new RuntimeException("invalid input");
        HashMap<String, Integer> statisResult = new HashMap<String, Integer>();

        for (char c: columnText) {
            switch (c) {
                case '.':
                case ',':
                case '>':
                case '<':
                    statisResult.put("matched", statisResult.getOrDefault("matched", 0) + 1);
                    break;
                case 'A':
                case 'a':
                    statisResult.put("A", statisResult.getOrDefault("A", 0) + 1);
                    break;
                case 'T':
                case 't':
                    statisResult.put("T", statisResult.getOrDefault("T", 0) + 1);
                    break;
                case 'C':
                case 'c':
                    statisResult.put("C", statisResult.getOrDefault("C", 0) + 1);
                    break;
                case 'G':
                case 'g':
                    statisResult.put("G", statisResult.getOrDefault("G", 0) + 1);
                    break;
                case 'N':
                case 'n':
                    statisResult.put("ambiguous", statisResult.getOrDefault("ambiguous", 0) + 1);
                    break;
            }
        }

        return statisResult;
    }


    /**
     * get max read of the 4 kind of nucleotide
     * @param nuclearReads reads
     * @return max reads number
     */
    private static int getMaxRead(int[] nuclearReads) {
        int max = nuclearReads[0];

        for (int i : nuclearReads) {
            if (i > max) {
                max = i;
            }
        }

        return max;
    }
}
