package AseSeqSimulator;


import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;

public class GeneExpDistribution {

    private GeneExpDistribution() {

    }

    public static GeneExpDistribution getInstance() {
        return new GeneExpDistribution();
    }

    /**
     * get gene expression mean value and standard deviation in different cell lines
     * @param cellLineExpFile downloading file from http://medicalgenomics.org/rna_seq_atlas/download
     * @return gene expression data
     */
    public HashMap<String, double[]> experimentalGeneExp(String cellLineExpFile) {
        HashMap<String, double[]> geneExp = new HashMap<>();
        geneExp = this.readFromFile(cellLineExpFile, geneExp);

        return geneExp;
    }

    private HashMap<String, double[]> readFromFile(String cellLineExpFile, HashMap<String, double[]> geneExp) {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(
                    new InputStreamReader(new FileInputStream(new File(cellLineExpFile)))
            );
            String line = "", geneName;
            int lineNum = 0;
            String[] lineInfo;
            ArrayList<Double> expData;

            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (lineNum == 0) {
                        lineNum ++;
                        continue;
                    }
                    expData = new ArrayList<>();
                    double sum = 0, mean, std = 0;
                    lineInfo = line.split("\t");
                    geneName = lineInfo[2];
                    for (int i = 5; i < lineInfo.length; i++) {
                        Double exp = Double.parseDouble(lineInfo[i]);
                        sum = sum + exp;
                        expData.add(exp);
                    }
                    mean = sum / expData.size();
                    for (Double data: expData) {
                        std = std + (data - mean) * (data - mean);
                    }
                    std = Math.sqrt(std / (expData.size() - 1));

                    geneExp.put(geneName, new double[]{mean, std});
                    expData = null;
                }
            }
            bfr.close();
        } catch (IOException ie) {
            ie.printStackTrace();
        } finally {
            if (bfr != null) {
                try {
                    bfr.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }

        return geneExp;
    }
}
