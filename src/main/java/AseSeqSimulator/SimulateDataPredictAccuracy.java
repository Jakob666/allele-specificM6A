package AseSeqSimulator;

import java.io.*;
import java.util.HashMap;
import java.util.HashSet;

/**
 * 将检测结果的peak定位到相应的基因
 */
public class PeakCoveredSNPLocateGene {
    private String significantTestResult, peakCoveredSNPFile, mutationFile, rpkmFile;
    private HashMap<String, HashSet<int[]>> significantChrPeaks = new HashMap<>();
    private HashMap<String, HashSet<Integer>> significantPeakCoveredSNPs = new HashMap<>();
    private HashSet<String> significantPeakCoveredGeneId = new HashSet<>();
    private HashMap<String, Double> significantPeakCoveredGeneAse = new HashMap<>();
    private double accuracy = 0;

    /**
     * Constructor
     * @param significantTestResult 显著性检验的结果文件
     * @param peakCoveredSNPFile peak covered SNP记录文件
     * @param mutationFile 模拟数据的突变位点文件
     */
    public PeakCoveredSNPLocateGene(String significantTestResult, String peakCoveredSNPFile, String mutationFile,
                                    String rpkmFile) {
        this.significantTestResult = significantTestResult;
        this.peakCoveredSNPFile = peakCoveredSNPFile;
        this.mutationFile = mutationFile;
        this.rpkmFile = rpkmFile;
    }

    /**
     * 获取显著性检验结果中每条染色体上ASE specific m6A signal的peak范围
     */
    private void parseSignificantTestResult() {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.significantTestResult))));
            String line = "";
            String[] info;
            String chrNum;
            double qVal;
            int peakStart, peakEnd;
            HashSet<int[]> peaks;
            while (line!=null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    qVal = Double.parseDouble(info[3]);
                    if (qVal-0.05 < 0.000001)
                        break;
                    chrNum = info[0];
                    peakStart = Integer.parseInt(info[1]);
                    peakEnd = Integer.parseInt(info[2]);
                    peaks = this.significantChrPeaks.getOrDefault(chrNum, new HashSet<>());
                    peaks.add(new int[]{peakStart, peakEnd});
                    this.significantChrPeaks.put(chrNum, peaks);
                }
            }
            bfr.close();
        } catch (IOException ie) {
            ie.printStackTrace();
        } finally {
            if (bfr != null) {
                try {
                    bfr.close();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
    }

    /**
     * 获取显著性的peak在每条染色体上覆盖的SNP位点的位置
     */
    private void parsePeakCoveredSNPFile() {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.peakCoveredSNPFile))));
            String line = "", chrNum;
            String[] info;
            HashSet<int[]> significantPeaks;
            HashSet<Integer> peakCoveredSNP;
            int snpPosition, peakStart, peakEnd;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    chrNum = info[0];
                    snpPosition = Integer.parseInt(info[2]);
                    peakStart = Integer.parseInt(info[3]);
                    peakEnd = Integer.parseInt(info[4]);

                    significantPeaks = this.significantChrPeaks.get(chrNum);
                    if (significantPeaks.contains(new int[]{peakStart, peakEnd})) {
                        peakCoveredSNP = this.significantPeakCoveredSNPs.getOrDefault(chrNum, new HashSet<>());
                        peakCoveredSNP.add(snpPosition);
                        this.significantPeakCoveredSNPs.put(chrNum, peakCoveredSNP);
                    }
                }
            }
            bfr.close();
        } catch (IOException ie) {
            ie.printStackTrace();
        } finally {
            if (bfr != null) {
                try {
                    bfr.close();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
    }

    /**
     * 获取显著性peak覆盖的SNP位点所在基因的GeneID
     */
    private void getGeneId() {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.mutationFile))));
            String line = "", chrNum, geneId;
            int lineNum = 1, genomePosition;
            String[] info;
            HashSet<Integer> peakCoveredSnp;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (lineNum == 1) {
                        lineNum++;
                        continue;
                    }
                    info = line.split("\t");
                    chrNum = info[0];
                    geneId = info[1];
                    genomePosition = Integer.parseInt(info[2]);

                    peakCoveredSnp = this.significantPeakCoveredSNPs.getOrDefault(chrNum, null);
                    if (peakCoveredSnp != null && peakCoveredSnp.contains(genomePosition)) {
                        this.significantPeakCoveredGeneId.add(geneId);
                    }
                }
            }
            bfr.close();
        } catch (IOException ie) {
            ie.printStackTrace();
        } finally {
            if (bfr != null) {
                try {
                    bfr.close();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
    }

    /**
     * 获取显著性peak覆盖的SNP所在基因的ASE ratio
     */
    private void getGeneAseRatio() {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.rpkmFile))));
            String line = "", geneID, geneName;
            String[] info;
            int lineNum = 1;
            double majorAlleleRatio;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (lineNum == 1) {
                        lineNum++;
                        continue;
                    }
                    info = line.split("\t");
                    geneID = info[1];
                    geneName = info[2];
                    majorAlleleRatio = Double.parseDouble(info[6]);
                    if (this.significantPeakCoveredGeneId.contains(geneID)) {
                        this.significantPeakCoveredGeneAse.put(geneName, majorAlleleRatio);
                    }
                }
            }
        } catch (IOException ie) {
            ie.printStackTrace();
        } finally {
            if (bfr != null) {
                try {
                    bfr.close();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
    }

    private void getAccuracy() {
        double geneNum = this.significantPeakCoveredGeneAse.size();
        for (Double ase: this.significantPeakCoveredGeneAse.values()) {
            if (ase > 0.5 && ase < 1.0)
                this.accuracy++;
        }
        this.accuracy = this.accuracy / geneNum;
    }
}
