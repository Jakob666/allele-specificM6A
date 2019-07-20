package AseSeqSimulator;

import java.io.*;
import java.util.HashMap;
import java.util.HashSet;

/**
 * 将检测结果的peak定位到相应的基因，检验流程的预测准确率
 */
public class SimulateDataPredictAccuracy {
    private String significantAseTestResult, significantAsmTestResult, peakCoveredSNPFile, mutationFile, rpkmFile, peakBedFile, simulateM6aFile,
            snpCallingFile;
    private int fragmentLength;
    public HashMap<String, HashSet<String>> significantAseChrPeaks = new HashMap<>(), nonsignificantAseChrPeaks = new HashMap<>(),
            significantAsmChrPeaks = new HashMap<>(), nonsignificantAsmChrPeaks = new HashMap<>();
    private HashMap<String, HashSet<Integer>> predSignificantPeakCoveredSNP = new HashMap<>(), predNonSignificantPeakCoveredSNP = new HashMap<>();
    private HashMap<String, HashSet<Integer>> simulatePeakCoveredSNP = new HashMap<>();
    public HashMap<String, HashSet<Double>> callingPeakCoveredM6a = new HashMap<>();
    private HashSet<String> significantAsePeakCoveredGeneId = new HashSet<>(), nonsignificantAsePeakCoveredGeneId = new HashSet<>();
    private HashMap<String, Double> significantPeakCoveredGeneAse = new HashMap<>(), nonsignificantPeakCoveredGeneAse = new HashMap<>();
    private HashMap<String, HashSet<int[]>> peakCallingIntervals = new HashMap<>(), simulatePeakIntervals = new HashMap<>();
    private HashMap<String, HashMap<Integer, Boolean>> simulatedM6aSites = new HashMap<>();
    private HashMap<String, HashSet<Integer>> snpCallingRes = new HashMap<>(), simulatedSnpSites = new HashMap<>();

    /**
     * Constructor
     * @param significantAseTestResult 显著性检验的结果文件
     * @param peakCoveredSNPFile peak covered SNP记录文件
     * @param mutationFile 模拟数据的突变位点文件
     * @param rpkmFile 基因表达值的记录文件
     * @param peakBedFile peak calling得到的文件
     */
    public SimulateDataPredictAccuracy(String significantAseTestResult, String significantAsmTestResult,
                                       String peakCoveredSNPFile, String mutationFile, String rpkmFile, String peakBedFile,
                                       String simulateM6aFile, String snpCallingFile, int fragmentLength) {
        this.significantAseTestResult = significantAseTestResult;
        this.significantAsmTestResult = significantAsmTestResult;
        this.peakCoveredSNPFile = peakCoveredSNPFile;
        this.mutationFile = mutationFile;
        this.rpkmFile = rpkmFile;
        this.peakBedFile = peakBedFile;
        this.simulateM6aFile = simulateM6aFile;
        this.snpCallingFile = snpCallingFile;
        this.fragmentLength = fragmentLength;
    }

    /**
     * 获取snp calling的结果
     * snpCallingRes -> {chrNum: (site1, site2,...)}
     */
    private void parseSnpCallingFile() {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.snpCallingFile))));
            String line = "", chrNum;
            String[] info;
            int position;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    chrNum = info[0];
                    position = Integer.parseInt(info[1]);
                    HashSet<Integer> chrSnpSites = this.snpCallingRes.getOrDefault(chrNum, new HashSet<>());
                    chrSnpSites.add(position);
                    this.snpCallingRes.put(chrNum, chrSnpSites);
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
     * 获取模拟的SNP结果
     * simulatedSnpSites -> {chrNum: (site1, site2,...)}
     */
    private void simulateSnpSites() {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.mutationFile))));
            String line = "", chrNum;
            String[] info;
            int position, lineNum = 1;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (lineNum == 1) {
                        lineNum++;
                        continue;
                    }
                    info = line.split("\t");
                    chrNum = info[0];
                    position = Integer.parseInt(info[2]);
                    HashSet<Integer> chrSnpSites = this.simulatedSnpSites.getOrDefault(chrNum, new HashSet<>());
                    chrSnpSites.add(position);
                    this.simulatedSnpSites.put(chrNum, chrSnpSites);
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
     * 查看SNP calling得到的VCF文件中多少个位点是模拟的SNP位点
     * @return VCF文件中正确预测的位点占总模拟位点数目的比例
     */
    public double[] getSNPCallingResult() {
        this.parseSnpCallingFile();
        this.simulateSnpSites();
        int totalSimulateSnp = 0, totalPredSnp = 0;
        double truePositiveSNP = 0, snpCallingRecall, snpCallingAcc;
        for (String chrNum: this.simulatedSnpSites.keySet()) {
            HashSet<Integer> chrSnpSites = this.simulatedSnpSites.get(chrNum);
            HashSet<Integer> callingSites = this.snpCallingRes.getOrDefault(chrNum, null);
            totalSimulateSnp += chrSnpSites.size();
            if (callingSites == null)
                continue;
            totalPredSnp += callingSites.size();
            chrSnpSites.retainAll(callingSites);
            truePositiveSNP += chrSnpSites.size();
        }
        snpCallingRecall = truePositiveSNP / totalSimulateSnp;
        snpCallingAcc = truePositiveSNP / totalPredSnp;

        return new double[]{snpCallingAcc, snpCallingRecall};
    }


    /**
     * 获取显著性检验结果中每条染色体上ASE specific m6A signal(q<0.05)及非特异(q>=0.05)及非特异(q>)的peak范围
     * significantAseChrPeaks -> {chrNum: (start1->end1, start2->end2,...)}
     */
    private void parseSignificantAseTestResult() {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.significantAseTestResult))));
            String line = "";
            String[] info;
            String chrNum, peakStart, peakEnd;
            double qVal;
            HashSet<String> peaks;
            while (line!=null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    qVal = Double.parseDouble(info[3]);
                    chrNum = info[0];
                    peakStart = info[1];
                    peakEnd = info[2];
                    if (qVal-0.05 > 0.000001) {
                        peaks = this.nonsignificantAseChrPeaks.getOrDefault(chrNum, new HashSet<>());
                        peaks.add(String.join("->", new String[]{peakStart, peakEnd}));
                        this.nonsignificantAseChrPeaks.put(chrNum, peaks);
                    } else {
                        peaks = this.significantAseChrPeaks.getOrDefault(chrNum, new HashSet<>());
                        peaks.add(String.join("->", new String[]{peakStart, peakEnd}));
                        this.significantAseChrPeaks.put(chrNum, peaks);
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
     * 获取ASE显著和非显著的peak在每条染色体上覆盖的SNP位点
     * predSignificantPeakCoveredSNP -> {chrNum: (site1, site2,...)}
     */
    private void parsePeakCoveredSNPFile() {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.peakCoveredSNPFile))));
            String line = "", chrNum, peakStart, peakEnd;
            String[] info;
            HashSet<String> significantPeaks, nonsignificantPeaks;
            HashSet<Integer> peakCoveredSNP;
            int snpPosition;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    chrNum = info[0];
                    snpPosition = Integer.parseInt(info[2]);
                    peakStart = info[3];
                    peakEnd = info[4];

                    significantPeaks = this.significantAseChrPeaks.getOrDefault(chrNum, null);
                    nonsignificantPeaks = this.nonsignificantAseChrPeaks.getOrDefault(chrNum, null);
                    if (significantPeaks == null && nonsignificantPeaks == null)
                        continue;

                    if (significantPeaks != null && significantPeaks.contains(String.join("->", new String[]{peakStart, peakEnd}))) {
                        peakCoveredSNP = this.predSignificantPeakCoveredSNP.getOrDefault(chrNum, new HashSet<>());
                        peakCoveredSNP.add(snpPosition);
                        this.predSignificantPeakCoveredSNP.put(chrNum, peakCoveredSNP);
                    } else if (nonsignificantPeaks != null && nonsignificantPeaks.contains(String.join("->", new String[]{peakStart, peakEnd}))) {
                        peakCoveredSNP = this.predNonSignificantPeakCoveredSNP.getOrDefault(chrNum, new HashSet<>());
                        peakCoveredSNP.add(snpPosition);
                        this.predNonSignificantPeakCoveredSNP.put(chrNum, peakCoveredSNP);
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
     * 获取显著和非显著的peak覆盖的SNP位点所在基因的GeneID
     * significantAsePeakCoveredGeneId -> (gene1, gene2,...)
     */
    private void getAseGeneId() {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.mutationFile))));
            String line = "", chrNum, geneId;
            int lineNum = 1, genomePosition;
            String[] info;
            HashSet<Integer> sigPeakCoveredSnp, nonsigPeakCoveredSnp;
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

                    sigPeakCoveredSnp = this.predSignificantPeakCoveredSNP.getOrDefault(chrNum, null);
                    if (sigPeakCoveredSnp != null && sigPeakCoveredSnp.contains(genomePosition)) {
                        this.significantAsePeakCoveredGeneId.add(geneId);
                    }
                    nonsigPeakCoveredSnp = this.predNonSignificantPeakCoveredSNP.getOrDefault(chrNum, null);
                    if (nonsigPeakCoveredSnp != null && nonsigPeakCoveredSnp.contains(genomePosition)) {
                        this.nonsignificantAsePeakCoveredGeneId.add(geneId);
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
     * 获取ASE显著的peak覆盖的SNP所在基因的ASE ratio
     * significantPeakCoveredGeneAse -> {gene1: majorAlleleRatio1, gene2: majorAlleleRatio2,...}
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
                    if (this.significantAsePeakCoveredGeneId.contains(geneID)) {
                        this.significantPeakCoveredGeneAse.put(geneName, majorAlleleRatio);
                    }
                    if (this.nonsignificantAsePeakCoveredGeneId.contains(geneID)) {
                        this.nonsignificantPeakCoveredGeneAse.put(geneName, majorAlleleRatio);
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

    /**
     * 得到模拟数据的ASE显著及非显著的m6A信号覆盖的基因是否为预设的ASE基因，并得到预测准确率
     * @return 预测准确率
     */
    public double getAseGenePredictResult() {
        this.parseSignificantAseTestResult();
        this.parsePeakCoveredSNPFile();
        this.getAseGeneId();
        this.getGeneAseRatio();
        double aseGeneAccuracy = 0;
        double geneNum = this.significantPeakCoveredGeneAse.size() + this.nonsignificantPeakCoveredGeneAse.size();
        for (Double ase: this.significantPeakCoveredGeneAse.values()) {
            if (ase > 0.59 && ase < 0.85)
                aseGeneAccuracy++;
        }
        for (Double ase: this.nonsignificantPeakCoveredGeneAse.values()) {
            if (ase == 0.5)
                aseGeneAccuracy ++;
        }
        aseGeneAccuracy = aseGeneAccuracy / geneNum;

        return aseGeneAccuracy;
    }


    /**
     * 获取ASE显著的m6A信号覆盖的基因上模拟peak的范围，即修饰位点前后 fragment length长度
     * simulatePeakIntervals -> {gene1: ([range1], [range2],...), gene2:([range1, range2,...)}
     */
    private void simulateM6aPeakInterval() {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.simulateM6aFile))));
            String line = "", geneId;
            String[] info;
            int position, lineNum = 1;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (lineNum == 1) {
                        lineNum++;
                        continue;
                    }
                    info = line.split("\t");
                    geneId = info[0];
                    if (!this.significantAsePeakCoveredGeneId.contains(geneId))
                        continue;
                    position = Integer.parseInt(info[2]);
                    HashSet<int[]> peakIntervals = this.simulatePeakIntervals.getOrDefault(geneId, new HashSet<>());
                    peakIntervals.add(new int[]{Math.max(0, position-this.fragmentLength), position+this.fragmentLength});
                    this.simulatePeakIntervals.put(geneId, peakIntervals);
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
     * 得到模拟数据中peak覆盖的SNP位点
     * simulatePeakCoveredSNP -> {chrNum: (site1, site2,...), ...}
     */
    private void realPeakCoveredSNPSite() {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.mutationFile))));
            String line = "", geneId, chrNum;
            String[] info;
            int lineNum = 1, position;
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
                    position = Integer.parseInt(info[2]);
                    HashSet<int[]> peakIntervals = this.simulatePeakIntervals.getOrDefault(geneId, null);
                    if (peakIntervals == null)
                        continue;
                    for (int[] range: peakIntervals) {
                        if (range[0] <= position && position <= range[1]) {
                            HashSet<Integer> peakCoveredSites = this.simulatePeakCoveredSNP.getOrDefault(chrNum, new HashSet<>());
                            peakCoveredSites.add(position);
                            this.simulatePeakCoveredSNP.put(chrNum, peakCoveredSites);
                        }
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
     * 预测结果文件中peak covered SNP中与模拟数据相同的部分，计算recall和accuracy
     * @return 准确率和召回率
     */
    public double[] getPeakCoveredSNPResult() {
        this.simulateM6aPeakInterval();
        this.realPeakCoveredSNPSite();
        int totalSimulateCoveredSNP = 0, totalPredCoveredSNP = 0;
        double truePositiveCoveredSNP = 0, peakCoverSNPRecall, peakCoverSNPAccuracy;
        for (String chrNum: this.simulatePeakCoveredSNP.keySet()) {
            HashSet<Integer> simulateCoveredSNP = this.simulatePeakCoveredSNP.get(chrNum);
            totalSimulateCoveredSNP += simulateCoveredSNP.size();
            HashSet<Integer> predCoveredSNP = this.predSignificantPeakCoveredSNP.getOrDefault(chrNum, null);
            if (predCoveredSNP == null)
                continue;
            totalPredCoveredSNP += predCoveredSNP.size();
            simulateCoveredSNP.retainAll(predCoveredSNP);
            truePositiveCoveredSNP += simulateCoveredSNP.size();
        }
        peakCoverSNPAccuracy = truePositiveCoveredSNP / totalPredCoveredSNP;
        peakCoverSNPRecall = truePositiveCoveredSNP / totalSimulateCoveredSNP;

        return new double[] {peakCoverSNPAccuracy, peakCoverSNPRecall};
    }


    /**
     * 获取peak calling得到的信息
     * peakCallingIntervals -> {gene1: ([range1], [range2],...), gene2:([range1, range2,...)}
     */
    public void parsePeakBedFile() {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.peakBedFile))));
            String line = "", geneId;
            String[] info, blockSizes, blockStarts;
            int peakStart, peakEnd;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    geneId = info[3];
                    peakStart = Integer.parseInt(info[1]);
                    peakEnd = Integer.parseInt(info[2]);
                    blockSizes = info[10].substring(0, info[10].length()-1).split(",");
                    blockStarts = info[11].split(",");

                    HashSet<int[]> peakIntervals = this.peakCallingIntervals.getOrDefault(geneId, new HashSet<>());
                    //                    for (int i=0; i<blockSizes.length; i++) {
                    //                        int startSite = peakStart+Integer.parseInt(blockStarts[i]);
                    //                        peakIntervals.add(new int[]{startSite, startSite+Integer.parseInt(blockSizes[i])});
                    //                    }
                    peakIntervals.add(new int[]{peakStart, peakEnd});
                    this.peakCallingIntervals.put(geneId, peakIntervals);
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
     * 获取显著性检验结果中每条染色体上ASM specific m6A signal(q<0.05)及非特异(q>=0.05)及非特异(q>)的peak范围
     * significantAsmChrPeaks -> {chrNum: (start1->end1, start2->end2,...)}
     */
    public void parseSignificantAsmTestResult() {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.significantAsmTestResult))));
            String line = "";
            String[] info;
            String chrNum, peakStart, peakEnd;
            double qVal;
            HashSet<String> peaks;
            while (line!=null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    qVal = Double.parseDouble(info[3]);
                    chrNum = info[0];
                    peakStart = info[1];
                    peakEnd = info[2];
                    if (qVal-0.05 > 0.000001) {
                        peaks = this.nonsignificantAsmChrPeaks.getOrDefault(chrNum, new HashSet<>());
                        peaks.add(String.join("->", new String[]{peakStart, peakEnd}));
                        this.nonsignificantAsmChrPeaks.put(chrNum, peaks);
                    } else {
                        peaks = this.significantAsmChrPeaks.getOrDefault(chrNum, new HashSet<>());
                        peaks.add(String.join("->", new String[]{peakStart, peakEnd}));
                        this.significantAsmChrPeaks.put(chrNum, peaks);
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
     * 将模拟的m6A位点定位到peak calling的peak中,并标记该peak中的m6A位点是否具有AS
     * {start->end -> (true, false,...), ...}
     */
    public void parseSimulateM6aFile() {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.simulateM6aFile))));
            String line = "", geneId, peak;
            String[] info;
            int lineNum = 1, position;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (lineNum == 1) {
                        lineNum++;
                        continue;
                    }
                    info = line.split("\t");
                    geneId = info[0];
                    position = Integer.parseInt(info[2]);

                    HashSet<int[]> genePeakRanges = this.peakCallingIntervals.getOrDefault(geneId, null);
                    if (genePeakRanges == null)
                        continue;
                    for (int[] range: genePeakRanges) {
                        if (range[0] <= position && position <= range[1]) {
                            peak = range[0] + "->" + range[1];
                            HashSet<Double> coveredM6a = this.callingPeakCoveredM6a.getOrDefault(peak, new HashSet<>());
                            coveredM6a.add(Double.parseDouble(info[3]));
                            this.callingPeakCoveredM6a.put(peak, coveredM6a);
                        }
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

    public double getAsmPeakPredResult() {
        this.parsePeakBedFile();
        this.parseSignificantAsmTestResult();
        this.parseSimulateM6aFile();
        int predSpecificPeak = 0, realSpecificPeak = 0;
        for (String chrNum: this.significantAsmChrPeaks.keySet()) {
            predSpecificPeak += this.significantAsmChrPeaks.get(chrNum).size();
        }

        for (String peak: this.callingPeakCoveredM6a.keySet()) {
            for (String chrNum: this.significantAsmChrPeaks.keySet()) {
                HashSet<String> significantPeak = this.significantAsmChrPeaks.get(chrNum);
                if (significantPeak.contains(peak)) {
                    HashSet<Double> majorAlleleSpecific = this.callingPeakCoveredM6a.get(peak);
                    for (Double ratio: majorAlleleSpecific) {
                        if (ratio != 0.5) {
                            realSpecificPeak++;
                            break;
                        }
                    }
                }
            }
        }
        return (double) (realSpecificPeak) / predSpecificPeak;
    }
}
