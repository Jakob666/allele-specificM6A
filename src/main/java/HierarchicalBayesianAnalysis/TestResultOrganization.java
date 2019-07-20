package HierarchicalBayesianAnalysis;

import java.io.*;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

/**
 * 对检验结果进行整理，依据ASE和ASM是否显著分为4类。
 */
public class TestResultOrganization {
    private String aseTestFile, asmTestFile, peakBedFile, inputVcfFile, ipVcfFile, asePeakCoveredSNPFile,
                   asmPeakCoveredSNPFile, outputFile;
    private HashMap<String, HashMap<String, Double>> aseM6aPeak, asmM6aPeak;
    public HashMap<String, String> peakCoveredGene = new HashMap<>();
    private HashMap<String, String> ipMajorAllele, inputMajorAllele;
    private HashMap<String, HashSet<String>> asePeakCoveredSNP, asmPeakCoveredSNP;
    public HashMap<String, String> significantPeak = new HashMap<>(), nonsignificantPeak = new HashMap<>();

    public TestResultOrganization(String aseTestFile, String asmTestFile, String peakBedFile, String inputVcfFile,
                                  String ipVcfFile, String asePeakCoveredSNPFile, String asmPeakCoveredSNPFile,
                                  String outputFile) {
        this.aseTestFile = aseTestFile;
        this.asmTestFile = asmTestFile;
        this.peakBedFile = peakBedFile;
        this.inputVcfFile = inputVcfFile;
        this.ipVcfFile = ipVcfFile;
        this.asePeakCoveredSNPFile = asePeakCoveredSNPFile;
        this.asmPeakCoveredSNPFile = asmPeakCoveredSNPFile;
        this.outputFile = outputFile;
    }

    public void organizePeaks() {
        this.aseM6aPeak = this.parseTestFile(this.aseTestFile);
        this.asmM6aPeak = this.parseTestFile(this.asmTestFile);
        this.parsePeakBedFile();
        this.ipMajorAllele = this.parseVcfFile(this.ipVcfFile);
        this.inputMajorAllele = this.parseVcfFile(this.inputVcfFile);
        this.asePeakCoveredSNP = this.parseCoveredSNPFile(this.asePeakCoveredSNPFile);
        this.asmPeakCoveredSNP = this.parseCoveredSNPFile(this.asmPeakCoveredSNPFile);
        this.classify();
        this.aseM6aPeak = this.parseTestFile(this.aseTestFile);
        this.asmM6aPeak = this.parseTestFile(this.asmTestFile);
        this.writeIn();
    }

    /**
     * 对Hierarchical model检验的输出文件进行解读
     * @param testFile 检验结果文件
     * @return 读取的结果 testResult -> {sig-> {chr:start:end-> qVal, ...}, unsig-> {chr:start:end-> qVal, ...}}
     */
    public HashMap<String, HashMap<String, Double>> parseTestFile(String testFile) {
        HashMap<String, HashMap<String, Double>> testResult = new HashMap<>();
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(testFile))));
            String line = "", chrNum, peakStart, peakEnd, label;
            String[] info;
            int lineNum = 1;
            double qVal;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (lineNum == 1) {
                        lineNum++;
                        continue;
                    }
                    info = line.split("\t");
                    chrNum = info[0];
                    peakStart = info[1];
                    peakEnd = info[2];
                    qVal = Double.parseDouble(info[3]);

                    label = String.join(":", new String[]{chrNum, peakStart, peakEnd});
                    if (qVal >= 0.05) {
                        HashMap<String, Double> unsignificantPeak = testResult.getOrDefault("unsig", new HashMap<>());
                        unsignificantPeak.put(label, qVal);
                        testResult.put("unsig", unsignificantPeak);
                    } else {
                        HashMap<String, Double> significantPeak = testResult.getOrDefault("sig", new HashMap<>());
                        significantPeak.put(label, qVal);
                        testResult.put("sig", significantPeak);
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

        return testResult;
    }

    /**
     * 将peak calling得到的bed文件中的peak与基因对应, {chr:peakStart:peakEnd -> geneId, chr:peakStart:peakEnd -> geneId,...}
     */
    public void parsePeakBedFile() {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.peakBedFile))));
            String line = "", chrNum, peakStart, peakEnd, geneId, label;
            String[] info;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    chrNum = info[0];
                    peakStart = info[1];
                    peakEnd = info[2];
                    geneId = info[3];
                    label = String.join(":", new String[]{chrNum, peakStart, peakEnd});
                    this.peakCoveredGene.put(label, geneId);
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
     * 从IP和INPUT样本SNP calling的VCF文件中得到各位点major allele的信息
     * @param vcfFile VCF文件
     * @return VCF文件中major allele的信息 {chr:position -> ref, chr:position -> alt, ...}
     */
    public HashMap<String, String> parseVcfFile(String vcfFile) {
        HashMap<String, String> majorAlleleRecord = new HashMap<>();
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(vcfFile))));
            String line = "", chrNum, position, dp4 = "", maj, label;
            String[] info;
            boolean refIsMajor;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    chrNum = info[0];
                    position = info[1];
                    label = chrNum + ":" + position;
                    for (String s: info[7].split(";")) {
                        if (s.startsWith("DP4"))
                            dp4 = s;
                    }
                    refIsMajor = this.parseDp4(dp4);
                    if (refIsMajor)
                        maj = "ref";
                    else
                        maj = "alt";
                    majorAlleleRecord.put(label, maj);
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

        return majorAlleleRecord;
    }

    /**
     * 从VCF文件的DP4信息推断major allele是 reference gene还是 alternative gene
     * @param dp4String VCF文件DP4信息
     * @return 如果reference是major allele则返回true
     */
    private boolean parseDp4(String dp4String) {
        String[] readsCount = dp4String.split("=")[1].split(",");
        int refCount = Integer.parseInt(readsCount[0]) + Integer.parseInt(readsCount[1]);
        int altCount = Integer.parseInt(readsCount[2]) + Integer.parseInt(readsCount[3]);

        return refCount >= altCount;
    }

    /**
     * 通过peak covered记录文件得到每个peak下覆盖的SNP位点
     * @param coveredSNPFile cover snp文件
     * @return 每个peak下覆盖位点的记录 {chr:peakStart:peakEnd -> (chr:position,...), ...}
     */
    public HashMap<String, HashSet<String>> parseCoveredSNPFile(String coveredSNPFile) {
        HashMap<String, HashSet<String>> coveredSNPRecord = new HashMap<>();
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(coveredSNPFile))));
            String line = "", chrNum, position, peakStart, peakEnd, label;
            String[] info;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    chrNum = info[0];
                    position = info[2];
                    peakStart = info[3];
                    peakEnd = info[4];

                    label = String.join(":", new String[]{chrNum, peakStart, peakEnd});
                    HashSet<String> snpSites = coveredSNPRecord.getOrDefault(label, new HashSet<>());
                    snpSites.add(chrNum+":"+position);
                    coveredSNPRecord.put(label, snpSites);
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

        return coveredSNPRecord;
    }

    /**
     * 依据ASE和ASM的检验结果将m6A peak分为四类，均显著、均不显著和其一显著
     */
    private void classify() {
        // ASE和ASM均显著的peak
        Set<String> aseSignificantPeaks = this.aseM6aPeak.get("sig").keySet();
        Set<String> asmSignificantPeaks = this.asmM6aPeak.get("sig").keySet();
        aseSignificantPeaks.retainAll(asmSignificantPeaks);
        for (String peak: aseSignificantPeaks) {
            HashSet<String> aseSNP = this.asePeakCoveredSNP.get(peak);
            HashSet<String> asmSNP = this.asmPeakCoveredSNP.get(peak);
            if (!this.aseEqualAsm(aseSNP, asmSNP))
                this.significantPeak.put(peak, this.peakCoveredGene.get(peak));
            else
                this.nonsignificantPeak.put(peak, this.peakCoveredGene.get(peak));
        }
        // ASE和ASM均不显著的peak
        Set<String> aseNonsignificantPeaks = this.aseM6aPeak.get("unsig").keySet();
        Set<String> asmNonsignificantPeaks = this.asmM6aPeak.get("unsig").keySet();
        aseNonsignificantPeaks.retainAll(asmNonsignificantPeaks);
        for (String peak: aseNonsignificantPeaks) {
            this.nonsignificantPeak.put(peak, this.peakCoveredGene.get(peak));
        }
        // ASE显著而ASM不显著
        aseSignificantPeaks = this.aseM6aPeak.get("sig").keySet();
        asmNonsignificantPeaks = this.asmM6aPeak.get("unsig").keySet();
        aseSignificantPeaks.retainAll(asmNonsignificantPeaks);
        for (String peak: aseSignificantPeaks) {
            this.nonsignificantPeak.put(peak, this.peakCoveredGene.get(peak));
        }
        // ASE不显著而ASM显著
        aseNonsignificantPeaks = this.aseM6aPeak.get("unsig").keySet();
        asmSignificantPeaks = this.asmM6aPeak.get("sig").keySet();
        aseNonsignificantPeaks.retainAll(asmSignificantPeaks);
        for (String peak: aseNonsignificantPeaks) {
            this.significantPeak.put(peak, this.peakCoveredGene.get(peak));
        }
    }

    /**
     * 判断某个peak覆盖范围内ASE和ASM的major allele是否相同
     * @param aseSNP peak覆盖的ASE位点
     * @param asmSNP peak覆盖的ASM位点
     * @return 如果major allele相同则返回true
     */
    private boolean aseEqualAsm(HashSet<String> aseSNP, HashSet<String> asmSNP) {
        String aseStrand, asmStrand;
        int ref = 0, alt = 0;
        for (String snpSite: aseSNP) {
            if (this.inputMajorAllele.get(snpSite).equals("ref"))
                ref++;
            else
                alt++;
        }
        if (ref > alt)
            aseStrand = "ref";
        else
            aseStrand = "alt";
        ref = 0;
        alt = 0;
        for (String snpSite: asmSNP) {
            if (this.ipMajorAllele.get(snpSite).equals("ref"))
                ref++;
            else
                alt++;
        }
        if (ref > alt)
            asmStrand = "ref";
        else
            asmStrand = "alt";

        return aseStrand.equals(asmStrand);
    }

    /**
     * 结果写入文件
     */
    private void writeIn() {
        BufferedWriter bfw = null;
        try {
            bfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(this.outputFile))));
            bfw.write("#chr\tpeakStart\tpeakEnd\tcoveredGene\tASE q-value\tASM q-value\tsignificant\n");
            String chrNum, peakStart, peakEnd, geneName, line;
            Double aseQVal, asmQVal;
            String[] info;
            for (String peak: this.significantPeak.keySet()) {
                info = peak.split(":");
                chrNum = info[0];
                peakStart = info[1];
                peakEnd = info[2];

                geneName = this.peakCoveredGene.get(peak);
                aseQVal = this.aseM6aPeak.get("sig").getOrDefault(peak, null);
                if (aseQVal ==null)
                    aseQVal = this.aseM6aPeak.get("unsig").get(peak);
                asmQVal = this.asmM6aPeak.get("sig").get(peak);
                line = String.join("\t", new String[]{chrNum, peakStart, peakEnd, geneName,
                        Double.toString(aseQVal), Double.toString(asmQVal), Boolean.toString(true)});
                bfw.write(line);
                bfw.newLine();
            }
            for (String peak: this.nonsignificantPeak.keySet()) {
                info = peak.split(":");
                chrNum = info[0];
                peakStart = info[1];
                peakEnd = info[2];
                geneName = this.peakCoveredGene.get(peak);
                aseQVal = this.aseM6aPeak.get("sig").getOrDefault(peak, null);
                if (aseQVal == null)
                    aseQVal = this.aseM6aPeak.get("unsig").get(peak);
                asmQVal = this.asmM6aPeak.get("sig").getOrDefault(peak, null);
                if (asmQVal == null)
                    asmQVal = this.asmM6aPeak.get("unsig").get(peak);
                bfw.write(String.join("\t", new String[]{chrNum, peakStart, peakEnd, geneName,
                        Double.toString(aseQVal), Double.toString(asmQVal), Boolean.toString(false)}));
                bfw.newLine();
            }
            bfw.close();
        } catch (IOException ie) {
            ie.printStackTrace();
        } finally {
            if (bfw != null) {
                try {
                    bfw.close();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
    }
}
