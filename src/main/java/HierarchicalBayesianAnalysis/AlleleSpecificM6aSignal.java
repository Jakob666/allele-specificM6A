package HierarchicalBayesianAnalysis;

import org.apache.log4j.Logger;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

public class AseSpecificM6aSignal {
    private File aseTestFile, asmTestFile, outputFile;
    private HashMap<String, String[]> aseSignificantGene = new HashMap<>(), aseNonsignificantGene = new HashMap<>(),

                                      significantResult = new HashMap<>(), nonSignificantResult = new HashMap<>();
    private HashMap<String, ArrayList<String[]>> asmSignificantPeak = new HashMap<>(), asmNonsignificantPeak = new HashMap<>();
    private Logger log;
    private DecimalFormat df = new DecimalFormat("0.0000");

    /**
     * Constructor
     * @param aseTestFile ASE检验结果文件
     * @param asmTestFile ASM检验结果文件
     * @param outputFile 输出文件
     * @param log Log4j Logger对象
     */
    public AseSpecificM6aSignal(String aseTestFile, String asmTestFile, String outputFile, Logger log) {
        this.aseTestFile = new File(aseTestFile);
        this.asmTestFile = new File(asmTestFile);
        this.outputFile = new File(outputFile);
        this.log = log;
    }

    /**
     * 检验ASE特异的m6A信号并将其写入文件
     */
    public void detect() {
        this.parseAseTestFile();
        this.parseAsmTestFile();
        this.mergeAseAsmTestResult();
        this.outputResult();
    }

    /**
     * 获取ASE Gene检验结果
     */
    private void parseAseTestFile() {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(this.aseTestFile)));
            String line = "", geneId, geneName, majorAlleleStrand;
            String[] info, record;
            double qVal;
            while (line != null) {
                line = bfr.readLine();
                if (line != null && !line.startsWith("#")) {
                    info = line.split("\t");
                    geneId = info[0];
                    geneName = info[1];
                    qVal = Double.parseDouble(info[2]);
                    majorAlleleStrand = info[5];
                    record = new String[] {geneName, info[2], majorAlleleStrand};
                    if (qVal - 0.05 < 0.000001)
                        this.aseSignificantGene.put(geneId, record);
                    else
                        this.aseNonsignificantGene.put(geneId, record);
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
     * 获取ASM Peak检验结果
     */
    private void parseAsmTestFile() {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(this.asmTestFile)));
            String line = "", chrNum, peakStart, peakEnd, geneId, majorAlleleStrand;
            String[] info, record;
            ArrayList<String[]> recordList;
            double qVal;
            while (line != null) {
                line = bfr.readLine();
                if (line != null && !line.startsWith("#")) {
                    info = line.split("\t");
                    chrNum = info[0];
                    peakStart = info[1];
                    peakEnd = info[2];
                    geneId = info[3];
                    qVal = Double.parseDouble(info[4]);
                    majorAlleleStrand = info[7];

                    record = new String[] {chrNum, peakStart, peakEnd, info[4], majorAlleleStrand};

                    if (qVal - 0.05 < 0.00001) {
                        recordList = this.asmSignificantPeak.getOrDefault(geneId, new ArrayList<>());
                        recordList.add(record);
                        this.asmSignificantPeak.put(geneId, recordList);
                    } else {
                        recordList = this.asmNonsignificantPeak.getOrDefault(geneId, new ArrayList<>());
                        recordList.add(record);
                        this.asmNonsignificantPeak.put(geneId, recordList);
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
     * 合并ASE和ASM检验的结果得到 ASE specific m6A peak
     */
    private void mergeAseAsmTestResult() {
        ArrayList<String[]> asmRecordList;
        String[] aseRecord, combineRecord;
        String aseMajorAllele, asmMajorAllele;
        Set<String> aseGeneId = this.aseSignificantGene.keySet();
        Set<String> asmGeneId = this.asmSignificantPeak.keySet();
        // ASE和ASM均显著的结果
        aseGeneId.retainAll(asmGeneId);
        for (String geneId: aseGeneId) {
            aseRecord = this.aseSignificantGene.get(geneId);
            asmRecordList = this.asmSignificantPeak.get(geneId);
            aseMajorAllele = aseRecord[2];
            for (String[] asmRecord: asmRecordList) {
                asmMajorAllele = asmRecord[4];
                // chrNum, geneName, peakStart, peakEnd, ASE q-value, ASM q-value
                combineRecord = new String[] {asmRecord[0], aseRecord[0], asmRecord[1], asmRecord[2], aseRecord[1], asmRecord[3]};
                if (aseMajorAllele.equals(asmMajorAllele))
                    this.nonSignificantResult.put(geneId, combineRecord);
                else
                    this.significantResult.put(geneId, combineRecord);
            }
        }

        // ASE显著而ASM不显著的结果
        aseGeneId = this.aseSignificantGene.keySet();
        Set<String> nonAsmGeneId = this.asmNonsignificantPeak.keySet();
        nonAsmGeneId.retainAll(aseGeneId);
        for (String geneId: nonAsmGeneId) {
            aseRecord = this.aseSignificantGene.get(geneId);
            asmRecordList = this.asmNonsignificantPeak.get(geneId);
            for (String[] asmRecord: asmRecordList) {
                combineRecord = new String[] {asmRecord[0], aseRecord[0], asmRecord[1], asmRecord[2], aseRecord[1], asmRecord[3]};
                this.nonSignificantResult.put(geneId, combineRecord);
            }
        }

        // ASE不显著而ASM显著的结果
        Set<String> nonAseGeneId = this.aseNonsignificantGene.keySet();
        asmGeneId.retainAll(nonAseGeneId);
        for (String geneId: asmGeneId) {
            aseRecord = this.aseNonsignificantGene.get(geneId);
            asmRecordList = this.asmSignificantPeak.get(geneId);
            for (String[] asmRecord: asmRecordList) {
                combineRecord = new String[] {asmRecord[0], aseRecord[0], asmRecord[1], asmRecord[2], aseRecord[1], asmRecord[3]};
                this.significantResult.put(geneId, combineRecord);
            }
        }

        // ASE和ASM均不显著的结果
        nonAseGeneId = this.aseNonsignificantGene.keySet();
        nonAsmGeneId = this.asmNonsignificantPeak.keySet();
        nonAseGeneId.retainAll(nonAsmGeneId);
        for (String geneId: nonAseGeneId) {
            aseRecord = this.aseNonsignificantGene.get(geneId);
            asmRecordList = this.asmNonsignificantPeak.get(geneId);
            for (String[] asmRecord: asmRecordList) {
                combineRecord = new String[] {asmRecord[0], aseRecord[0], asmRecord[1], asmRecord[2], aseRecord[1], asmRecord[3]};
                this.nonSignificantResult.put(geneId, combineRecord);
            }
        }
    }

    /**
     * 将合并后的结果写入到文件
     */
    private void outputResult() {
        BufferedWriter bfw = null;
        try {
            bfw = new BufferedWriter(
                    new OutputStreamWriter(new FileOutputStream(this.outputFile))
            );

            // 依据ASE的p值从小到大排序
            ArrayList<Map.Entry<String, String[]>> sigResult = new ArrayList<>(this.significantResult.entrySet());
            Collections.sort(sigResult, new Comparator<Map.Entry<String, String[]>>() {
                public int compare(Map.Entry<String, String[]> o1,
                                   Map.Entry<String, String[]> o2) {
                    Double val1 = Double.parseDouble(o1.getValue()[4]);
                    Double val2 = Double.parseDouble(o2.getValue()[4]);
                    return val1.compareTo(val2);
                }
            });

            ArrayList<Map.Entry<String, String[]>> nonsigResult = new ArrayList<>(this.nonSignificantResult.entrySet());
            Collections.sort(nonsigResult, new Comparator<Map.Entry<String, String[]>>() {
                public int compare(Map.Entry<String, String[]> o1,
                                   Map.Entry<String, String[]> o2) {
                    Double val1 = Double.parseDouble(o1.getValue()[4]);
                    Double val2 = Double.parseDouble(o2.getValue()[4]);
                    return val1.compareTo(val2);
                }
            });

            String geneId, chrNum, geneName, peakStart, peakEnd, aseQVal, asmQVal, line;
            String[] record;
            bfw.write("#chr\tgeneName\tgeneId\tpeakStart\tpeakEnd\tASE q-value\tASM q-value\tsignificant\n");
            for (Map.Entry<String, String[]> result: sigResult) {
                geneId = result.getKey();
                // chrNum, geneName, peakStart, peakEnd, ASE q-value, ASM q-value
                record = result.getValue();
                chrNum = record[0];
                geneName = record[1];
                peakStart = record[2];
                peakEnd = record[3];
                aseQVal = record[4];
                asmQVal = record[5];
                line = String.join("\t", new String[] {chrNum, geneName, geneId, peakStart, peakEnd, aseQVal,
                                                                  asmQVal, Boolean.toString(true)});
                bfw.write(line);
                bfw.newLine();
            }
            for (Map.Entry<String, String[]> result: nonsigResult) {
                geneId = result.getKey();
                // chrNum, geneName, peakStart, peakEnd, ASE q-value, ASM q-value
                record = result.getValue();
                chrNum = record[0];
                geneName = record[1];
                peakStart = record[2];
                peakEnd = record[3];
                aseQVal = record[4];
                asmQVal = record[5];
                line = String.join("\t", new String[] {chrNum, geneName, geneId, peakStart, peakEnd, aseQVal,
                        asmQVal, Boolean.toString(false)});
                bfw.write(line);
                bfw.newLine();
            }
            bfw.close();
        } catch (IOException ie) {
            this.log.error("can not write ASE m6a peak result file");
            this.log.error(ie.getMessage());
        } finally {
            if (bfw != null) {
                try {
                    bfw.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }
}
