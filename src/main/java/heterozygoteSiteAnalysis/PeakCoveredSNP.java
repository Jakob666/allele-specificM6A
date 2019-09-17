package heterozygoteSiteAnalysis;

import org.apache.log4j.Logger;

import java.io.*;
import java.util.HashMap;

public class PeakCoveredSNP {
    private File vcfFile, peakCallingRes, outputFile;
    int readsInfimum;
    private Logger logger;

    /**
     * Constructor
     * @param vcfRecordFile vcf record file
     * @param peakCallingRes m6a peak calling result bed file
     * @param outputFile 输出文件路径
     * @param readsInfimum SNP reads 对peak覆盖的SNP位点进行记录时reads数目下确界
     * @param logger Logger instance
     */
    public PeakCoveredSNP(String vcfRecordFile, String peakCallingRes, String outputFile, int readsInfimum, Logger logger) {
        this.logger = logger;
        this.vcfFile = new File(vcfRecordFile);
        if (!vcfFile.exists() || !this.vcfFile.isFile()) {
            this.logger.error("vcf file not exists");
            System.exit(2);
        }
        this.peakCallingRes = new File(peakCallingRes);
        if (!this.peakCallingRes.exists() || !this.peakCallingRes.isFile()) {
            this.logger.error("peak calling result bed file not exists");
            System.exit(2);
        }
        this.readsInfimum = readsInfimum;
        this.outputFile = new File(outputFile);
    }

    /**
     * 筛选覆盖SNP的peak并记录SNP信息，输出到文件
     */
    public void filterSNPAndPeak() {
        HashMap<String, HashMap<String, IntervalTree>> m6aTreeMap = this.getM6aPeakTree();
        BufferedReader bfr = null;
        BufferedWriter bfw = null;
        try {
            bfr = new BufferedReader(
                    new InputStreamReader(new FileInputStream(this.vcfFile))
            );
            bfw = new BufferedWriter(
                    new OutputStreamWriter(new FileOutputStream(this.outputFile))
            );

            String line = "", writeOut;
            String[] info;
            String chrNum, refNc, altNc, majorCount, minorCount, majorAlleleStrand, majorNc, minorNc;
            int[] refAndAltCount;
            int position;
            bfw.write("#chr\tSNP strand\tposition\tpeakStart\tpeakEnd\tmajorAlleleStrand\tmajorNc\tminorNc\tmajorCount\tminorCount\n");
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    // 获取一条记录中的SNP的相关信息
                    info = line.split("\t");
                    chrNum = info[0];
                    refNc = info[3];
                    altNc = info[4];
                    HashMap<String, IntervalTree> chrTree = m6aTreeMap.getOrDefault(chrNum, null);
                    if (chrTree == null)
                        continue;
                    // 只保留单核苷酸突变
                    if (refNc.length() > 1 | altNc.length() > 1)
                        continue;
                    position = Integer.parseInt(info[1]);
                    refAndAltCount = this.getReadsCountViaDp4(info[7]);
                    // 设置阈值防止出现假阳性位点
                    if (refAndAltCount[1] == 0) // refAndAltCount[0] == 0 |
                        continue;
                    // 通过设定的阈值对SNP进行筛选，减少误差
                    if (Math.max(refAndAltCount[0], refAndAltCount[1]) < this.readsInfimum)
                        continue;
                    if (refAndAltCount[0] >= refAndAltCount[1]) {
                        majorCount = Integer.toString(refAndAltCount[0]);
                        minorCount = Integer.toString(refAndAltCount[1]);
                        majorAlleleStrand = "ref";
                        majorNc = refNc;
                        minorNc = altNc;
                    } else {
                        majorCount = Integer.toString(refAndAltCount[1]);
                        minorCount = Integer.toString(refAndAltCount[0]);
                        majorAlleleStrand = "alt";
                        majorNc = altNc;
                        minorNc = refNc;
                    }

                    // 获取对应的染色体的区间树并在正负链的peak上搜寻是否被覆盖。如果被覆盖则写入文件
                    IntervalTree posStrandTree = chrTree.get("+");
                    IntervalTree negStrandTree = chrTree.get("-");
                    IntervalTreeNode posStrandSearchResult = posStrandTree.search(posStrandTree.root, position);
                    IntervalTreeNode negStrandSearchResult = negStrandTree.search(negStrandTree.root, position);
                    if (posStrandSearchResult != null) {
                        String[] newLine = new String[]{chrNum, "+", info[1], Integer.toString(posStrandSearchResult.peakStart),
                                                        Integer.toString(posStrandSearchResult.peakEnd), majorAlleleStrand, majorNc, minorNc, majorCount, minorCount};
                        writeOut = String.join("\t", newLine);
                        bfw.write(writeOut);
                        bfw.newLine();
                    }else if (negStrandSearchResult != null) {
                        String[] newLine = new String[]{chrNum, "-", info[1], Integer.toString(negStrandSearchResult.peakStart),
                                                        Integer.toString(negStrandSearchResult.peakEnd), majorAlleleStrand, majorNc, minorNc, majorCount, minorCount};
                        writeOut = String.join("\t", newLine);
                        bfw.write(writeOut);
                        bfw.newLine();
                    }
                }
            }
            bfw.flush();
        } catch (IOException ie) {
            this.logger.error("load peak calling result information failed.");
            this.logger.error(ie.getMessage());
        } finally {
            if (bfr != null) {
                try {
                    bfr.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
            if (bfw != null) {
                try {
                    bfw.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    /**
     * 从peak calling的结果得到m6A信号的区间树
     * @return 每条染色体正负链上的peak分别建立区间树构成的哈希表 [chrNum: [Strand: interval tree]]
     */
    private HashMap<String, HashMap<String, IntervalTree>> getM6aPeakTree() {
        PeakIntervalTree pit = new PeakIntervalTree(this.peakCallingRes.getAbsolutePath(), this.logger);
        HashMap<String, HashMap<String, IntervalTree>> treeMap = pit.getPeakTrees();
        pit = null;

        return treeMap;
    }

    /**
     * 从VCF文件的DP4内容中提取出ref和alt的reads count
     * @param record VCF内容列
     * @return [refReadsCount, altReadsCount]
     */
    private int[] getReadsCountViaDp4(String record) {
        // record的形式大致为 DP=5;SGB=-0.379885;RPB=1;MQB=1;MQSB=1;BQB=1;MQ0F=0;ICB=1;HOB=0.5;AC=0;AN=0;DP4=1,3,0,1;MQ=60
        String[] data = record.split(";");
        String key = null, value = null;
        String[] kAndv;
        for (String kv: data) {
            kAndv = kv.split("=");
            key = kAndv[0];
            if (key.equals("DP4")) {
                value = kAndv[1];
                break;
            }
        }
        if (value == null) {
            this.logger.error("can not find DP4 field in VCF file");
            System.exit(2);
        }
        String[] reads = value.split(",");
        int refReadsCount = Integer.parseInt(reads[0]) + Integer.parseInt(reads[1]);
        int altReadsCount = Integer.parseInt(reads[2]) + Integer.parseInt(reads[3]);

        return new int[]{refReadsCount, altReadsCount};
    }
}
