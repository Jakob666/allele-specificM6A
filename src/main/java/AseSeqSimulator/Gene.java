package AseSeqSimulator;

import GTFComponent.ElementRecord;
import GTFComponent.TranscriptRecord;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.UniformIntegerDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;

/**
 * random fragment a gene transcript sequence, the fragment length obeys normal distribution
 */
public class Gene {
    private String geneId, strand, geneName;
    private TranscriptRecord longestTranscriptRecord;
    private double rnaRPKM, dnaRPKM;
    private String chr, exonSeq;
    private int rnaReadsCount, dnaReadsCount;
    private double[] pmRange = new double[]{0.8, 0.9};
    private ElementRecord exonList = null;
    private HashMap<Integer, ArrayList<int[]>> m6aSiteFragments = new HashMap<>(), m6aSiteMutateFragments = new HashMap<>();
    // 对应的read覆盖SNP位点的fragment
    private HashMap<Integer, ArrayList<int[]>> inputMutateFragments = new HashMap<>(), ipMutateFragments = new HashMap<>();

    // 位置相同的reads的fragment在exon上面终止位点
    private ArrayList<int[]> inputFragment = new ArrayList<>(), ipFragment = new ArrayList<>(), m6aSiteNormalFragments = new ArrayList<>();

    public Gene(String geneId, String geneName, String strand, String chr) {
        this.geneId = geneId;
        this.strand = strand;
        this.geneName = geneName;
        this.chr = chr;
    }

    /**
     * get the longest transcript of a certain gene, and set it a property for a Gene instance
     * @param longestTranscriptRecord TranscriptRecord instance with the longest sequence length
     */
    public void setLongestTranscriptRecord(TranscriptRecord longestTranscriptRecord) {
        this.longestTranscriptRecord = longestTranscriptRecord;
    }

    public void setPmRange(double lower, double upper) {
        this.pmRange = new double[]{lower, upper};
    }

    public void setRnaRPKM(double RPKM) {
        this.rnaRPKM = RPKM;
    }

    public void setDnaRPKM(double RPKM) {
        this.dnaRPKM = RPKM;
    }

    public String getGeneId() {
        return this.geneId;
    }

    public String getGeneName() {
        return this.geneName;
    }

    public String getStrand() {
        return this.strand;
    }

    public String getChr() {
        return this.chr;
    }

    public int getRnaReadsCount() {
        return this.rnaReadsCount;
    }

    public int getDnaReadsCount() {
        return this.dnaReadsCount;
    }

    public double getDnaRPKM() {
        return this.dnaRPKM;
    }

    public ElementRecord getExonList() {
        return this.exonList;
    }

    public String getExonSeq() {
        return this.exonSeq;
    }

    public TranscriptRecord getLongestTranscriptRecord() {
        return this.longestTranscriptRecord;
    }

    /**
     * 将基因的外显子区域拼合得到完整的外显子序列
     * @param twoBit TwoBitParser instance
     */
    public void splicing(TwoBitParser twoBit) {
        HashMap<String, ElementRecord> elementsOnTranscript = longestTranscriptRecord.getElementList();
        ElementRecord transcriptExon = elementsOnTranscript.getOrDefault("exon", null);
        StringBuilder sb = new StringBuilder();
        if (transcriptExon != null)
            this.exonList = transcriptExon;
        try {
            while (transcriptExon != null) {
                String strand = transcriptExon.getStrand();
                int exonStart = transcriptExon.getElementStart();
                int exonEnd = transcriptExon.getElementEnd();
                twoBit.reset();
                String exonSeq = twoBit.loadFragment(exonStart-1, (exonEnd - exonStart + 1));
                if (strand.equals("-")) {
                    sb.append(CommonMethod.AntiChain(exonSeq));
                } else {
                    sb.append(exonSeq);
                }
                transcriptExon = transcriptExon.getNextElement();
            }
        } catch (IOException ie) {
            ie.printStackTrace();
        }
        this.exonSeq = sb.toString();
        sb = null;
    }

    /**
     * 依据RPKM值和文库大小计算该基因的reads count
     * @param librarySize library size
     */
    public void calculateReadsCountViaLibrarySize(long librarySize) {
        this.rnaReadsCount = (int) ((this.exonSeq.length() * librarySize * this.rnaRPKM) / Math.pow(10, 9));
    }

    /**
     * 通过设定的测序深度和文库大小计算该基因的reads count
     * @param depth sequencing depth
     * @param readLength read length
     */
    public void calculateReadsCountViaSequencingDepth(int depth, int readLength) {
        this.dnaReadsCount = this.exonSeq.length() * depth / readLength;
    }

    /**
     * 依据m6A的富集率计算IP样本的reads count
     * @param pm m6A富集率
     * @param backgroundSize 背景大小
     */
    private int m6aReadsCount(double pm, int backgroundSize) {
        return  (int) (backgroundSize * (pm / (1 - pm)));
    }

    /**
     * 随机抽取reads count条fragment，将其富集到INPUT组
     * @param fragmentMean fragment平均长度
     * @param fragmentTheta fragment长度标准差
     */
    public void enrichInputFragment(int fragmentMean, int fragmentTheta, int readLength, Set<Integer> mutatePositions,
                                    String type) {
        int readsCount = (type.equals("rna"))? this.rnaReadsCount: this.dnaReadsCount;
        List<Integer> mutations = null;
        if (mutatePositions != null) {
            mutations = new ArrayList<>(mutatePositions);
            Collections.sort(mutations);
        }
        // 用于随机生成fragment的长度
        NormalDistribution nordi = new NormalDistribution(fragmentMean, fragmentTheta);
        int curReadsCount = 0, fragmentLength, breakPoint, endPoint;
        boolean cover;
        while (curReadsCount < readsCount) {
            cover = false;
            // 随机生成fragment的长度并确定fragment在外显子序列上的起始终止位点
            fragmentLength = Math.abs((int) nordi.sample());
            int[] breakAndEndPoint = this.getBreakEndPoint(fragmentLength);
            breakPoint = breakAndEndPoint[0];
            endPoint = breakAndEndPoint[1];

            // 检验fragment对应的read是否覆盖ASE位点，如果覆盖则需对其进行记录
            if (mutations != null)
                cover = this.ifCoverMutateStie(breakPoint, endPoint, readLength, mutations, this.inputMutateFragments);
            if (!cover)
                this.inputFragment.add(breakAndEndPoint);
            curReadsCount++;
        }
        if (mutations != null)
            mutations.clear();
        nordi = null;
    }

    /**
     * 随机抽取reads count条fragment，将其富集到IP组
     * @param fragmentMean fragment平均长度
     * @param fragmentTheta fragment长度标准差
     * @param geneM6aSites 该基因上m6A修饰位点
     */
    public void enrichIpFragment(int fragmentMean, int fragmentTheta, int readLength, Set<Integer> geneM6aSites,
                                 Set<Integer> mutatePositions, String type) {
        int readsCount = (type.equals("rna"))? rnaReadsCount: dnaReadsCount;
        List<Integer> mutations = null;
        if (mutatePositions != null) {
            mutations = new ArrayList<>(mutatePositions);
            Collections.sort(mutations);
        }
        List<Integer> m6aSites = new ArrayList<>(geneM6aSites);
        Collections.sort(m6aSites);
        NormalDistribution nordi = new NormalDistribution(fragmentMean, fragmentTheta);
        int curReadsCount = 0, fragmentLength, break_point, end_point;
        boolean cover;
        while (curReadsCount < readsCount) {
            cover = false;
            // 随机生成fragment的长度并确定fragment在外显子序列上的起始终止位点
            fragmentLength = Math.abs((int) nordi.sample());
            int[] breakAndEndPoint = this.getBreakEndPoint(fragmentLength);
            break_point = breakAndEndPoint[0];
            end_point = breakAndEndPoint[1];

            // 检查fragment对应的read是否覆盖ASE位点，如果覆盖则需对其进行记录
            if (mutations != null)
                cover = this.ifCoverMutateStie(break_point, end_point, readLength, mutations, this.ipMutateFragments);
            if (!cover)
                this.ipFragment.add(breakAndEndPoint);

            // 检查Fragment是否覆盖了甲基化位点，如果覆盖则富集到IP样本中；反之，则抛弃
            if (geneM6aSites.size() != 0)
                this.ifCoverM6aSite(break_point, end_point, m6aSites);
            curReadsCount++;
        }
        m6aSites.clear();
        if (mutations != null)
            mutations.clear();
        nordi = null;
    }

    /**
     * 生成测序样本
     * @param mateFile1 测序文件
     * @param mateFile2 测序文件(只有在pair-end时候指定，否则为null)
     * @param readLength 测序read length
     * @param fragmentMean fragment长度均值
     * @param fragmentTheta fragment长度标准差
     * @param direct 测序方式 single-end 或 pair-end
     * @param mutExonSeq 当基因具有ASE位点时传入突变的序列，不存在时为null
     * @param refProp major allele的频率
     * @param seqErrorModel 测序误差模型
     * @param sample ip或input
     * @param geneM6aAsm 如果基因上存在ASE位点，则传入各位点的ASM ratio
     */
    public void generateReads(BufferedWriter mateFile1, BufferedWriter mateFile2, int readLength, int fragmentMean,
                              int fragmentTheta, double refProp, String mutExonSeq, String direct, SequencingError seqErrorModel,
                              String sample, HashMap<Integer, Double> geneM6aAsm, HashMap<Integer, Boolean> minorBias) {
        ArrayList<int[]> fragmentRanges;
        ArrayList<Integer> m6aSites = null;
        if (geneM6aAsm != null) {
            m6aSites = new ArrayList<>(geneM6aAsm.keySet());
            Collections.sort(m6aSites);
        }
        HashMap<Integer, ArrayList<int[]>> mutateFragmentRanges;
        fragmentRanges = (sample.equals("ip"))? this.ipFragment : this.inputFragment;
        mutateFragmentRanges = (sample.equals("ip"))? this.ipMutateFragments : this.inputMutateFragments;

        try {
            // 生成一般fragment的reads(reads不覆盖SNP位点)
            this.writeIn(fragmentRanges, direct, readLength, seqErrorModel, mateFile1, mateFile2, this.exonSeq, "ref");

            // 生成覆盖突变位点的fragment的reads
            ArrayList<int[]> mutateFragments;
            int count, majorAlleleCount;
            List<int[]> majorAlleleFragmentRanges, minorAlleleFragmentRanges;
            for (Integer mutateSite: mutateFragmentRanges.keySet()) {
                // 获取每个ASE位点覆盖的fragment的范围
                mutateFragments = mutateFragmentRanges.get(mutateSite);
                if (sample.equals("ip") && m6aSites != null) {
                     int site = this.ifInM6aPeakRange(mutateSite, m6aSites, fragmentTheta+2*fragmentTheta);
                     if (site != -1)
                         refProp = geneM6aAsm.get(site);
                }
                count = mutateFragments.size();
                if (Math.abs(refProp - 0.5) < 0.00001) {
                    if (count % 2 == 0)
                        majorAlleleCount = count / 2;
                    else
                        majorAlleleCount = (count - 1) / 2;
                } else {
                    majorAlleleCount = (int) (count * refProp);
                }
                boolean bias = false;
                if (sample.equals("ip") && m6aSites != null) {
                    int site = this.ifInM6aPeakRange(mutateSite, m6aSites, fragmentMean + 2 * fragmentTheta);
                    if (site != -1) {
                        refProp = geneM6aAsm.get(site);
                        bias = minorBias.get(site);
                        if (!bias)
                            majorAlleleCount = count - majorAlleleCount;
                    }
                }
                Collections.shuffle(mutateFragments);
                majorAlleleFragmentRanges = mutateFragments.subList(0, majorAlleleCount);
                minorAlleleFragmentRanges = mutateFragments.subList(majorAlleleCount, mutateFragments.size());
                this.writeIn(majorAlleleFragmentRanges, direct, readLength, seqErrorModel, mateFile1, mateFile2, this.exonSeq, "maj");
                majorAlleleFragmentRanges = null;
                this.writeIn(minorAlleleFragmentRanges, direct, readLength, seqErrorModel, mateFile1, mateFile2, mutExonSeq, "min");
                minorAlleleFragmentRanges = null;
            }

        } catch (Exception io) {
            io.printStackTrace();
        }
    }

    /**
     * m6a peak 的reads富集
     * @param mateFile1 测序文件
     * @param mateFile2 测序文件(只有在pair-end时候指定，否则为null)
     * @param readLength 测序read length
     * @param fragmentMean fragment长度均值
     * @param fragmentTheta fragment长度标准差
     * @param direct 测序方式 single-end 或 pair-end
     * @param mutExonSeq 当基因具有ASE位点时传入突变的序列，不存在时为null
     * @param refProp major allele的频率
     * @param seqErrorModel 测序误差模型
     */
    public void peakFragmentFromBackground(BufferedWriter mateFile1, BufferedWriter mateFile2, int readLength, int fragmentMean,
                                           int fragmentTheta, double refProp, String mutExonSeq, String direct,
                                           SequencingError seqErrorModel, Set<Integer> geneMutateSites,
                                           HashMap<Integer, Double> geneM6aAsm, HashMap<Integer, Boolean> minorBias) {
        List<Integer> mutations = null;
        if (geneMutateSites != null) {
            mutations = new ArrayList<>(geneMutateSites);
            Collections.sort(mutations);
        }
        ArrayList<Integer> m6aSites = null;
        if (geneM6aAsm != null) {
            m6aSites = new ArrayList<>(geneM6aAsm.keySet());
            Collections.sort(m6aSites);
        }
        int peakReadsCount, peakFragmentCount;
        UniformRealDistribution urd = new UniformRealDistribution(this.pmRange[0], this.pmRange[1]);
        ArrayList<int[]> fragmentRanges;
        List<int[]> randomFragmentRanges = new ArrayList<>();
        try {
            for (Integer m6aSite: this.m6aSiteFragments.keySet()) {
                fragmentRanges = this.m6aSiteFragments.getOrDefault(m6aSite, null);
                if (fragmentRanges == null)
                    continue;
                peakFragmentCount = fragmentRanges.size();
                // 计算每个m6A位点富集的reads数目并从所有覆盖甲基化的fragment中抽取
                peakReadsCount = this.m6aReadsCount(urd.sample(), peakFragmentCount);
                for (int i = 0; i < peakReadsCount; i++) {
                    int[] range = fragmentRanges.get((int) (Math.random() * fragmentRanges.size()));
                    randomFragmentRanges.add(range);
                }

                // 检验随机抽取的富集到peak的Fragment对应read是否覆盖ASE位点，如果覆盖则需对其进行记录
                boolean cover;
                for (int[] fragmentRange: randomFragmentRanges) {
                    cover = false;
                    if (mutations != null)
                        cover = this.ifCoverMutateStie(fragmentRange[0], fragmentRange[1], readLength, mutations,
                                                       this.m6aSiteMutateFragments);
                    if (!cover)
                        this.m6aSiteNormalFragments.add(fragmentRange);
                }
                randomFragmentRanges.clear();
                // 生成未覆盖ASE位点的reads
                this.writeIn(this.m6aSiteNormalFragments, direct, readLength, seqErrorModel, mateFile1, mateFile2, this.exonSeq, "ref");
                // 生成覆盖ASE位点的reads
                ArrayList<int[]> mutateFragments;
                int count, majorAlleleCount;
                List<int[]> majorAlleleFragmentRanges, minorAlleleFragmentRanges;
                for (Integer mutateSite: this.m6aSiteMutateFragments.keySet()) {
                    mutateFragments = this.m6aSiteMutateFragments.get(mutateSite);
                    boolean bias = false;
                    // 查看突变位点是否被m6A修饰位点区间覆盖，区间范围[site-fragmentMean-2*theta, site+fragmentMean+2*theta]
                    if (m6aSites != null) {
                        int site = this.ifInM6aPeakRange(mutateSite, m6aSites, fragmentMean+2*fragmentTheta);
                        // 如果在区域范围内，则该突变位点的reads数目依据修饰位点的ASM ratio进行生成，并查看是否具有major allele偏向性
                        if (site != -1) {
                            refProp = geneM6aAsm.get(site);
                            bias = minorBias.get(site);
                        }
                    }
                    count = mutateFragments.size();
                    if (bias)
                        majorAlleleCount = (int) (count * refProp);
                    else
                        majorAlleleCount = count - (int) (count * refProp);
                    Collections.shuffle(mutateFragments);
                    majorAlleleFragmentRanges = mutateFragments.subList(0, majorAlleleCount);
                    minorAlleleFragmentRanges = mutateFragments.subList(majorAlleleCount, mutateFragments.size());
                    this.writeIn(majorAlleleFragmentRanges, direct, readLength, seqErrorModel, mateFile1, mateFile2, this.exonSeq, "maj");
                    majorAlleleFragmentRanges = null;
                    this.writeIn(minorAlleleFragmentRanges, direct, readLength, seqErrorModel, mateFile1, mateFile2, mutExonSeq, "min");
                    minorAlleleFragmentRanges = null;
                }
                this.m6aSiteNormalFragments.clear();
                this.m6aSiteMutateFragments.clear();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        this.m6aSiteFragments.clear();
    }

    /**
     * 检查fragment对应的reads是否覆盖m6A位点
     * @param start fragment起始位点
     * @param end fragment终止位点
     * @param sites 该基因m6A位点的集合
     */
    private void ifCoverM6aSite(int start, int end, List<Integer> sites) {
        ArrayList<int[]> fragmentRange;
        for (Integer site: sites) {
            if (start <= site && site <= end) {
                fragmentRange = this.m6aSiteFragments.getOrDefault(site, new ArrayList<>());
                fragmentRange.add(new int[]{start, end});
                this.m6aSiteFragments.put(site, fragmentRange);
                break;
            }
        }
    }

    /**
     * 检查fragment对应的read是否覆盖突变位点
     * @param start fragment起始位点
     * @param end fragment终止位点
     * @param sites 该基因突变位点集合
     * @return 是否覆盖
     */
    private boolean ifCoverMutateStie(int start, int end, int readLength, List<Integer> sites,
                                      HashMap<Integer, ArrayList<int[]>> fragments) {
        ArrayList<int[]> fragmentRange;
        int readEnd = ((start + readLength -1) > this.exonSeq.length())? end: (start + readLength);
        for (Integer site: sites) {
            if (start <= site && site <= readEnd) {
                fragmentRange = fragments.getOrDefault(site, new ArrayList<>());
                fragmentRange.add(new int[]{start, end});
                fragments.put(site, fragmentRange);
                return true;
            }
        }
        return false;
    }

    /**
     * 判断reads上的突变位点是否位于甲基化peak内
     * @param mutatePosition 突变位点位置
     * @param m6aSites 甲基化位点集合
     * @param fragmentLength fragment长度
     * @return 是否在peak范围内
     */
    private int ifInM6aPeakRange(int mutatePosition, ArrayList<Integer> m6aSites, int fragmentLength) {
        for (Integer site: m6aSites) {
            if (Math.max(0, site-fragmentLength) <= mutatePosition && mutatePosition <= site + fragmentLength)
                return site;
        }
        return -1;
    }

    /**
     * get break and end point of fragment on exon sequence, according to the random generated position
     * @param fragmentLength fragment length
     * @return break and end point
     */
    private int[] getBreakEndPoint(int fragmentLength) {
        int break_point, end_point;
        int exonSeqLength = this.exonSeq.length();
        if (fragmentLength > exonSeqLength) {
            break_point = 0;
            end_point = exonSeqLength;
        } else {
            UniformIntegerDistribution uniformdi = new UniformIntegerDistribution(0, exonSeqLength - fragmentLength);
            break_point = uniformdi.sample();
            end_point = break_point + fragmentLength;
            uniformdi = null;
        }

        return new int[]{break_point, end_point};
    }

    private void writeIn(List<int[]> fragmentRanges, String direct, int readLength, SequencingError seqErrorModel,
                         BufferedWriter mateFile1, BufferedWriter mateFile2, String exonSequence, String type) {
        String fragmentString;
        Fragmentation fragment;
        int break_point, endPoint;
        // 生成一般fragment的reads
        for (int[] fragmentRange: fragmentRanges) {
            // 获取fragment的起始、终止位点
            break_point = fragmentRange[0];
            endPoint = break_point + readLength - 1;
            // 如果基因含有ASE位点，则随机选取refCount个fragment作为major allele，其余的是minor allele
            if (endPoint > exonSequence.length())
                endPoint = exonSequence.length() - 1;
            fragmentString = exonSequence.substring(break_point, endPoint+1);
            fragment = this.getSequencingReads(fragmentString, break_point, endPoint, readLength, direct);
            if (direct.equals("SE"))
                this.writeReadInFile(fragment, seqErrorModel, mateFile1, readLength, break_point, endPoint, type);
            else
                this.pairReadToFile(fragment, seqErrorModel, mateFile1, mateFile2, readLength, break_point, endPoint);
            fragment = null;
            fragmentString = null;
        }
    }

    /**
     * 单端测序的read写入文件
     * @param fragment Fragmentation 对象
     * @param seqErrorModel sequencing error model对象
     * @param fw BufferedWriter
     * @param readLength read长度
     * @param break_point fragment 在外显子序列上的起始位点
     * @param end_point fragment 在外显子序列上的终止位点
     */
    private void writeReadInFile(Fragmentation fragment, SequencingError seqErrorModel, BufferedWriter fw,
                                 int readLength, int break_point, int end_point, String type) {
        try {
            String sequencingRead = fragment.getSingleEndRead();
            sequencingRead = seqErrorModel.pcrErrorReads(sequencingRead);
            if (sequencingRead.length() != readLength) {
                int baseNNum = readLength - sequencingRead.length();
                String baseNSeq = this.fillBaseN(baseNNum);
                sequencingRead = sequencingRead + baseNSeq;
            }

            fw.write(">chr" + this.chr + "_" + this.geneId + "_" +break_point + "_" + end_point + " " + type);
            fw.newLine();
            fw.write(sequencingRead);
            fw.newLine();
        } catch (IOException | StringIndexOutOfBoundsException e) {
            e.printStackTrace();
        }
    }

    /**
     * 双端测序的read写入文件
     * @param fragment Fragmentation 对象
     * @param seqErrorModel sequencing error model对象
     * @param mateFile1 输出文件1
     * @param mateFile2 输出文件2
     * @param readLength read长度
     * @param break_point fragment 在外显子序列上的起始位点
     * @param end_point fragment 在外显子序列上的终止位点
     */
    private void pairReadToFile(Fragmentation fragment, SequencingError seqErrorModel, BufferedWriter mateFile1,
                                BufferedWriter mateFile2, int readLength, int break_point, int end_point) {
        try {
            String[] sequencingRead = fragment.getPairEndRead();
            String mate1 = sequencingRead[0];
            String mate2 = sequencingRead[1];
            mate1 = seqErrorModel.pcrErrorReads(mate1);
            mate2 = seqErrorModel.pcrErrorReads(mate2);
            if (mate1.length() != readLength) {
                int baseNNum = readLength - mate1.length();
                String baseNSeq = this.fillBaseN(baseNNum);
                mate1 = mate1 + baseNSeq;
                mate2 = baseNSeq + mate2;
            }

            mateFile1.write(">chr" + this.chr + "_" + this.geneId + "_" +break_point + "_" + end_point);
            mateFile1.newLine();
            mateFile1.write(mate1);
            mateFile1.newLine();

            mateFile2.write(">chr" + this.chr + "_" + this.geneId + "_"+break_point + "_" + end_point);
            mateFile2.newLine();
            mateFile2.write(mate2);
            mateFile2.newLine();
        } catch (IOException | StringIndexOutOfBoundsException e) {
            e.printStackTrace();
        }
    }

    /**
     * use base N to fill the reads, when read length less than 50
     * @param baseNNum number of base N to fill
     * @return base N sequence
     */
    private String fillBaseN(int baseNNum) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < baseNNum; i++) {
            sb.append("N");
        }
        return sb.toString();
    }

    /**
     * generate fragment and form sequencing reads
     * @param break_point start position of the fragment on exon sequence
     * @param end_point end position of the fragment on exon sequence
     * @param readLength sequencing read length
     * @param direction "SE" single-end, "PE" pair-end
     * @return fragment instance
     */
    private Fragmentation getSequencingReads(String fragmentString, int break_point, int end_point, int readLength, String direction) {
        Fragmentation fragment = new Fragmentation(fragmentString, break_point, end_point);
        double randNum = Math.random();
        String strand = (randNum < 0.5)? "+" : "-";
        if (direction.equals("SE")) {   // single-end sequencing
            fragment.singleEndSequencingRead(strand, readLength);
        } else {    // pair-end sequencing
            fragment.pairEndSequencingRead(strand, readLength);
        }

        return fragment;
    }

    /**
     * 释放内存
     */
    public void release() {
        this.ipFragment.clear();
        this.inputFragment.clear();
        this.inputMutateFragments.clear();
        this.ipMutateFragments.clear();
    }
}

