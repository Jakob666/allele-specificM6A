package AseSeqSimulator;

import GTFComponent.ElementRecord;
import GTFComponent.TranscriptRecord;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.UniformIntegerDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * random fragment a gene transcript sequence, the fragment length obeys normal distribution
 */
public class Gene {
    private String geneId, strand, geneName;
    private TranscriptRecord longestTranscriptRecord;
    private double RPKM, pm;
    // exonSeq是基因最长的一个转录本的外显子序列，exomeSeq是基因全外显子组序列
    private String chr, exonSeq, exomeSeq = null;
    private int readsCount;
    private int dnaReadsCount;
    private double[] pmRange = new double[]{0.05, 0.9};
    // exonList是基因最长的一个转录本的外显子的区间记录，exomeList是基因全外显子组区间记录(已去除重叠区域)
    private ElementRecord exonList = null;
    private ArrayList<int[]> exomeList;
    // m6aSiteFragments记录IP样本中Gene在每个m6A位点上覆盖该位点的fragment
    private HashMap<Integer, ArrayList<int[]>> m6aSiteFragments = new HashMap<>(), m6aSiteMutateFragments = new HashMap<>();
    // 对应的read覆盖SNP位点的fragment
    private HashMap<Integer, ArrayList<int[]>> inputMutateFragments = new HashMap<>(), ipMutateFragments = new HashMap<>();
    // inputFragment和ipFragment记录没有覆盖SNP位点的fragment，
    private ArrayList<int[]> inputFragment = new ArrayList<>(), ipFragment = new ArrayList<>(), m6aSiteNormalFragments = new ArrayList<>();
    private ArrayList<int[]> exomeFragment = new ArrayList<>();
    private HashMap<Integer, ArrayList<int[]>> exomeMutateFragments = new HashMap<>();
    private Object obj;

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

    public void setObj(Object obj) {
        this.obj = obj;
    }

    public void setPmRange(double lower, double upper) {
        this.pmRange = new double[]{lower, upper};
    }

    public void setRPKM(double RPKM) {
        this.RPKM = RPKM;
    }

    public void setExomeSeq(String exomeSeq) {
        this.exomeSeq = exomeSeq;
    }

    public void setExomeList(ArrayList<int[]> exomeList) {
        this.exomeList = exomeList;
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

    public int getReadsCount() {
        return this.readsCount;
    }

    public ElementRecord getExonList() {
        return this.exonList;
    }

    public ArrayList<int[]> getExomeList() {
        return this.exomeList;
    }

    public String getExonSeq() {
        return this.exonSeq;
    }

    public String getExomeSeq() {
        return this.exomeSeq;
    }

    public TranscriptRecord getLongestTranscriptRecord() {
        return this.longestTranscriptRecord;
    }

    /**
     * 将基因最长转录本所对应的的外显子区域拼合得到完整的外显子序列
     * @param twoBit TwoBitParser instance
     */
    public void splicing(TwoBitParser twoBit) {
        HashMap<String, ElementRecord> elementsOnTranscript = this.longestTranscriptRecord.getElementList();
        String strand = this.longestTranscriptRecord.getStrand();
        ElementRecord transcriptExon = elementsOnTranscript.getOrDefault("exon", null);
        StringBuilder sb = new StringBuilder();
        if (transcriptExon != null)
            this.exonList = transcriptExon;
        try {
            while (transcriptExon != null) {
                int exonStart = transcriptExon.getElementStart();
                int exonEnd = transcriptExon.getElementEnd();
                twoBit.reset();
                String exonSeq = twoBit.loadFragment(exonStart-1, (exonEnd - exonStart+1));
                if (strand.equals("-")) {
                    sb.insert(0, CommonMethod.AntiChain(exonSeq));
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
     * 依据RPKM值和文库大小计算该基因的reads count，RNA-seq的reads count
     * @param librarySize library size
     */
    public void calculateReadsCountViaLibrarySize(long librarySize) {
        this.readsCount = (int) ((this.exonSeq.length() * librarySize * this.RPKM) / Math.pow(10, 9));
    }

    /**
     * 通过设定的测序深度和文库大小计算该基因的reads count，如果设置了depth相当于DNA-seq测序
     * @param depth sequencing depth
     * @param readLength read length
     */
    public void calculateReadsCountViaSequencingDepth(int depth, int readLength) {
        this.dnaReadsCount = this.exomeSeq.length() * depth / readLength;
    }

    /**
     * 依据m6A的富集率计算IP样本的reads count
     * @param backgroundSize 背景大小
     */
    private int m6aReadsCount(int backgroundSize) {
        return  (int) (backgroundSize * (pm / (1 - pm)));
    }

    /**
     * 设置基因的甲基化富集度
     */
    public void setGenePM() {
        UniformRealDistribution urd = new UniformRealDistribution(0.85, 0.95);
        this.pm = urd.sample();
    }

    public double getPm() {
        return this.pm;
    }

    /**
     * 随机抽取reads count条fragment，将其富集到INPUT组
     * @param fragmentMean fragment平均长度
     * @param fragmentTheta fragment长度标准差
     */
    public void enrichInputFragment(int fragmentMean, int fragmentTheta, int readLength, Set<Integer> mutatePositions) {
        List<Integer> mutations = null;
        if (mutatePositions != null) {
            mutations = new ArrayList<>(mutatePositions);
            Collections.sort(mutations);
        }
        // 用于随机生成fragment的长度
        NormalDistribution nordi = new NormalDistribution(fragmentMean, fragmentTheta);
        int curReadsCount = 0, fragmentLength, breakPoint, endPoint;
        boolean cover;
        while (curReadsCount < this.readsCount) {
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
        mutations = null;
        nordi = null;
    }

    /**
     * 随机抽取reads count条fragment，将其富集到IP组
     * @param fragmentMean fragment平均长度
     * @param fragmentTheta fragment长度标准差
     * @param geneM6aSites 该基因上m6A修饰位点
     */
    public void enrichIpFragment(int fragmentMean, int fragmentTheta, int readLength, Set<Integer> geneM6aSites,
                                 Set<Integer> mutatePositions) {
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
        while (curReadsCount < this.readsCount) {
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
        m6aSites = null;
        mutations = null;
        nordi = null;
    }

    /**
     * 生成测序样本
     * @param mateFile1 测序文件
     * @param mateFile2 测序文件(只有在pair-end时候指定，否则为null)
     * @param readLength 测序read length
     * @param multiple 实验重复次数
     * @param direct 测序方式 single-end 或 pair-end
     * @param mutExonSeq 当基因具有ASE位点时传入突变的序列，不存在时为null
     * @param refProp major allele的频率
     * @param seqErrorModel 测序误差模型
     * @param sample ip或input
     */
    public void generateReads(BufferedWriter mateFile1, BufferedWriter mateFile2, int readLength, int multiple,
                                   double refProp, String mutExonSeq, String direct, SequencingError seqErrorModel,
                                   String sample, ExecutorService service) {
        ArrayList<int[]> fragmentRanges;
        HashMap<Integer, ArrayList<int[]>> mutateFragmentRanges;
        fragmentRanges = (sample.equals("ip"))? this.ipFragment : this.inputFragment;
        mutateFragmentRanges = (sample.equals("ip"))? this.ipMutateFragments : this.inputMutateFragments;

        try {
            // 生成一般fragment的reads(reads不覆盖SNP位点)
            Runnable runNormal = this.writeIn(fragmentRanges, direct, readLength, seqErrorModel, mateFile1, mateFile2, this.exonSeq, "ref", this.chr, this.geneId, this.obj);
            service.execute(runNormal);
            // 生成覆盖突变位点的fragment的reads
            ArrayList<int[]> mutateFragments;
            int count, majorAlleleCount;
            List<int[]> majorAlleleFragmentRanges, minorAlleleFragmentRanges;
            Runnable runMajor, runMinor;
            for (Integer mutateSite: mutateFragmentRanges.keySet()) {
                // 获取每个包含ASE位点的fragment的范围
                mutateFragments = mutateFragmentRanges.get(mutateSite);
                count = mutateFragments.size();
                majorAlleleCount = (int) (count * refProp);
                Collections.shuffle(mutateFragments);
                majorAlleleFragmentRanges = mutateFragments.subList(0, majorAlleleCount);
                minorAlleleFragmentRanges = mutateFragments.subList(majorAlleleCount, mutateFragments.size());
                runMajor = this.writeIn(majorAlleleFragmentRanges, direct, readLength, seqErrorModel, mateFile1, mateFile2, this.exonSeq, "major", this.chr, this.geneId, this.obj);
                service.execute(runMajor);
                majorAlleleFragmentRanges = null;
                runMinor = this.writeIn(minorAlleleFragmentRanges, direct, readLength, seqErrorModel, mateFile1, mateFile2, mutExonSeq,"minor", this.chr, this.geneId, this.obj);
                service.execute(runMinor);
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
     * @param multiple 实验重复次数
     * @param direct 测序方式 single-end 或 pair-end
     * @param mutExonSeq 当基因具有ASE位点时传入突变的序列，不存在时为null
     * @param refProp major allele的频率
     * @param seqErrorModel 测序误差模型
     */
    public void peakFragmentFromBackground(BufferedWriter mateFile1, BufferedWriter mateFile2, int readLength, int multiple,
                                           double refProp, String mutExonSeq, String direct, SequencingError seqErrorModel,
                                           Set<Integer> geneMutateSites, ExecutorService service) {
        List<Integer> mutations = null;
        if (geneMutateSites != null) {
            mutations = new ArrayList<>(geneMutateSites);
            Collections.sort(mutations);
        }
        int peakReadsCount, peakFragmentCount;
        // 模拟数据将甲基化的富集度暂时设定一个较高的值
//        UniformRealDistribution urd = new UniformRealDistribution(this.pmRange[0], this.pmRange[1]);

        ArrayList<int[]> fragmentRanges;
        List<int[]> randomFragmentRanges = new ArrayList<>();
        try {
            for (Integer m6aSite: this.m6aSiteFragments.keySet()) {
                fragmentRanges = this.m6aSiteFragments.getOrDefault(m6aSite, null);
                if (fragmentRanges == null)
                    continue;
                peakFragmentCount = fragmentRanges.size();
                // 计算每个m6A位点富集的reads数目并从所有覆盖甲基化的fragment中抽取
                peakReadsCount = this.m6aReadsCount(peakFragmentCount);
                for (int i = 0; i < peakReadsCount; i++) {
                    Collections.shuffle(fragmentRanges);
                    int[] range = fragmentRanges.get(0);
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
                Runnable runNormal = this.writeIn(this.m6aSiteNormalFragments, direct, readLength, seqErrorModel, mateFile1,
                                                  mateFile2, this.exonSeq, "ref", this.chr, this.geneId, this.obj);
                service.execute(runNormal);
                // 生成覆盖ASE位点的reads
                ArrayList<int[]> mutateFragments;
                int count, majorAlleleCount;
                List<int[]> majorAlleleFragmentRanges, minorAlleleFragmentRanges;
                Runnable runMajor, runMinor;
                for (Integer mutateSite: this.m6aSiteMutateFragments.keySet()) {
                    mutateFragments = this.m6aSiteMutateFragments.get(mutateSite);
                    count = mutateFragments.size();
                    majorAlleleCount = (int) (count * refProp);
                    Collections.shuffle(mutateFragments);
                    majorAlleleFragmentRanges = mutateFragments.subList(0, majorAlleleCount);
                    minorAlleleFragmentRanges = mutateFragments.subList(majorAlleleCount, mutateFragments.size());
                    runMajor = this.writeIn(majorAlleleFragmentRanges, direct, readLength, seqErrorModel, mateFile1, mateFile2,
                                            this.exonSeq, "major", this.chr, this.geneId, this.obj);
                    service.execute(runMajor);
                    majorAlleleFragmentRanges = null;
                    runMinor = this.writeIn(minorAlleleFragmentRanges, direct, readLength, seqErrorModel, mateFile1, mateFile2,
                                            mutExonSeq,"minor", this.chr, this.geneId, this.obj);
                    service.execute(runMinor);
                    minorAlleleFragmentRanges = null;
                }
                this.m6aSiteNormalFragments.clear();
                this.m6aSiteMutateFragments.clear();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        this.m6aSiteFragments = null;
    }

    /**
     * 随机抽取reads count条fragment，用于生成 exome sequencing数据
     * @param fragmentMean fragment平均长度
     * @param fragmentTheta fragment长度标准差
     */
    public void enrichExomeFragment(int fragmentMean, int fragmentTheta, int readLength, Set<Integer> mutatePositions) {
        List<Integer> mutations = null;
        if (mutatePositions != null) {
            mutations = new ArrayList<>(mutatePositions);
            Collections.sort(mutations);
        }
        // 用于随机生成fragment的长度
        NormalDistribution nordi = new NormalDistribution(fragmentMean, fragmentTheta);
        int curReadsCount = 0, fragmentLength, breakPoint, endPoint;
        boolean cover;
        while (curReadsCount < this.dnaReadsCount) {
            cover = false;
            // 随机生成fragment的长度并确定fragment在外显子序列上的起始终止位点
            fragmentLength = Math.abs((int) nordi.sample());
            int[] breakAndEndPoint = this.getBreakEndPoint(fragmentLength);
            breakPoint = breakAndEndPoint[0];
            endPoint = breakAndEndPoint[1];

            // 检验fragment对应的read是否覆盖ASE位点，如果覆盖则需对其进行记录
            if (mutations != null)
                cover = this.ifCoverMutateStie(breakPoint, endPoint, readLength, mutations, this.exomeMutateFragments);
            if (!cover)
                this.exomeFragment.add(breakAndEndPoint);
            curReadsCount++;
        }
        mutations = null;
        nordi = null;
    }

    /**
     * 生成外显子组测序样本
     * @param mateFile1 测序文件
     * @param mateFile2 测序文件(只有在pair-end时候指定，否则为null)
     * @param readLength 测序read length
     * @param multiple 实验重复次数
     * @param direct 测序方式 single-end 或 pair-end
     * @param mutExomeSeq 当基因具有ASE位点时传入突变的序列，不存在时为null
     * @param refProp major allele的频率
     * @param seqErrorModel 测序误差模型
     */
    public void generateExomeSeqReads(BufferedWriter mateFile1, BufferedWriter mateFile2, int readLength, int multiple,
                                      double refProp, String mutExomeSeq, String direct, SequencingError seqErrorModel,
                                      ExecutorService service) {
        try {
            // 生成一般fragment的reads(reads不覆盖SNP位点)
            Runnable runNormal = this.writeIn(this.exomeFragment, direct, readLength, seqErrorModel, mateFile1, mateFile2,
                                              this.exomeSeq, "ref", this.chr, this.geneId, this.obj);
            service.execute(runNormal);

            // 生成覆盖突变位点的fragment的reads
            ArrayList<int[]> mutateFragments;
            int count, majorAlleleCount;
            List<int[]> majorAlleleFragmentRanges, minorAlleleFragmentRanges;
            Runnable runMajor, runMinor;
            for (Integer mutateSite: this.exomeMutateFragments.keySet()) {
                // 获取每个包含ASE位点的fragment的范围
                mutateFragments = this.exomeMutateFragments.get(mutateSite);
                count = mutateFragments.size();
                majorAlleleCount = (int) (count * refProp);
                Collections.shuffle(mutateFragments);
                majorAlleleFragmentRanges = mutateFragments.subList(0, majorAlleleCount);
                minorAlleleFragmentRanges = mutateFragments.subList(majorAlleleCount, mutateFragments.size());
                runMajor = this.writeIn(majorAlleleFragmentRanges, direct, readLength, seqErrorModel,
                                        mateFile1, mateFile2, this.exomeSeq, "major", this.chr,
                                        this.geneId, obj);
                service.execute(runMajor);
                majorAlleleFragmentRanges = null;
                runMinor = this.writeIn(minorAlleleFragmentRanges, direct, readLength, seqErrorModel,
                                        mateFile1, mateFile2, mutExomeSeq,"minor", this.chr, this.geneId, this.obj);
                service.execute(runMinor);
                minorAlleleFragmentRanges = null;
            }
        } catch (Exception io) {
            io.printStackTrace();
        }
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

    private Runnable writeIn(List<int[]> fragmentRanges, String direct, int readLength, SequencingError seqErrorModel,
                             BufferedWriter mateFile1, BufferedWriter mateFile2, String exonSequence,  String type,
                             String chr, String geneId, Object obj) {
        // 此处将List转化为Vector因为Vector是线程安全的
        Vector<int[]> ranges = new Vector<>(fragmentRanges);
        fragmentRanges = null;
        return new Runnable() {
            @Override
            public void run() {
                String fragmentString;
                Fragmentation fragment;
                int break_point, endPoint;
                // 生成一般fragment的reads
                for (int[] fragmentRange: ranges) {
                    synchronized (obj) {
                        // 获取fragment的起始、终止位点
                        break_point = fragmentRange[0];
                        endPoint = fragmentRange[1];
                        // 如果基因含有ASE位点，则随机选取refCount个fragment作为major allele，其余的是minor allele
                        fragmentString = exonSequence.substring(break_point, endPoint);
                        fragment = getSequencingReads(fragmentString, readLength, direct);
                        if (direct.equals("SE"))
                            writeReadInFile(fragment, seqErrorModel, mateFile1, readLength, break_point, endPoint, type,
                                            chr, geneId);
                        else
                            pairReadToFile(fragment, seqErrorModel, mateFile1, mateFile2, readLength, break_point,
                                           endPoint, type, chr, geneId);
                        fragment = null;
                        fragmentString = null;
                    }
                }
            }
        };

    }

    /**
     * 单端测序的read写入文件
     * @param fragment Fragmentation 对象
     * @param seqErrorModel sequencing error model对象
     * @param fw BufferedWriter
     * @param readLength read长度
     * @param break_point fragment 在外显子序列上的起始位点
     * @param end_point fragment 在外显子序列上的终止位点
     * @param type major or minor
     */
    private static void writeReadInFile(Fragmentation fragment, SequencingError seqErrorModel, BufferedWriter fw,
                                        int readLength, int break_point, int end_point, String type, String chr,
                                        String geneId) {
        try {
            String sequencingRead = fragment.getSingleEndRead();
            String strand = fragment.getReadStrand();
            sequencingRead = seqErrorModel.pcrErrorReads(sequencingRead);
            if (sequencingRead.length() != readLength) {
                int baseNNum = readLength - sequencingRead.length();
                String baseNSeq = fillBaseN(baseNNum);
                sequencingRead = sequencingRead + baseNSeq;
            }

            fw.write(">chr" + chr + "_" + geneId + "_" +type + "_"+break_point + ":" + end_point +"\t" + strand); // +"_"+readsStart+":"+readsEnd
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
     * @param type major or minor
     */
    private static void pairReadToFile(Fragmentation fragment, SequencingError seqErrorModel, BufferedWriter mateFile1,
                                BufferedWriter mateFile2, int readLength, int break_point, int end_point, String type,
                                       String chr, String geneId) {
        try {
            String[] sequencingRead = fragment.getPairEndRead();
            String mate1 = sequencingRead[0];
            String mate2 = sequencingRead[1];
            String strand = fragment.getReadStrand();
            mate1 = seqErrorModel.pcrErrorReads(mate1);
            mate2 = seqErrorModel.pcrErrorReads(mate2);
            if (mate1.length() != readLength) {
                int baseNNum = readLength - mate1.length();
                String baseNSeq = fillBaseN(baseNNum);
                mate1 = mate1 + baseNSeq;
                mate2 = baseNSeq + mate2;
            }

            String label = ">chr" + chr + "_" + geneId + "_" +type + "_"+break_point + ":" + end_point +"\t" + strand;
            mateFile1.write(label);
            mateFile1.newLine();
            mateFile1.write(mate1);
            mateFile1.newLine();

            mateFile2.write(label);
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
    private static String fillBaseN(int baseNNum) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < baseNNum; i++) {
            sb.append("N");
        }
        return sb.toString();
    }

    /**
     * generate fragment and form sequencing reads
     * @param readLength sequencing read length
     * @param direction "SE" single-end, "PE" pair-end
     * @return fragment instance
     */
    private static Fragmentation getSequencingReads(String fragmentString, int readLength, String direction) {
        Fragmentation fragment = new Fragmentation(fragmentString);
        double randNum = Math.random();
        if (direction.equals("SE")) {   // single-end sequencing
            fragment.singleEndSequencingRead(randNum, readLength);
        } else {    // pair-end sequencing
            fragment.pairEndSequencingRead(randNum, readLength);
        }

        return fragment;
    }

    /**
     * 释放内存
     */
    public void release() {
        this.ipFragment = null;
        this.inputFragment = null;
        this.exomeFragment = null;
        this.inputMutateFragments = null;
        this.ipMutateFragments = null;
        this.exomeMutateFragments = null;
    }
}

