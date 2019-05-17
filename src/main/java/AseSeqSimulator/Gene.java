package AseSeqSimulator;

import GTFComponent.ElementRecord;
import GTFComponent.TranscriptRecord;
import PeakSimulator.Peak;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.UniformIntegerDistribution;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;

/**
 * random fragment a gene transcript sequence, the fragment length obeys normal distribution
 */
public class Gene {
    private String geneId, strand, geneName;
    private TranscriptRecord longestTranscriptRecord;
    private int geneStart, geneEnd;
    private double RPKM;
    private String chr, exonSeq;
    private int readsCount, refReadsCount = 0, altReadsCount = 0;
    private double[] pmRange;
    private ElementRecord exonList = null;
    private LinkedList<Peak> peakList;
    public Gene(String geneId, int geneStart, String geneName, int geneEnd, String strand, String chr) {
        this.geneId = geneId;
        this.geneStart = geneStart;
        this.geneEnd = geneEnd;
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
        pmRange = new double[]{lower, upper};
    }

    public void setRPKM(double RPKM) {
        this.RPKM = RPKM;
    }

    public void setPeakList(LinkedList<Peak> peakList) {
        this.peakList = peakList;
    }

    public int getGeneStart() {
        return this.geneStart;
    }

    public int getGeneEnd() {
        return this.geneEnd;
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

    public int getRefReadsCount() {
        return this.refReadsCount;
    }

    public int getAltReadsCount() {
        return this.altReadsCount;
    }

    public double getRPKM() {
        return this.RPKM;
    }

    public double[] getPmRange() {
        return this.pmRange;
    }

    public ElementRecord getExonList() {
        return this.exonList;
    }

    public LinkedList<Peak> getPeakList() {
        return this.peakList;
    }

    public String getExonSeq() {
        return this.exonSeq;
    }

    public TranscriptRecord getLongestTranscriptRecord() {
        return this.longestTranscriptRecord;
    }

    /**
     * join all exon region of a gene together and form the exon sequence
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
     * calculate sequencing reads count base on library size
     * @param librarySize library size
     */
    public void calculateReadsCountViaLibrarySize(long librarySize) {
        this.readsCount = (int) ((this.exonSeq.length() * librarySize * this.RPKM) / Math.pow(10, 9));
    }

    /**
     * calculate sequencing reads count base on sequencing depth
     * @param depth sequencing depth
     * @param readLength read length
     */
    public void calculateReadsCountViaSequencingDepth(int depth, int readLength, long librarySize) {
        this.readsCount = this.exonSeq.length() * depth / readLength;
        this.RPKM = Math.pow(10.0, 9) * this.readsCount / this.exonSeq.length() / librarySize;
    }

    /**
     * generate reads in exonic regions and form fragment sequence, the length of which obeys a normal distribution
     * and the reads count is calculated based on the library size exonic region length and RPKM value
     * @param fragmentMean the mean length of fragments
     * @param fragmentTheta fragment length std
     * @param mateFile1 BufferedWriter instance
     * @param mateFile2 BufferedWriter instance
     * @param readLength sequencing read length
     * @param multiple replication of the experiment
     */
    public void generateInputReads(int fragmentMean, int fragmentTheta, BufferedWriter mateFile1, BufferedWriter mateFile2,
                                   int readLength, int multiple, double refProp, String mutExonSeq, String direct,
                                   SequencingError seqErrorModel) {
        // 位置相同的reads的数目
        HashMap<String, Integer> readsDistribution = new HashMap<>();
        // 位置相同的reads的fragment在exon上面起始、终止位点
        HashMap<String, ArrayList<Integer>> fragmentEndSites = new HashMap<>();
        // used to randomly generate fragment length
        NormalDistribution nordi = new NormalDistribution(fragmentMean, fragmentTheta);
        int curReadsCount = 0;
        String refFragmentString, altFragmentString;
        Fragmentation refFragment, altFragment;
        try {
            int break_point, end_point, fragmentLength;
            int readStart, readEnd;
            String readRecord;

            while (curReadsCount < this.readsCount) {
                // 随机生成fragment的长度
                fragmentLength = Math.abs((int) nordi.sample());
                // 通过fragment长度确定确定fragment在外显子序列上的起始终止位点
                int[] breakAndEndPoint = this.getBreakEndPoint(fragmentLength);
                break_point = breakAndEndPoint[0];
                end_point = breakAndEndPoint[1];

                readStart = break_point;
                readEnd = (end_point - break_point >= readLength)? (break_point+readLength):end_point;
                readRecord = readStart + "-" + readEnd;

                Integer countRecord = readsDistribution.getOrDefault(readRecord, 0);
                ArrayList<Integer> fragmentRecord = fragmentEndSites.getOrDefault(readRecord, new ArrayList<>());
                fragmentRecord.add(end_point);
                fragmentEndSites.put(readRecord, fragmentRecord);
                readsDistribution.put(readRecord, countRecord+1);
                curReadsCount++;
            }

            for (String readRange: readsDistribution.keySet()) {
                String[] points = readRange.split("-");
                // 获取reads的起始、终止位点
                break_point = Integer.parseInt(points[0]);
                readEnd = Integer.parseInt(points[1]);
                // 获取reads对应的fragment的终止位点
                ArrayList<Integer> endPoints = fragmentEndSites.get(readRange);
                Collections.shuffle(endPoints);

                // 相同位置的reads数目，计算存在ASE时，major allele的reads数目
                int count = readsDistribution.get(readRange);
                int refCount = (int) (count * refProp);

                // 如果基因含有ASE位点
                if (mutExonSeq != null) {
                    List<Integer> majorAlleleFragmentEnds = endPoints.subList(0, refCount);
                    List<Integer> minorAlleleFragmentEnds = endPoints.subList(refCount, endPoints.size());
                    // 生成major allele的reads
                    for (Integer endPoint: majorAlleleFragmentEnds) {
                        refFragmentString = this.exonSeq.substring(break_point, endPoint);
                        refFragment = this.getSequencingReads(refFragmentString, break_point, endPoint, readLength, direct);
                        if (direct.equals("SE"))
                            this.writeReadInFile(refFragment, seqErrorModel, mateFile1, readLength, break_point, endPoint, "ref");
                        else
                            this.pairReadToFile(refFragment, seqErrorModel, mateFile1, mateFile2, readLength, break_point, endPoint, "ref");
                        refFragment = null;
                        refFragmentString = null;
                    }
                    // 生成minor allele的reads
                    for (Integer endPoint: minorAlleleFragmentEnds) {
                        altFragmentString = mutExonSeq.substring(break_point, endPoint);
                        altFragment = this.getSequencingReads(altFragmentString, break_point, endPoint, readLength, direct);
                        if (direct.equals("SE"))
                            this.writeReadInFile(altFragment, seqErrorModel, mateFile1, readLength, break_point, endPoint, "ref");
                        else
                            this.pairReadToFile(altFragment, seqErrorModel, mateFile1, mateFile2, readLength, break_point, endPoint, "ref");
                        altFragment = null;
                        altFragmentString = null;
                    }
                } else { // 如果基因不包含SNP位点
                    for (Integer endPoint: endPoints) {
                        refFragmentString = this.exonSeq.substring(break_point, endPoint);
                        refFragment = this.getSequencingReads(refFragmentString, break_point, endPoint, readLength, direct);
                        if (direct.equals("SE"))
                            this.writeReadInFile(refFragment, seqErrorModel, mateFile1, readLength, break_point, endPoint, "ref");
                        else
                            this.pairReadToFile(refFragment, seqErrorModel, mateFile1, mateFile2, readLength, break_point, endPoint, "ref");

                        // release memory
                        refFragment = null;
                        refFragmentString = null;
                    }
                }
            }

        } catch (Exception io) {
            io.printStackTrace();
        }
    }

    /**
     * if a sequencing reads covers SNP site on exon sequence
     * @param start start position
     * @param end end position
     * @param sites set of sites
     * @return boolean
     */
    private boolean ifCoverSite(int start, int end, Set<Integer> sites) {
        for (Integer site: sites) {
            if (start <= site && site <= end)
                return true;
        }
        return false;
    }

    public void generateIpReads(int fragmentMean, int fragmentTheta, BufferedWriter mateFile1, BufferedWriter mateFile2,
                                int readLength, Set<Integer> geneM6aSites, int multiple, double refProp, String mutExonSeq,
                                String direct, SequencingError seqErrorModel) {
        // 如果基因外显子区域不存在m6A修饰位点，则跳过该基因
        if (geneM6aSites.size() == 0)
            return;
        // 记录位置相同的reads的数目
        HashMap<String, Integer> readsDistribution = new HashMap<>();
        // 位置相同的reads的fragment在exon上面起始、终止位点
        HashMap<String, ArrayList<Integer>> fragmentEndSites = new HashMap<>();
        // 用于随机生成fragment的长度
        NormalDistribution nordi = new NormalDistribution(fragmentMean, fragmentTheta);
        int curReadsCount = 0;
        String refFragmentString, altFragmentString;
        Fragmentation refFragment, altFragment;

        try {
            int break_point, end_point, fragmentLength;
            int readStart, readEnd;
            String readRecord;
            boolean cover;
            while (curReadsCount < this.readsCount) {
                // 随机生成一个fragment长度
                fragmentLength = Math.abs((int) nordi.sample());
                // 根据生成的fragment长度得到该fragment在外显子序列上的起始、终止位点
                int[] breakAndEndPoint = this.getBreakEndPoint(fragmentLength);
                break_point = breakAndEndPoint[0];
                end_point = breakAndEndPoint[1];
                // 检查这段区域是否覆盖m6A位点，如果没有覆盖则跳过该fragment
                cover = this.ifCoverSite(break_point, end_point, geneM6aSites);
                if (!cover) {
                    curReadsCount++;
                    continue;
                }

                // 得到fragment对应read的起始、终止位点
                readStart = break_point;
                readEnd = (end_point - break_point >= readLength)? (break_point+readLength):end_point;
                readRecord = readStart + "-" + readEnd;

                Integer countRecord = readsDistribution.getOrDefault(readRecord, 0);
                ArrayList<Integer> fragmentRecord = fragmentEndSites.getOrDefault(readRecord, new ArrayList<>());
                fragmentRecord.add(end_point);
                fragmentEndSites.put(readRecord, fragmentRecord);
                readsDistribution.put(readRecord, countRecord+1);
                curReadsCount++;
            }

            for (String readRange: readsDistribution.keySet()) {
                String[] points = readRange.split("-");
                // 获取reads的起始、终止位点
                break_point = Integer.parseInt(points[0]);
                readEnd = Integer.parseInt(points[1]);
                // 获取reads对应的fragment的终止位点
                ArrayList<Integer> endPoints = fragmentEndSites.get(readRange);
                Collections.shuffle(endPoints);

                // 相同位置的reads数目，计算存在ASE时，major allele的reads数目
                int count = readsDistribution.get(readRange);
                int refCount = (int) (count * refProp);

                // 如果基因含有ASE位点
                if (mutExonSeq != null) {
                    List<Integer> majorAlleleFragmentEnds = endPoints.subList(0, refCount);
                    List<Integer> minorAlleleFragmentEnds = endPoints.subList(refCount, endPoints.size());
                    // 生成major allele的reads
                    for (Integer endPoint: majorAlleleFragmentEnds) {
                        refFragmentString = this.exonSeq.substring(break_point, endPoint);
                        refFragment = this.getSequencingReads(refFragmentString, break_point, endPoint, readLength, direct);
                        if (direct.equals("SE"))
                            this.writeReadInFile(refFragment, seqErrorModel, mateFile1, readLength, break_point, endPoint, "ref");
                        else
                            this.pairReadToFile(refFragment, seqErrorModel, mateFile1, mateFile2, readLength, break_point, endPoint, "ref");
                        refFragment = null;
                        refFragmentString = null;
                    }
                    // 生成minor allele的reads
                    for (Integer endPoint: minorAlleleFragmentEnds) {
                        altFragmentString = mutExonSeq.substring(break_point, endPoint);
                        altFragment = this.getSequencingReads(altFragmentString, break_point, endPoint, readLength, direct);
                        if (direct.equals("SE"))
                            this.writeReadInFile(altFragment, seqErrorModel, mateFile1, readLength, break_point, endPoint, "ref");
                        else
                            this.pairReadToFile(altFragment, seqErrorModel, mateFile1, mateFile2, readLength, break_point, endPoint, "ref");
                        altFragment = null;
                        altFragmentString = null;
                    }
                } else { // 如果基因不包含SNP位点
                    for (Integer endPoint: endPoints) {
                        refFragmentString = this.exonSeq.substring(break_point, endPoint);
                        refFragment = this.getSequencingReads(refFragmentString, break_point, endPoint, readLength, direct);
                        if (direct.equals("SE"))
                            this.writeReadInFile(refFragment, seqErrorModel, mateFile1, readLength, break_point, endPoint, "ref");
                        else
                            this.pairReadToFile(refFragment, seqErrorModel, mateFile1, mateFile2, readLength, break_point, endPoint, "ref");

                        // release memory
                        refFragment = null;
                        refFragmentString = null;
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
//                    for (Peak peak: this.peakList) {
//                        int peakStart = peak.getPeak_start();
//                        int peakEnd = peak.getPeak_end();
//                        int centre = (peakStart + peakEnd)/2;
//                        if (break_point <= centre && centre <= end_point) {
//                            peak.addFragments(fragment);
//                        }
//                    }
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

    /**
     * write generated reads into file
     * @param fragment Fragmentation instance
     * @param seqErrorModel sequencing error model
     * @param fw BufferedWriter
     * @param readLength readLength
     * @param break_point fragment start position on exon sequence
     * @param end_point fragment end position on exon sequence
     * @param type ref or alt
     */
    private void writeReadInFile(Fragmentation fragment, SequencingError seqErrorModel, BufferedWriter fw,
                                 int readLength, int break_point, int end_point, String type) {
        try {
            String sequencingRead = fragment.getSingleEndRead();
            String strand = fragment.getReadStrand();
            sequencingRead = seqErrorModel.pcrErrorReads(sequencingRead);
            if (sequencingRead.length() != readLength) {
                int baseNNum = readLength - sequencingRead.length();
                String baseNSeq = this.fillBaseN(baseNNum);
                sequencingRead = sequencingRead + baseNSeq;
            }

            fw.write(">chr" + this.chr + "_" + this.geneId + "_" +type + "_"+break_point + ":" + end_point +"\t" + strand); // +"_"+readsStart+":"+readsEnd
            fw.newLine();
            fw.write(sequencingRead);
            fw.newLine();
        } catch (IOException | StringIndexOutOfBoundsException e) {
            e.printStackTrace();
        }
    }

    private void pairReadToFile(Fragmentation fragment, SequencingError seqErrorModel, BufferedWriter mateFile1,
                                BufferedWriter mateFile2, int readLength, int break_point, int end_point, String type) {
        try {
            String[] sequencingRead = fragment.getPairEndRead();
            String mate1 = sequencingRead[0];
            String mate2 = sequencingRead[1];
            String strand = fragment.getReadStrand();
            mate1 = seqErrorModel.pcrErrorReads(mate1);
            mate2 = seqErrorModel.pcrErrorReads(mate2);
            if (mate1.length() != readLength) {
                int baseNNum = readLength - mate1.length();
                String baseNSeq = this.fillBaseN(baseNNum);
                mate1 = mate1 + baseNSeq;
                mate2 = baseNSeq + mate2;
            }

            mateFile1.write(">chr" + this.chr + "_" + this.geneId + "_" +type + "_"+break_point + ":" + end_point +"\t" + strand);
            mateFile1.newLine();
            mateFile1.write(mate1);
            mateFile1.newLine();

            mateFile2.write(">chr" + this.chr + "_" + this.geneId + "_" +type + "_"+break_point + ":" + end_point +"\t" + strand);
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

}

