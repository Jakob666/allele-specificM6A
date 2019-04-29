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
    private String geneId, strand;
    private TranscriptRecord longestTranscriptRecord;
    private int geneStart, geneEnd;
    private double RPKM;
    private String chr, exonSeq;
    private int readsCount, refReadsCount = 0, altReadsCount = 0;
    private double[] pmRange;
    private ElementRecord exonList = null;
    private LinkedList<Peak> peakList;
    public Gene(String geneId, int geneStart, int geneEnd, String strand, String chr) {
        this.geneId = geneId;
        this.geneStart = geneStart;
        this.geneEnd = geneEnd;
        this.strand = strand;
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
    public void calculateReadsCount(long librarySize) {
        this.readsCount = (int) ((this.exonSeq.length() * librarySize * this.RPKM) / Math.pow(10, 9));
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
        this.refReadsCount = (int) (this.readsCount * refProp / multiple);
        this.altReadsCount = this.readsCount / multiple - refReadsCount;
        int maxReadCount = Math.max(this.refReadsCount, this.altReadsCount);
        // used to randomly generate fragment length
        NormalDistribution nordi = new NormalDistribution(fragmentMean, fragmentTheta);
        int curReadsCount = 0;
        String fragmentString;
        Fragmentation fragment;
        try {
            int break_point, end_point, fragmentLength;
            // generate alternative reads
            while (curReadsCount < maxReadCount) {
                // randomly generate a fragment length and form exon fragment
                fragmentLength = Math.abs((int) nordi.sample());
                // measure the start and end position of fragment on exon sequence
                int[] breakAndEndPoint = this.getBreakEndPoint(fragmentLength);
                break_point = breakAndEndPoint[0];
                end_point = breakAndEndPoint[1];
                // generate reference reads
                if (curReadsCount < this.refReadsCount) {
                    fragmentString = this.exonSeq.substring(break_point, end_point);
                    fragment = this.getSequencingReads(fragmentString, break_point, end_point, readLength, direct);
                    if (direct.equals("SE"))
                        this.writeReadInFile(fragment, seqErrorModel, mateFile1, readLength, break_point, end_point, "ref");
                    else
                        this.pairReadToFile(fragment, seqErrorModel, mateFile1, mateFile2, readLength, break_point, end_point, "ref");

                }
                // generate alternative reads
                if (curReadsCount < this.altReadsCount) {
                    fragmentString = mutExonSeq.substring(break_point, end_point);
                    fragment = this.getSequencingReads(fragmentString, break_point, end_point, readLength, direct);
                    if (direct.equals("SE"))
                        this.writeReadInFile(fragment, seqErrorModel, mateFile1, readLength, break_point, end_point, "alt");
                    else
                        this.pairReadToFile(fragment, seqErrorModel, mateFile1, mateFile2, readLength, break_point, end_point, "alt");

                }
                fragment = null;

                curReadsCount++;
            }
        } catch (Exception io) {
            io.printStackTrace();
        }
    }

    public void generateIpReads(int fragmentMean, int fragmentTheta, BufferedWriter mateFile1, BufferedWriter mateFile2,
                                int readLength, int multiple, double refProp, String mutExonSeq, String direct,
                                SequencingError seqErrorModel) {
        this.refReadsCount = (int) (this.readsCount * refProp / multiple);
        this.altReadsCount = this.readsCount / multiple - refReadsCount;
        int maxReadCount = Math.max(this.refReadsCount, this.altReadsCount);
        // used to randomly generate fragment length
        NormalDistribution nordi = new NormalDistribution(fragmentMean, fragmentTheta);
        int curReadsCount = 0;
        String fragmentString;
        Fragmentation fragment;
        try {
            int break_point, end_point, fragmentLength;
            // generate alternative reads
            while (curReadsCount < maxReadCount) {
                // randomly generate a fragment length and form exon fragment
                fragmentLength = Math.abs((int) nordi.sample());
                // measure the start and end position of fragment on exon sequence
                int[] breakAndEndPoint = this.getBreakEndPoint(fragmentLength);
                break_point = breakAndEndPoint[0];
                end_point = breakAndEndPoint[1];
                // generate reference reads
                if (curReadsCount < this.refReadsCount) {
                    fragmentString = this.exonSeq.substring(break_point, end_point);
                    fragment = this.getSequencingReads(fragmentString, break_point, end_point, readLength, direct);
                    if (direct.equals("SE"))
                        this.writeReadInFile(fragment, seqErrorModel, mateFile1, readLength, break_point, end_point, "ref");
                    else
                        this.pairReadToFile(fragment, seqErrorModel, mateFile1, mateFile2, readLength, break_point, end_point, "ref");
//                    for (Peak peak: this.peakList) {
//                        int peakStart = peak.getPeak_start();
//                        int peakEnd = peak.getPeak_end();
//                        int centre = (peakStart + peakEnd)/2;
//                        if (break_point <= centre && centre <= end_point) {
//                            peak.addFragments(fragment);
//                        }
//                    }
                }
                // generate alternative reads
                if (curReadsCount < this.altReadsCount) {
                    fragmentString = mutExonSeq.substring(break_point, end_point);
                    fragment = this.getSequencingReads(fragmentString, break_point, end_point, readLength, direct);
                    if (direct.equals("SE"))
                        this.writeReadInFile(fragment, seqErrorModel, mateFile1, readLength, break_point, end_point, "alt");
                    else
                        this.pairReadToFile(fragment, seqErrorModel, mateFile1, mateFile2, readLength, break_point, end_point, "alt");
                }

                curReadsCount++;
            }
        } catch (Exception ie) {
            ie.printStackTrace();
        }
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

