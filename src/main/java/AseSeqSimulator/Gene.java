package AseSeqSimulator;

import GTFComponent.ElementRecord;
import GTFComponent.TranscriptRecord;
import PeakSimulator.Peak;
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
    private String geneId, strand;
    private TranscriptRecord longestTranscriptRecord;
    private int geneStart, geneEnd;
    private double RPKM;
    private String chr, exonSeq;
    private int readsCount, refReadsCount = 0, altReadsCount = 0;
    private double[] pmRange;
    private LinkedList<ElementRecord> exonList = new LinkedList<ElementRecord>();
    private LinkedList<Fragmentation> inputFragmentList = new LinkedList<Fragmentation>(), ipFragmentList = new LinkedList<Fragmentation>();
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

    public void setExonSeq(String exonSeq) {
        this.exonSeq = exonSeq;
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

    public LinkedList getInputFragmentList() {
        return this.inputFragmentList;
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

    public LinkedList<ElementRecord> getExonList() {
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
        String exon_splicing = "";
        ElementRecord transcriptExon = elementsOnTranscript.getOrDefault("exon", null);
        try {
            while (transcriptExon != null) {
                String strand = transcriptExon.getStrand();
                int exonStart = transcriptExon.getElementStart();
                int exonEnd = transcriptExon.getElementEnd();
                twoBit.reset();
                String exonSeq = twoBit.loadFragment(exonStart-1, (exonEnd - exonStart + 1));
                this.exonList.add(transcriptExon);
                if (strand.equals("-")) {
                    exon_splicing = CommonMethod.AntiChain(exonSeq) + exon_splicing;
                } else {
                    exon_splicing = exon_splicing + exonSeq;
                }
                transcriptExon = transcriptExon.getNextElement();
            }
        } catch (IOException ie) {
            ie.printStackTrace();
        }
        this.exonSeq = exon_splicing;
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
     * @param fw FileWriter instance
     * @param readLength sequencing read length
     * @param multiple replication of the experiment
     */
    public void generateInputReads(int fragmentMean, int fragmentTheta, BufferedWriter fw, int readLength, int multiple,
                                   double refProp, String mutExonSeq, HashSet<Integer> geneMutationPositions,
                                   SequencingError seqErrorModel) {
        ArrayList<Integer> mutPositions = new ArrayList<>();
        if (geneMutationPositions != null)
            mutPositions = new ArrayList<>(geneMutationPositions);
        this.refReadsCount = (int) (this.readsCount * refProp / multiple);
        this.altReadsCount = this.readsCount / multiple - refReadsCount;
        int maxReadCount = Math.max(this.refReadsCount, this.altReadsCount);
        // used to randomly generate fragment length
        NormalDistribution nordi = new NormalDistribution(fragmentMean, fragmentTheta);
        int curReadsCount = 0;
        Fragmentation fragment;
        UniformIntegerDistribution uniformdi;
        String strand;
        try {
            String fragmentString;
            int break_point, end_point, fragmentLength;
            // generate alternative reads
            while (curReadsCount < maxReadCount) {
                // randomly generate a fragment length and form exon fragment
                fragmentLength = Math.abs((int) nordi.sample());
                // measure the start and end position of fragment on exon sequence
                if(this.exonSeq.length() <= fragmentLength){
                    break_point = 0;
                    end_point = mutExonSeq.length();
                }else if (geneMutationPositions == null) { // if is not mutated gene
                    uniformdi = new UniformIntegerDistribution(0, this.exonSeq.length() - fragmentLength);
                    break_point = Math.abs(uniformdi.sample());
                    end_point = break_point + fragmentLength;
                } else{ // if is mutated gene
                    Collections.shuffle(mutPositions);
                    int pos = mutPositions.get(0);
                    uniformdi = new UniformIntegerDistribution(0, readLength/2);
                    int distance = uniformdi.sample();
                    break_point = (pos - distance > 0)? pos-distance : pos;
                    end_point = (break_point + fragmentLength < mutExonSeq.length())? break_point + fragmentLength: mutExonSeq.length();
                }
                // generate reference reads
                if (curReadsCount < this.refReadsCount) {
                    fragmentString = this.exonSeq.substring(break_point, end_point);
                    fragment = new Fragmentation(fragmentString, break_point, end_point);
                    if (Math.random() < 0.5) {
                        strand = "+";
                        fragment.getSequencingRead("+", readLength);
                    } else {
                        strand = "-";
                        fragment.getSequencingRead("-", readLength);
                    }
                    this.writeReadInFile(fragment, seqErrorModel, fw, readLength, break_point, end_point, strand,"ref");
                    this.inputFragmentList.add(fragment);
                }
                // generate alternative reads
                if (curReadsCount < this.altReadsCount) {
                    fragmentString = mutExonSeq.substring(break_point, end_point);
                    fragment = new Fragmentation(fragmentString, break_point, end_point);
                    if (Math.random() < 0.5) {
                        strand = "+";
                        fragment.getSequencingRead("+", readLength);
                    } else {
                        strand = "-";
                        fragment.getSequencingRead("-", readLength);
                    }
                    this.writeReadInFile(fragment, seqErrorModel, fw, readLength, break_point, end_point, strand, "alt");
                    this.inputFragmentList.add(fragment);
                }

                curReadsCount++;
            }
        } catch (Exception io) {
            io.printStackTrace();
        }
    }

    public void generateIpReads(int fragmentMean, int fragmentTheta, BufferedWriter fw, int readLength, int multiple,
                                double refProp, String mutExonSeq, HashSet<Integer> geneMutationPositions,
                                SequencingError seqErrorModel) {
        ArrayList<Integer> mutPositions = new ArrayList<>();
        if (geneMutationPositions != null)
            mutPositions = new ArrayList<>(geneMutationPositions);
        this.refReadsCount = (int) (this.readsCount * refProp / multiple);
        this.altReadsCount = this.readsCount / multiple - refReadsCount;
        int maxReadCount = Math.max(this.refReadsCount, this.altReadsCount);
        // used to randomly generate fragment length
        NormalDistribution nordi = new NormalDistribution(fragmentMean, fragmentTheta);
        int curReadsCount = 0;
        Fragmentation fragment;
        UniformIntegerDistribution uniformdi;
        String strand;
        try {
            String fragmentString;
            int break_point, end_point, fragmentLength;
            // generate alternative reads
            while (curReadsCount < maxReadCount) {
                // randomly generate a fragment length and form exon fragment
                fragmentLength = Math.abs((int) nordi.sample());
                // measure the start and end position of fragment on exon sequence
                if(this.exonSeq.length() <= fragmentLength){
                    break_point = 0;
                    end_point = mutExonSeq.length();
                }else if (geneMutationPositions == null) { // if is not mutated gene
                    uniformdi = new UniformIntegerDistribution(0, this.exonSeq.length() - fragmentLength);
                    break_point = Math.abs(uniformdi.sample());
                    end_point = break_point + fragmentLength;
                } else{ // if is mutated gene
                    Collections.shuffle(mutPositions);
                    int pos = mutPositions.get(0);
                    uniformdi = new UniformIntegerDistribution(0, readLength/2);
                    int distance = uniformdi.sample();
                    break_point = (pos - distance > 0)? pos-distance : pos;
                    end_point = (break_point + fragmentLength < mutExonSeq.length())? break_point + fragmentLength: mutExonSeq.length();
                }
                // generate reference reads
                if (curReadsCount < this.refReadsCount) {
                    fragmentString = this.exonSeq.substring(break_point, end_point);
                    fragment = new Fragmentation(fragmentString, break_point, end_point);
                    if (Math.random() < 0.5) {
                        strand = "+";
                        fragment.getSequencingRead("+", readLength);
                    } else {
                        strand = "-";
                        fragment.getSequencingRead("-", readLength);
                    }
                    this.writeReadInFile(fragment, seqErrorModel, fw, readLength, break_point, end_point, strand,"ref");
                    this.ipFragmentList.add(fragment);

                    for (Peak peak: this.peakList) {
                        int peakStart = peak.getPeak_start();
                        int peakEnd = peak.getPeak_end();
                        int centre = (peakStart + peakEnd)/2;
                        if (break_point <= centre && centre <= end_point) {
                            peak.addFragments(fragment);
                        }
                    }
                    this.ipFragmentList.add(fragment);
                }
                // generate alternative reads
                if (curReadsCount < this.altReadsCount) {
                    fragmentString = mutExonSeq.substring(break_point, end_point);
                    fragment = new Fragmentation(fragmentString, break_point, end_point);
                    if (Math.random() < 0.5) {
                        strand = "+";
                        fragment.getSequencingRead("+", readLength);
                    } else {
                        strand = "-";
                        fragment.getSequencingRead("-", readLength);
                    }
                    this.writeReadInFile(fragment, seqErrorModel, fw, readLength, break_point, end_point, strand, "alt");
                    this.ipFragmentList.add(fragment);

                    for (Peak peak: this.peakList) {
                        int peakStart = peak.getPeak_start();
                        int peakEnd = peak.getPeak_end();
                        int centre = (peakStart + peakEnd)/2;
                        if (break_point <= centre && centre <= end_point) {
                            peak.addFragments(fragment);
                        }
                    }
                    this.ipFragmentList.add(fragment);
                }

                curReadsCount++;
            }
        } catch (Exception ie) {
            ie.printStackTrace();
        }
    }

    /**
     * write generated reads into file
     * @param fragment Fragmentation instance
     * @param seqErrorModel sequencing error model
     * @param fw BufferedWriter
     * @param readLength readLength
     * @param break_point fragment start position on exon sequence
     * @param end_point fragment end position on exon sequence
     * @param strand read strand
     * @param type ref or alt
     */
    private void writeReadInFile(Fragmentation fragment, SequencingError seqErrorModel, BufferedWriter fw,
                                 int readLength, int break_point, int end_point, String strand, String type) {
        try {
            String sequencingRead = fragment.getReadSeq();
            sequencingRead = seqErrorModel.pcrErrorReads(sequencingRead);
            if (sequencingRead.length() != readLength) {
                int baseNNum = readLength - sequencingRead.length();
                String baseNSeq = this.fillBaseN(baseNNum);
                sequencingRead = sequencingRead + baseNSeq;
            }

            int readsStart = fragment.getFragmentStart();
            int readsEnd = (readLength < fragment.getFragmentLength()) ? readsStart + readLength - 1 : fragment.getFragmentEnd();
            fw.write(">chr" + this.chr + "_" + this.geneId + "_" +type + "_"+break_point + ":" + end_point +"\t" + strand); // +"_"+readsStart+":"+readsEnd
            fw.newLine();
            fw.write(sequencingRead);
            fw.newLine();
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

}

