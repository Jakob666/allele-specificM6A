package AseSeqSimulator;

import GTFComponent.ElementRecord;
import GTFComponent.TranscriptRecord;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;

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
    private LinkedList<Fragmentation> inputFragmentList = new LinkedList<Fragmentation>();

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

    public LinkedList<ElementRecord> getExonList() {
        return this.exonList;
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
                int feature_start = transcriptExon.getElementStart();
                int feature_end = transcriptExon.getElementEnd();
                twoBit.reset();
                String featureSeq = twoBit.loadFragment(feature_start-1, (feature_end - feature_start + 1));
                this.exonList.add(transcriptExon);
                if (strand.equals("-")) {
                    exon_splicing = exon_splicing + CommonMethod.AntiChain(featureSeq);
                } else {
                    exon_splicing = exon_splicing + featureSeq;
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
    public void generateReads(int fragmentMean, int fragmentTheta, BufferedWriter fw, int readLength, int multiple,
                              HashSet<Integer> geneMutationPositions, SequencingError seqErrorModel) {

        NormalDistribution nordi = new NormalDistribution(fragmentMean, fragmentTheta);
        int curReadsCount = 0;
        Fragmentation fragment;
        UniformRealDistribution uniformdi;
        try {
            String fragmentString;
            int break_point, end_point, fragmentLength;
            while (curReadsCount < this.readsCount/multiple) {
                // first generate fragment from exon sequence
                fragmentLength = Math.abs((int) nordi.sample());
                if(this.exonSeq.length() <= fragmentLength){
                    break_point = 0;
                    end_point = this.exonSeq.length();
                }else{
                    uniformdi = new UniformRealDistribution(0, this.exonSeq.length() - fragmentLength);
                    break_point= Math.abs((int)uniformdi.sample());
                    end_point = break_point + fragmentLength;
                }
                fragmentString = this.exonSeq.substring(break_point, end_point);
                fragment = new Fragmentation(fragmentString, break_point, end_point);

                if (geneMutationPositions != null && !this.isAltRead(geneMutationPositions, break_point, end_point))
                    this.altReadsCount ++;
                else
                    this.refReadsCount ++;

                // generate a read from fragment
                if (Math.random() < 0.5) {
                    fragment.getSequencingRead("+", readLength);
                } else {
                    fragment.getSequencingRead("-", readLength);
                }
                this.inputFragmentList.add(fragment);

                // sequencing reads may contains sequencing error
                String sequencingRead = fragment.getReadSeq();
                sequencingRead = seqErrorModel.pcrErrorReads(sequencingRead);

                int readsStart = fragment.getFragmentStart();
                int readsEnd = (readLength < fragment.getFragmentLength()) ? readsStart + readLength - 1: fragment.getFragmentEnd();
                // write read in fasta format ">chrNum_geneId_fragmentStart_fragmentEnd"
                fw.write(">chr" + this.chr + "_" + this.geneId + "_" + break_point + ":" + end_point+"_"+readsStart+":"+readsEnd);
                fw.newLine();
                fw.write(sequencingRead);
                fw.newLine();

                curReadsCount++;
            }
        } catch (Exception io) {
            io.printStackTrace();
        }
    }

    private boolean isAltRead(HashSet<Integer> mutationPositions, int breakPoint, int endPoint) {
        boolean alt = false;
        for (Integer position: mutationPositions) {
            if (position >= breakPoint && position <= endPoint) {
                alt = true;
                break;
            }
        }
        return alt;
    }
}

