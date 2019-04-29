package AseSeqSimulator;

/**
 * fragment is a sequence cut from a transcript sequence. Generate sequencing reads from it
 */
public class Fragmentation {
    private int fragmentStart, fragmentEnd;
    private String transcriptFragmentSeq, mate1, mate2, mate1Strand = "+";
    private int fragmentLength;

    /**
     * constructor
     * @param transcriptFragmentSeq fragment sequence extract from exon region
     * @param fragmentStart fragment start point on exon sequence
     * @param fragmentEnd fragment end point on exon sequence
     */
    public Fragmentation(String transcriptFragmentSeq, int fragmentStart, int fragmentEnd){
        this.transcriptFragmentSeq = transcriptFragmentSeq;
        this.fragmentStart = fragmentStart;
        this.fragmentEnd = fragmentEnd;
        this.fragmentLength = transcriptFragmentSeq.length();
    }

    public int getFragmentStart(){
        return this.fragmentStart;
    }

    public int getFragmentEnd(){
        return this.fragmentEnd;
    }

    public String getTranscriptSeq(){
        return this.transcriptFragmentSeq;
    }

    public String getSingleEndRead() {
        return this.mate1;
    }

    public String getReadStrand() {
        return this.mate1Strand;
    }

    public String[] getPairEndRead() {
        return new String[]{this.mate1, this.mate2};
    }

    public int getFragmentLength() {
        return this.fragmentLength;
    }

    /**
     * generate one sequencing reads based on fragment, the length of which is equals to user input reads
     * @param strand positive or negative strand
     * @param readLength read length
     */
    public void singleEndSequencingRead(String strand, int readLength){
        if(this.transcriptFragmentSeq.length() > readLength) {
            this.mate1 = this.transcriptFragmentSeq.substring(0, readLength);
        }else{
            this.mate1 = this.transcriptFragmentSeq;
        }
        if (strand.equals("-")) {
            this.mate1 = CommonMethod.AntiChain(this.mate1);
            this.mate1Strand = "-";
        }
    }

    public void pairEndSequencingRead(String strand, int readLength){
        if(this.transcriptFragmentSeq.length() > readLength) {
            this.mate1 = this.transcriptFragmentSeq.substring(0, readLength);
            this.mate2 = this.transcriptFragmentSeq.substring(this.transcriptFragmentSeq.length() - readLength);
        }else{
            this.mate1 = this.transcriptFragmentSeq;
            this.mate2 = this.transcriptFragmentSeq;
        }
        if (strand.equals("-")) {
            this.mate1 = CommonMethod.AntiChain(this.mate1);
            this.mate2 = CommonMethod.AntiChain(this.mate2);
            this.mate1Strand = "-";
        }
    }
}
