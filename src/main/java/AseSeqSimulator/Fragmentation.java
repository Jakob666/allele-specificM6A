package AseSeqSimulator;

/**
 * fragment is a sequence cut from a transcript sequence. Generate sequencing reads from it
 */
public class Fragmentation {
    private int fragmentStart, fragmentEnd;
    private String transcriptFragmentSeq, readSeq;
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

    public String getReadSeq(){
        return this.readSeq;
    }

    public int getFragmentLength() {
        return this.fragmentLength;
    }

    /**
     * generate one sequencing reads based on fragment, the length of which is equals to user input reads
     * @param strand positive or negative strand
     * @param readLength read length
     */
    public void getSequencingRead(String strand, int readLength){
        if(this.transcriptFragmentSeq.length() > readLength) {
            this.readSeq = this.transcriptFragmentSeq.substring(0, readLength);
        }else{
            this.readSeq = this.transcriptFragmentSeq;
        }
        if (strand.equals("-"))
            this.readSeq = CommonMethod.AntiChain(this.readSeq);
    }
}
