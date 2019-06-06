package AseSeqSimulator;

/**
 * fragment is a sequence cut from a transcript sequence. Generate sequencing reads from it
 */
public class Fragmentation {
    // mate1Strand值为 same表示与reads与转录本同链，anti则表示与转录本异链
    private String transcriptFragmentSeq, mate1, mate2, mate1Strand = "same";

    /**
     * constructor
     * @param transcriptFragmentSeq fragment sequence extract from exon region
     */
    public Fragmentation(String transcriptFragmentSeq){
        this.transcriptFragmentSeq = transcriptFragmentSeq;
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

    /**
     * generate one sequencing reads based on fragment, the length of which is equals to user input reads
     * @param randNum random number to determine strand
     * @param readLength read length
     */
    public void singleEndSequencingRead(double randNum, int readLength){
        if(this.transcriptFragmentSeq.length() > readLength) {
            this.mate1 = this.transcriptFragmentSeq.substring(0, readLength);
        }else{
            this.mate1 = this.transcriptFragmentSeq;
        }
        if (randNum > 0.5) {
            this.mate1 = CommonMethod.AntiChain(this.mate1);
            this.mate1Strand = "anti";
        }
    }

    public void pairEndSequencingRead(double randNum, int readLength){
        if(this.transcriptFragmentSeq.length() > readLength) {
            this.mate1 = this.transcriptFragmentSeq.substring(0, readLength);
            this.mate2 = this.transcriptFragmentSeq.substring(this.transcriptFragmentSeq.length() - readLength);
        }else{
            this.mate1 = this.transcriptFragmentSeq;
            this.mate2 = this.transcriptFragmentSeq;
        }
        if (randNum > 0.5) {
            this.mate1 = CommonMethod.AntiChain(this.mate1);
            this.mate2 = CommonMethod.AntiChain(this.mate2);
            this.mate1Strand = "anti";
        }
    }
}
