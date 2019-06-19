package AseSeqSimulator;

import GTFComponent.GeneRecord;
import GTFComponent.TranscriptRecord;

import java.util.*;

public class CommonMethod {

    /**
     * generate anti-chain for a particular sequence
     * @param Seq sequence in direction 5' -> 3'
     * @return anti-chain sequence
     */
    public static String AntiChain(String Seq) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < Seq.length(); i++) {
            char nucleotide = Seq.charAt(i);
            if (nucleotide == 'A' || nucleotide == 'a') {
                sb.insert(0, 'T');
            } else if (nucleotide == 'T' || nucleotide == 't') {
                sb.insert(0, 'A');
            } else if (nucleotide == 'C' || nucleotide == 'c') {
                sb.insert(0, 'G');
            } else {
                sb.insert(0, 'C');
            }
        }
        return sb.toString();
    }

    /**
     * find the transcript with longest sequence length of a certain gene
     * @param geneRecord GeneRecord instance
     * @return TranscriptRecord instance
     */
    public static TranscriptRecord findLongestTranscript(GeneRecord geneRecord) {
        HashMap<String, TranscriptRecord> transcriptIsoforms = geneRecord.getTranscriptIsoform();
        int max_Length = 0;
        TranscriptRecord transcriptRecord = new TranscriptRecord();
        Set<String> isoformIds = transcriptIsoforms.keySet();
        int transcriptStart, transcriptEnd;
        for (String isoformId: isoformIds) {
            TranscriptRecord record = transcriptIsoforms.get(isoformId);
            transcriptStart = record.getTranscriptStart();
            transcriptEnd = record.getTranscriptEnd();
            int transcript_length = transcriptEnd - transcriptStart;
            if (transcript_length > max_Length) {
                transcriptRecord = record;
                max_Length = transcript_length;
            }
        }
        return transcriptRecord;
    }


}
