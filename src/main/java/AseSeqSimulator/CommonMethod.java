package AseSeqSimulator;

import GTFComponent.ElementRecord;
import GTFComponent.GeneRecord;
import GTFComponent.TranscriptRecord;

import java.io.IOException;
import java.util.*;

public class CommonMethod {

    /**
     * generate anti-chain for a particular sequence
     * @param Seq sequence in direction 5' -> 3'
     * @return anti-chain sequence
     */
    public static String AntiChain(String Seq) {
        StringBuffer sb = new StringBuffer();
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

    /**
     * 获取基因的外显子组
     * @param geneRecord GeneRecord对象
     * @param twoBit TwoBit Parser对象
     * @return 外显子区间(合并重叠区域)
     */
    public static ArrayList<int[]> getGeneExome(GeneRecord geneRecord, TwoBitParser twoBit) {
        HashMap<String, TranscriptRecord> transcriptIsoforms = geneRecord.getTranscriptIsoform();
        ArrayList<int[]> exonRanges = new ArrayList<>();
        for (String isoformId: transcriptIsoforms.keySet()) {
            TranscriptRecord record = transcriptIsoforms.get(isoformId);
            ElementRecord exon = record.getElementList().getOrDefault("exon", null);
            while (exon != null) {
                int elementStart = exon.getElementStart();
                int elementEnd = exon.getElementEnd();
                exonRanges.add(new int[]{elementStart, elementEnd});
                exon = exon.getNextElement();
            }
        }
        // 将外显子按顺序从前到后排列
        Collections.sort(exonRanges, new Comparator<int[]>() {
            @Override
            public int compare(int[] o1, int[] o2) {
                return Integer.compare(o1[0], o2[0]);
            }
        });
        // 不同的transcript的外显子之可能存在重叠将重叠区域合并
        ArrayList<int[]> mergeOverlapExons = new ArrayList<>();
        for (int[] curExon: exonRanges) {
            if (mergeOverlapExons.size() == 0)
                mergeOverlapExons.add(curExon);
            else {
                int[] preExon = mergeOverlapExons.get(mergeOverlapExons.size()-1);
                if (preExon[1] > curExon[0]) {
                    preExon[1] = Math.max(preExon[1], curExon[1]);
                }  else {
                    mergeOverlapExons.add(curExon);
                }
            }
        }
        return mergeOverlapExons;
    }
}
