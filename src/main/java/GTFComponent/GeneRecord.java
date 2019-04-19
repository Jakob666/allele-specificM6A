package GTFComponent;

import java.util.HashMap;

public class GeneRecord {
    private String geneName = "unknown",strand = "unknown", geneId = "unknown", bioType = "unknown";
    private int start, end;
    private HashMap<String, TranscriptRecord> transcriptIsoform;

    public GeneRecord() {
        this.transcriptIsoform = new HashMap<>();
    }

    public void setGeneName(String name) {
        this.geneName = name;
    }

    public String getGeneName() {
        return this.geneName;
    }

    public void setGeneRange(int start, int end) {
        this.start = start;
        this.end = end;
    }

    public int getGeneStart() {
        return this.start;
    }

    public int getGeneEnd() {
        return this.end;
    }

    public void setStrand(String strand) {
        this.strand = strand;
    }

    public String getStrand() {
        return this.strand;
    }

    public void setGeneId(String geneId) {
        this.geneId = geneId;
    }

    public String getGeneId() {
        return this.geneId;
    }

    public void setBioType(String bioType) {
        this.bioType = bioType;
    }

    public String getBioType() {
        return this.bioType;
    }

    public void renewTranscriptIsoform(String isoformName, TranscriptRecord transcript) {
        this.transcriptIsoform.put(isoformName, transcript);
    }

    public HashMap<String, TranscriptRecord> getTranscriptIsoform() {
        return this.transcriptIsoform;
    }
}
