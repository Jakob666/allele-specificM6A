package GTFComponent;

import java.util.HashMap;

public class TranscriptRecord {
    private String transcriptId = "unknown", transcriptName = "unknown", strand = "unknown", bioType = "unknown";
    private int start, end;
    private HashMap<String, ElementRecord> elementList;

    public TranscriptRecord() {
        this.elementList = new HashMap<>();
    }

    public void setTranscriptName(String name) {
        this.transcriptName = name;
    }

    public String getTranscriptName() {
        return this.transcriptName;
    }

    public void setTranscriptRange(int start, int end) {
        this.start = start;
        this.end = end;
    }

    public int getTranscriptStart() {
        return this.start;
    }

    public int getTranscriptEnd() {
        return this.end;
    }

    public void setStrand(String strand) {
        this.strand = strand;
    }

    public String getStrand() {
        return this.strand;
    }

    public void setTranscriptId(String transcriptId) {
        this.transcriptId = transcriptId;
    }

    public String getTranscriptId() {
        return this.transcriptId;
    }

    public void setBioType(String bioType) {
        this.bioType = bioType;
    }

    public String getBioType() {
        return this.bioType;
    }

    public void renewElementList(String elementType, ElementRecord elementRecord) {
        ElementRecord firstEr = this.elementList.getOrDefault(elementType, null);
        ElementRecord curEr = firstEr;
        if (firstEr == null)
            firstEr = elementRecord;
        else {
            while (curEr.getNextElement() != null) {
                curEr = curEr.getNextElement();
            }
            curEr.setNextElement(elementRecord);
        }
        this.elementList.put(elementType, firstEr);
    }

    public HashMap<String, ElementRecord> getElementList() {
        return this.elementList;
    }
}
