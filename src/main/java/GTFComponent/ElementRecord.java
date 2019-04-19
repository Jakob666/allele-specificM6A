package GTFComponent;

public class ElementRecord {
    private String elementType, strand;
    private int start, end;
    private ElementRecord nextElement = null;
    private int frame;

    public ElementRecord(String elementType, String strand, int start, int end) {
        this.elementType = elementType;
        this.strand = strand;
        this.start = start;
        this.end = end;
    }

    public int getElementStart() {
        return this.start;
    }

    public int getElementEnd() {
        return this.end;
    }

    public String getElementName() {
        return this.elementType;
    }

    public String getStrand() {
        return this.strand;
    }

    public void setNextElement(ElementRecord elementRecord) {
        this.nextElement = elementRecord;
    }

    public ElementRecord getNextElement() {
        return this.nextElement;
    }

    public void setFrame(int frame) {
        this.frame = frame;
    }

    public int getFrame() {
        return this.frame;
    }
}
