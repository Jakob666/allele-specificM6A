package AseSeqSimulator;

public class VcfRecord {
    public int vcfSite;
    public String chrNum, vcfId, ref, alt;

    public VcfRecord(String chrNum, int vcfSite, String vcfId, String ref, String alt) {
        this.chrNum = chrNum;
        this.vcfSite = vcfSite;
        this.vcfId = vcfId;
        this.ref= ref;
        this.alt =alt;
    }
}
