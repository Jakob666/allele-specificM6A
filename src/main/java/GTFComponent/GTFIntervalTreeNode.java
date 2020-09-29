package GTFComponent;


import heterozygoteSiteAnalysis.IntervalTreeNode;

public class GTFIntervalTreeNode extends IntervalTreeNode {
    public String geneName, geneId, strand;

    /**
     * Constructor
     * @param intervalStart block range start
     * @param intervalEnd block range end
     */
    public GTFIntervalTreeNode(int intervalStart, int intervalEnd, int peakStart, int peakEnd, String geneName, String geneId, String strand) {
        super(intervalStart, intervalEnd, peakStart, peakEnd);
        this.geneName = geneName;
        this.geneId = geneId;
        this.strand = strand;
    }
}
