package heterozygoteSiteAnalysis;

public class IntervalTreeNode {
    public IntervalTreeNode parent;
    public IntervalTreeNode rightChild, leftChild;
    public int center, intervalStart, intervalEnd, max, peakStart, peakEnd;
    // color 1 denotes black, 0 represents red
    public int color;

    /**
     * Constructor
     * @param intervalStart block range start
     * @param intervalEnd block range end
     * @param peakStart peak range start
     * @param peakEnd peak range end
     */
    public IntervalTreeNode(int intervalStart, int intervalEnd, int peakStart, int peakEnd) {
        this.intervalStart = intervalStart;
        this.intervalEnd = intervalEnd;
        this.peakStart = peakStart;
        this.peakEnd = peakEnd;
        this.center = (this.intervalStart + this.intervalEnd) / 2;
        // property max is used to record the max range for a Interval tree
        this.max = intervalEnd;
    }
}
