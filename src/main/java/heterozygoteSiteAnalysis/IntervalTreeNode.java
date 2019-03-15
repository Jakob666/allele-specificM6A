package heterozygoteSiteAnalysis;

public class IntervalTreeNode {
    public IntervalTreeNode parent;
    public IntervalTreeNode rightChild, leftChild;
    public int center, intervalStart, intervalEnd, max;
    // color 1 denotes black, 0 represents red
    public int color;

    /**
     * Constructor
     * @param intervalStart range start
     * @param intervalEnd range end
     */
    public IntervalTreeNode(int intervalStart, int intervalEnd) {
        this.intervalStart = intervalStart;
        this.intervalEnd = intervalEnd;
        this.center = (this.intervalStart + this.intervalEnd) / 2;
        // property max is used to record the max range for a Interval tree
        this.max = intervalEnd;
    }
}
