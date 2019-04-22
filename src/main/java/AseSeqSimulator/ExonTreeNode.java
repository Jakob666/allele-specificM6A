package AseSeqSimulator;

import heterozygoteSiteAnalysis.IntervalTreeNode;

public class ExonTreeNode extends IntervalTreeNode {
    public String chrNum, geneName;
    public ExonTreeNode parent, leftChild, rightChild;

    public ExonTreeNode(int start, int end, String geneName, String chrNum) {
        super(start, end, start, end);
        this.geneName = geneName;
        this.chrNum = chrNum;
    }
}
