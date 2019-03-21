package SamtoolsPileupSNPCalling;

public class VcfTreeNode {
    public int position;
    public String id, refNucleotide, altNucleotide;
    public VcfTreeNode leftChild, rightChild, parent;
    public int color;

    public VcfTreeNode(int position, String id, String refNc, String altNc) {
        this.position = position;
        this.id = id;
        this.refNucleotide = refNc;
        this.altNucleotide = altNc;
    }
}
