import heterozygoteSiteAnalysis.IntervalTree;
import heterozygoteSiteAnalysis.IntervalTreeNode;
import org.junit.Test;

public class IntervalTreeTest {

    @Test
    public void intervalTreeTest() {
        IntervalTreeNode[] nodes = new IntervalTreeNode[]{new IntervalTreeNode(16, 21),
                new IntervalTreeNode(8, 9), new IntervalTreeNode(5, 8),
                new IntervalTreeNode(25, 30), new IntervalTreeNode(15, 23),
                new IntervalTreeNode(17, 19), new IntervalTreeNode(26, 26),
                new IntervalTreeNode(0, 3), new IntervalTreeNode(6, 10),
                new IntervalTreeNode(19, 20)};

        IntervalTree it = new IntervalTree();
        for (IntervalTreeNode itn : nodes) {
            it = it.insertNode(it, itn);
        }
        System.out.println("root: " + it.root.intervalStart + "->" + it.root.intervalEnd);
        System.out.println(it.root.leftChild.intervalStart + "->" + it.root.leftChild.intervalEnd);

        it.walkAlongTheTree(it.root);

        IntervalTreeNode res = it.search(it.root, 7);

        boolean onTree = it.ifOnTree(it, 7);
        System.out.println(onTree);
        System.out.println("7 in interval [" + res.intervalStart + ", " + res.intervalEnd + "]");
    }
}
