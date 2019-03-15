package heterozygoteSiteAnalysis;

public class IntervalTree {
    public IntervalTreeNode root;

    public IntervalTree() {}

    public IntervalTree(IntervalTreeNode root) {
        this.root = root;
    }

    /**
     * insert a node into an intervaltree
     * @param node the latest insert node of the interval tree
     * @param tree interval tree instance
     */
    public IntervalTree insertNode(IntervalTree tree, IntervalTreeNode node) {
        if (node == null)
            return tree;

        // define a tag to record parent node of the latest insert node
        IntervalTreeNode tag = null;
        IntervalTreeNode treeRoot = tree.root;
        while (treeRoot != null) {
            treeRoot.max = Math.max(treeRoot.max, node.max);
            tag = treeRoot;
            if (node.intervalStart < treeRoot.intervalStart) {
                treeRoot = treeRoot.leftChild;
            } else {
                treeRoot = treeRoot.rightChild;
            }
        }

        node.parent = tag;
        if (tag == null) {
            tree.root = node;
        } else if (node.intervalStart < tag.intervalStart) {
            tag.leftChild = node;
        } else {
            tag.rightChild = node;
        }

        // set the latest insert node color red and fix the tree balance
        node.color = 0;
        tree = this.insertFix(tree, node);

        return tree;
    }

    /**
     * right rotate to maintain the RB tree balance
     * @param tree Interval tree instance
     * @param node the latest insert node
     * @return balanced interval tree instance
     */
    private IntervalTree rightRotate(IntervalTree tree, IntervalTreeNode node) {
        if (tree == null | node == null)
            return tree;
        // after rotating, node become the right child of its current left child
        IntervalTreeNode curLeftChild = node.leftChild;
        node.leftChild = curLeftChild.rightChild;
        if (curLeftChild.rightChild != null) {
            curLeftChild.rightChild.parent = node;
        }

        curLeftChild.parent = node.parent;

        // if node is the tree node
        if (node.parent == null) {
            tree.root = curLeftChild;
        } else if (node == node.parent.leftChild) {
            node.parent.leftChild = curLeftChild;
        } else {
            node.parent.rightChild = curLeftChild;
        }

        curLeftChild.rightChild = node;
        node.parent = curLeftChild;

        // maintain the additional information
        curLeftChild.max = node.max;
        int newLeftChildMax = node.leftChild == null ? Integer.MIN_VALUE : node.leftChild.intervalEnd;
        int newRightChildMax = node.rightChild == null ? Integer.MIN_VALUE : node.rightChild.intervalEnd;
        node.max = Math.max(node.intervalEnd, Math.max(newLeftChildMax, newRightChildMax));

        return tree;
    }

    /**
     * left rotate to maintain the RB tree balance
     * @param tree Interval tree instance
     * @param node the latest insert node
     * @return balanced interval tree instance
     */
    private IntervalTree leftRotate(IntervalTree tree, IntervalTreeNode node) {
        if (tree == null | node == null)
            return tree;
        // after rotating, node become the left child of its current right child
        IntervalTreeNode curRightChild = node.rightChild;
        node.rightChild = curRightChild.leftChild;
        if (curRightChild.leftChild != null) {
            curRightChild.leftChild.parent = node;
        }

        curRightChild.parent = node.parent;
        // if node is the tree root
        if (node.parent == null) {
            tree.root = curRightChild;
        } else if (node == node.parent.leftChild) {
            node.parent.leftChild = curRightChild;
        } else {
            node.parent.rightChild = curRightChild;
        }

        node.parent = curRightChild;
        curRightChild.leftChild = node;
        // maintain additional information
        curRightChild.max = node.max;
        int newLeftChildMax = node.leftChild == null ? Integer.MIN_VALUE : node.leftChild.intervalEnd;
        int newRightChildMax = node.rightChild == null ? Integer.MIN_VALUE : node.rightChild.intervalEnd;
        node.max = Math.max(node.intervalEnd, Math.max(newLeftChildMax, newRightChildMax));

        return tree;
    }

    /**
     * fix the RB-tree to maintain the properties. The implementation of this operation according to the node color of
     * the latest insert node. The main thought of insert fix is to pull the new insert node upward.
     * @param tree Interval tree instance
     * @param node the latest insert node
     * @return balanced interval tree instance
     */
    private IntervalTree insertFix(IntervalTree tree, IntervalTreeNode node) {
        while (node.parent != null && node.parent.color == 0) {
            IntervalTreeNode grandFather = node.parent.parent;
            // if parent node is the left child of grandfather
            if (node.parent == grandFather.leftChild) {
                IntervalTreeNode uncleNode = grandFather.rightChild;
                // AttributeError raise when one of the child is None, the leave node in RB-tree is in black
                int uncleNodeColor = uncleNode != null ? uncleNode.color : 1;

                //both father and uncle nodes are in red
                if (node.parent.color == uncleNodeColor) {
                    grandFather.color = 0;
                    node.parent.color = 1;
                    uncleNode.color = 1;
                    node = grandFather;
                } else {
                    if (node == node.parent.rightChild) {
                        node = node.parent;
                        tree = this.leftRotate(tree, node);
                    }
                    node.parent.color = 1;
                    node.parent.parent.color = 0;
                    tree = this.rightRotate(tree, node.parent.parent);
                }
            } else if (node.parent == grandFather.rightChild) {
                // if parent node is the right child of grandfather
                IntervalTreeNode uncleNode = grandFather.leftChild;
                int uncleNodeColor = uncleNode != null ? uncleNode.color : 1;
                if (node.parent.color == uncleNodeColor) {
                    grandFather.color = 0;
                    node.parent.color = 1;
                    uncleNode.color = 1;
                    node = grandFather;
                } else {
                    if (node == node.parent.leftChild) {
                        node = node.parent;
                        tree = this.rightRotate(tree, node);
                    }
                    node.parent.color = 1;
                    node.parent.parent.color = 0;
                    tree = this.leftRotate(tree, node.parent.parent);
                }
            }
        }
        // root must be in black
        tree.root.color = 1;

        return tree;
    }

    /**
     * search nodes on a RB-tree which their interval contains a certain value
     * @param treeRoot the root node of an IntervalTree instance
     * @param searchValue a certain value
     * @return the interval tree node which contains the value or null
     */
    public IntervalTreeNode search(IntervalTreeNode treeRoot, int searchValue) {
        IntervalTreeNode searchResult = null;
        if (treeRoot == null)
            return searchResult;
        int treeRootStart = treeRoot.intervalStart;
        int treeRootEnd = treeRoot.intervalEnd;
        if (searchValue <= treeRootEnd & searchValue >= treeRootStart)
            searchResult = treeRoot;
        else if (treeRoot.leftChild != null && searchValue < treeRoot.leftChild.max)
            searchResult = this.search(treeRoot.leftChild, searchValue);
        else
            searchResult = this.search(treeRoot.rightChild, searchValue);

        return searchResult;
    }

    /**
     * whether a certain value is on the interval tree
     * @param tree IntervalTree instance
     * @param searchValue a certain value
     * @return true if the value on the tree, otherwise false
     */
    public boolean ifOnTree(IntervalTree tree, int searchValue) {
        boolean onTree = false;
        if (tree == null)
            return onTree;
        IntervalTreeNode searchRes = this.search(tree.root, searchValue);
        if (searchRes != null)
            onTree = true;

        return onTree;
    }

    /**
     * output the  interval tree structure.
     * @param treeRoot IntervalNode object
     */
    public void walkAlongTheTree(IntervalTreeNode treeRoot) {
        if (treeRoot != null) {
            this.walkAlongTheTree(treeRoot.leftChild);
            String nodeColor = treeRoot.color == 1 ? "black" : "red";
            String nodeInterval = "[" + treeRoot.intervalStart + ", " + treeRoot.intervalEnd +"]";
            System.out.println("node: " + nodeInterval + ", color: " + nodeColor);
            this.walkAlongTheTree(treeRoot.rightChild);
        }
    }
}
