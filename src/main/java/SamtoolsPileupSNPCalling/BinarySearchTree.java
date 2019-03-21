package SamtoolsPileupSNPCalling;



public class BinarySearchTree {
    public VcfTreeNode root;

    public BinarySearchTree() {}

    public BinarySearchTree(VcfTreeNode root) {
        this.root = root;
    }

    public BinarySearchTree insertNode(BinarySearchTree tree, VcfTreeNode newNode) {
        if (newNode == null)
            return tree;

        VcfTreeNode treeRoot = tree.root;
        VcfTreeNode tag = null;

        while (treeRoot != null) {
            tag = treeRoot;
            if (treeRoot.position > newNode.position) {
                treeRoot = treeRoot.leftChild;
            } else {
                treeRoot = treeRoot.rightChild;
            }
        }

        newNode.parent = tag;

        if (tag == null) {
            tree.root = newNode;
        } else if (tag.position > newNode.position) {
            tag.leftChild = newNode;
        } else {
            tag.rightChild = newNode;
        }

        // set the latest insert node color red and fix the tree balance
        newNode.color = 0;
        tree = this.fixTree(tree, newNode);

        return tree;
    }

    private BinarySearchTree fixTree(BinarySearchTree tree, VcfTreeNode node) {
        while (node.parent != null && node.parent.color == 0) {
            VcfTreeNode grandFather = node.parent.parent;
            // if parent node is the left child of grandfather
            if (node.parent == grandFather.leftChild) {
                VcfTreeNode uncleNode = grandFather.rightChild;
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
                VcfTreeNode uncleNode = grandFather.leftChild;
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

    private BinarySearchTree leftRotate(BinarySearchTree tree, VcfTreeNode node) {
        if (tree == null | node == null)
            return tree;
        // after rotating, node become the left child of its current right child
        VcfTreeNode curRightChild = node.rightChild;
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

        return tree;
    }

    private BinarySearchTree rightRotate(BinarySearchTree tree, VcfTreeNode node) {
        if (tree == null | node == null)
            return tree;

        // after rotating, node become the right child of its current left child
        VcfTreeNode curLeftChild = node.leftChild;
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

        return tree;
    }

    public VcfTreeNode search(VcfTreeNode root, int position) {
        if (root == null)
            return null;
        if (root.position == position) {
            return root;
        } else if (root.position > position) {
            return this.search(root.leftChild, position);
        } else {
            return this.search(root.rightChild, position);
        }
    }

    /**
     * output the  interval tree structure.
     * @param treeRoot IntervalNode object
     */
    public void walkAlongTheTree(VcfTreeNode treeRoot) {
        if (treeRoot != null) {
            this.walkAlongTheTree(treeRoot.leftChild);
            String nodeColor = treeRoot.color == 1 ? "black" : "red";
            System.out.println("node: " + treeRoot.position + ", color: " + nodeColor);
            this.walkAlongTheTree(treeRoot.rightChild);
        }
    }
}
