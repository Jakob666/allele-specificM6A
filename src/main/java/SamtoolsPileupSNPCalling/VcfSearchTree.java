package SamtoolsPileupSNPCalling;

import org.apache.log4j.Logger;

import java.io.*;
import java.util.HashMap;

/**
 * deprecated, use VcfSearch instead
 */
public class VcfSearchTree {
    private File vcfFile;
    private Logger log;

    public VcfSearchTree(File vcfFilePath, Logger logger) {
        this.vcfFile = vcfFilePath;
        this.log = logger;
    }

    public HashMap<String, BinarySearchTree> buildTree() {

        HashMap<String, BinarySearchTree> vcfSearchTree = new HashMap<String, BinarySearchTree>();
        try {
            BufferedReader bfr = new BufferedReader(
                    new InputStreamReader(new FileInputStream(this.vcfFile))
            );

            String line = "";
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;

                    String[] lineInfo = line.split("\t");
                    String chrNum = lineInfo[0];
                    // create new node
                    VcfTreeNode newNode = new VcfTreeNode(Integer.parseInt(lineInfo[1]), lineInfo[2], lineInfo[3], lineInfo[4]);
                    // get and renew the tree
                    BinarySearchTree bt = vcfSearchTree.getOrDefault(chrNum, new BinarySearchTree());
                    bt = bt.insertNode(bt, newNode);
                    vcfSearchTree.put(chrNum, bt);
                }
            }
            bfr.close();
        } catch (IOException ie) {
            this.log.error(ie.getMessage());
            System.exit(2);
        }

        return vcfSearchTree;
    }
}
