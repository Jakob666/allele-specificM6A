package SamtoolsPileupSNPCalling;

import java.io.*;

public class VcfSearchTree {
    private File vcfFile;

    public VcfSearchTree(File vcfFilePath) {
        this.vcfFile = vcfFilePath;
    }

    public BinarySearchTree buildTree() {
        BinarySearchTree bt = new BinarySearchTree();

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

                    VcfTreeNode newNode = new VcfTreeNode(Integer.parseInt(lineInfo[1]), lineInfo[2], lineInfo[3], lineInfo[4]);
                    bt = bt.insertNode(bt, newNode);
                }
            }

            bfr.close();

        } catch (FileNotFoundException fne) {
            fne.printStackTrace();
            System.exit(2);
        } catch (IOException ie) {
            ie.printStackTrace();
            System.exit(3);
        }

        return bt;
    }

}
