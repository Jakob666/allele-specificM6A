package GTFComponent;

import heterozygoteSiteAnalysis.IntervalTree;

import java.io.*;
import java.util.HashMap;

public class GeneExonIntervalTree {
    private String gtfFile;
    private HashMap<String, HashMap<String, IntervalTree>> geneExonIntervalTree;

    public GeneExonIntervalTree(String gtfFile) {
        this.gtfFile = gtfFile;
    }

    public void generateExonTree() {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.gtfFile))));
            this.geneExonIntervalTree = new HashMap<>();
            String line = "", chrNum, geneId, geneName, strand;
            String[] info, geneInfo;
            int exonStart, exonEnd;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    if (!info[2].equals("exon"))
                        continue;
                    chrNum = info[0];
                    geneInfo = this.getGeneInfo(info[8]);
                    geneId = geneInfo[0];
                    geneName = geneInfo[1];
                    exonStart = Integer.valueOf(info[3]);
                    exonEnd = Integer.valueOf(info[4]);
                    strand = info[6];
                    GTFIntervalTreeNode gitn = new GTFIntervalTreeNode(exonStart, exonEnd, 0, 0, geneName, geneId, strand);

                    HashMap<String, IntervalTree> chrGenes = this.geneExonIntervalTree.getOrDefault(chrNum, new HashMap<>());
                    IntervalTree it = chrGenes.getOrDefault(geneId, new IntervalTree());
                    it = it.insertNode(it, gitn);
                    chrGenes.put(geneId, it);
                    this.geneExonIntervalTree.put(chrNum, chrGenes);
                }
                info = null;
                geneInfo = null;
            }
        } catch (IOException ie) {
            ie.printStackTrace();
            System.exit(2);
        } finally {
            if (bfr != null) {
                try {
                    bfr.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    /**
     * get gene name
     * @param recordInfo GTF information
     * @return gene name
     */
    private String[] getGeneInfo(String recordInfo) {
        String[] info = recordInfo.split("; ");
        String geneName = null, geneId = null;
        for (String s: info) {
            if (s.startsWith("gene_id")) {
                String[] name = s.split(" ");
                geneId = name[1].substring(1, name[1].length() -1);
            }
            if (s.startsWith("gene_name")) {
                String[] name = s.split(" ");
                geneName = name[1].substring(1, name[1].length() -1);
            }
        }

        return new String[] {geneId, geneName};
    }

    public HashMap<String, HashMap<String, IntervalTree>> getGeneExonIntervalTree() {
        return this.geneExonIntervalTree;
    }
}
