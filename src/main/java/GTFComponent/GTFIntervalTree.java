package GTFComponent;

import heterozygoteSiteAnalysis.IntervalTree;

import java.io.*;
import java.util.HashMap;

/**
 * 构建GTF文件的区间树
 */
public class GTFIntervalTree {
    private String gtfFile;
    private HashMap<String, IntervalTree> gtfIntervalTrees = new HashMap<>();

    /**
     * Constructor
     * @param gtfFile GTF文件的文件路径
     */
    public GTFIntervalTree(String gtfFile) {
        this.gtfFile = gtfFile;
    }

    /**
     * 解析GTF文件中的基因信息，并给每条染色体建立基因信息的区间树
     */
    public void parseGTFFile() {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.gtfFile))));
            String line = "", chrNum, geneName, geneId, componentType;
            String[] info, geneInfo;
            int start, end;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    chrNum = info[0];
                    componentType = info[2];
                    if (!componentType.equals("gene"))
                        continue;
                    start = Integer.parseInt(info[3]);
                    end = Integer.parseInt(info[4]);
                    geneInfo = this.getGeneInfo(info[8]);
                    geneId = geneInfo[0];
                    geneName = geneInfo[1];

                    GTFIntervalTreeNode gitn = new GTFIntervalTreeNode(start, end, 0, 0, geneName, geneId);
                    IntervalTree it = this.gtfIntervalTrees.getOrDefault(chrNum, new IntervalTree());
                    it = it.insertNode(it, gitn);
                    this.gtfIntervalTrees.put(chrNum, it);
                }
            }
            bfr.close();
        } catch (IOException ie) {
            ie.printStackTrace();
        } finally {
            if (bfr != null) {
                try {
                    bfr.close();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
    }

    /**
     * 从信息文本中获取基因的名称
     * @param recordInfo 传入信息
     * @return 基因名称
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

    /**
     * 获取GTF文件对应的区间树
     * @return 区间树
     */
    public HashMap<String, IntervalTree> getGtfIntervalTrees() {
        return this.gtfIntervalTrees;
    }
}
