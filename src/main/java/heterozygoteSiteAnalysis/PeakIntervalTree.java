package heterozygoteSiteAnalysis;

import org.apache.log4j.Logger;

import java.io.*;
import java.math.BigDecimal;
import java.util.HashMap;


public class PeakIntervalTree {
    private File peakBedFile;
    private Logger logger;

    public PeakIntervalTree(String resultBedFile, Logger logger) {
        this.logger = logger;
        this.peakBedFile = new File(resultBedFile);
    }

    /**
     * 获取bed文件对应的 peak构成的区间树
     * @return 每条染色体正负链上的peak分别建立区间树构成的HashMap
     */
    public HashMap<String, HashMap<String, IntervalTree>> getPeakTrees() {
        HashMap<String, HashMap<String, IntervalTree>> treeMap = this.buildPeakIntervalTree();
        return treeMap;
    }

    /**
     * 解析bed文件头的内容
     */
    private HashMap<String, Integer> getFieldIndex() {
        HashMap<String, Integer> fieldIdx = new HashMap<String, Integer>();
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(
                    new InputStreamReader(new FileInputStream(this.peakBedFile))
            );
            String line = bfr.readLine();
            String[] columns = line.split("\t");
            int idx = 0;
            // the first row is field names
            for (String col: columns) {
                switch (col) {
                    case "# chr": fieldIdx.put("chr", idx);
                    case "chromStart": fieldIdx.put("chromstart", idx);
                    case "chromEnd": fieldIdx.put("chromend", idx);
                    case "strand": fieldIdx.put("strand", idx);
                    case "blockCount": fieldIdx.put("blockcount", idx);
                    case "blockSizes": fieldIdx.put("blocksizes", idx);
                    case "blockStarts": fieldIdx.put("blockstarts", idx);
                }
                idx += 1;
            }
            bfr.close();
        } catch (IOException ie) {
            this.logger.error("read bed file failed");
            this.logger.error(ie.getMessage());
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

        return fieldIdx;
    }

    /**
     * 对每条染色体正负链上的peak分别建立区间树
     * @return a nested HashMap, the keys of outer layer are chromosome number, and inner layer keys are pos and neg strand.
     */
    private HashMap<String, HashMap<String, IntervalTree>> buildPeakIntervalTree() {
        // 获取个标题对应的索引
        HashMap<String, Integer> fieldIdx = this.getFieldIndex();
        // 为每个染色体构建peak的区间树
        String preChr = null;
        HashMap<String, HashMap<String, IntervalTree>> peakTreeMap = new HashMap<>();
        BufferedReader bfr = null;
        String chrNum, start, end, strand, blockCount, blockSize, blockStart;
        try {
            bfr = new BufferedReader(
                    new InputStreamReader(new FileInputStream(this.peakBedFile))
            );
            String line = "";
            String[] lineInfo;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    lineInfo = line.split("\t");
                    chrNum = lineInfo[fieldIdx.get("chr")];
                    // 遇到新的染色体
                    if (!peakTreeMap.containsKey(chrNum)) {
                        IntervalTree posStrandTree = new IntervalTree();
                        IntervalTree negStrandTree = new IntervalTree();
                        HashMap<String, IntervalTree> peakTree = new HashMap<>();
                        peakTree.put("+", posStrandTree);
                        peakTree.put("-", negStrandTree);

                        peakTreeMap.put(chrNum, peakTree);
                    }
                    strand = lineInfo[fieldIdx.get("strand")];
                    start = lineInfo[fieldIdx.get("chromstart")];
                    end = lineInfo[fieldIdx.get("chromend")];
                    blockCount = lineInfo[fieldIdx.get("blockcount")];
                    blockSize = lineInfo[fieldIdx.get("blocksizes")];
                    blockStart = lineInfo[fieldIdx.get("blockstarts")];
                    // 如果block size大于1，则需要对该记录进行解析；反之，则直接加入树即可
                    if (Integer.parseInt(blockCount) > 1) {
                        IntervalTreeNode[] nodes = this.multiBlock(start, end, blockSize, blockStart);
                        for (IntervalTreeNode node : nodes) {
                            IntervalTree it = peakTreeMap.get(chrNum).get(strand);
                            it = it.insertNode(it, node);
                            peakTreeMap.get(chrNum).put(strand, it);
                        }
                    } else {
                        int peakStart, peakEnd;
                        try {
                            peakStart = Integer.parseInt(start);
                            peakEnd = Integer.parseInt(end);
                        } catch (NumberFormatException nfe) {
                            peakStart = Integer.parseInt(new BigDecimal(start).toPlainString());
                            peakEnd = Integer.parseInt(new BigDecimal(end).toPlainString());
                        }

                        IntervalTreeNode newNode = new IntervalTreeNode(peakStart, peakEnd, peakStart, peakEnd);
                        IntervalTree it = peakTreeMap.get(chrNum).get(strand);
                        it = it.insertNode(it, newNode);
                        peakTreeMap.get(chrNum).put(strand, it);
                    }
                    preChr = chrNum;
                }
            }

            bfr.close();

        } catch (IOException ie) {
            this.logger.error("build m6a peak tree failed");
            this.logger.error(ie.getMessage());
        } finally {
            if (bfr != null) {
                try {
                    bfr.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }

        return peakTreeMap;
    }

    /**
     * 当一个peak的block size大于1时，对该记录进行解析
     * @param start chromosome start site
     * @param blockSize block size records
     * @param blockStarts block start records
     * @return interval tree node list
     */
    private IntervalTreeNode[] multiBlock(String start, String end, String blockSize, String blockStarts) {
        String[] sizes = blockSize.split(",");
        String[] starts = blockStarts.split(",");
        int peakStart, peakEnd;
        try {
            peakStart = Integer.parseInt(start);
            peakEnd = Integer.parseInt(end);
        } catch (NumberFormatException nfe) {
            peakStart = Integer.parseInt(new BigDecimal(start).toPlainString());
            peakEnd = Integer.parseInt(new BigDecimal(end).toPlainString());
        }

        IntervalTreeNode[] nodes = new IntervalTreeNode[starts.length];
        for (int i = 0; i < starts.length; i++) {
            int intervalStart = peakStart + Integer.parseInt(starts[i]);
            int intervalEnd = intervalStart + Integer.parseInt(sizes[i]);
            IntervalTreeNode newNode = new IntervalTreeNode(intervalStart, intervalEnd, peakStart, peakEnd);
            nodes[i] = newNode;
        }

        return nodes;
    }
}
