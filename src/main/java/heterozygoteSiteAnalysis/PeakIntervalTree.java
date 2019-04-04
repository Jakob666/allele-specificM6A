package heterozygoteSiteAnalysis;

import org.apache.log4j.Logger;

import java.io.*;
import java.util.HashMap;


public class PeakIntervalTree {
    private File peakBedFile;
    private Logger logger;

    public PeakIntervalTree(String resultBedFile, Logger logger) {
        this.logger = logger;
        this.peakBedFile = new File(resultBedFile);
    }

    public HashMap<String, HashMap<String, IntervalTree>> getPeakTrees() {
        HashMap<String, HashMap<String, IntervalTree>> treeMap = this.buildPeakIntervalTree();
        return treeMap;
    }

    /**
     * get field name index of the peak calling result file
     */
    private HashMap<String, Integer> getFieldIndex() {
        HashMap<String, Integer> fieldIdx = new HashMap<String, Integer>();
        BufferedReader bfr;
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
        } catch (IOException ie) {
            this.logger.error("read bed file failed");
            this.logger.error(ie.getMessage());
            System.exit(2);
        }

        return fieldIdx;
    }

    /**
     * build peak interval tree for each strand of each chromosome
     * @return a nested HashMap, the keys of outer layer are chromosome number, and inner layer keys are pos and neg strand.
     */
    private HashMap<String, HashMap<String, IntervalTree>> buildPeakIntervalTree() {
        HashMap<String, Integer> fieldIdx = this.getFieldIndex();
        String preChr = null;
        HashMap<String, HashMap<String, IntervalTree>> peakTreeMap = new HashMap<>();
        BufferedReader bfr;
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

                    if (Integer.parseInt(blockCount) > 1) {
                        IntervalTreeNode[] nodes = this.multiBlock(start, end, blockSize, blockStart);
                        for (IntervalTreeNode node : nodes) {
                            IntervalTree it = peakTreeMap.get(chrNum).get(strand);
                            it = it.insertNode(it, node);
                            peakTreeMap.get(chrNum).put(strand, it);
                        }
                    } else {
                        int peakStart = Integer.parseInt(start);
                        int peakEnd = Integer.parseInt(end);
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
        }

        return peakTreeMap;
    }

    /**
     * interval tree nodes when a peak contains multiple blocks
     * @param start chromosome start site
     * @param blockSize block size records
     * @param blockStarts block start records
     * @return interval tree node list
     */
    private IntervalTreeNode[] multiBlock(String start, String end, String blockSize, String blockStarts) {
        String[] sizes = blockSize.split(",");
        String[] starts = blockStarts.split(",");
        int peakStart = Integer.parseInt(start);
        int peakEnd = Integer.parseInt(end);

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
