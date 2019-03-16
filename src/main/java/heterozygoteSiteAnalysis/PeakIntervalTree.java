package heterozygoteSiteAnalysis;

import jxl.Cell;
import jxl.Sheet;
import jxl.Workbook;
import jxl.read.biff.BiffException;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;


public class PeakIntervalTree {
    private Workbook workbook;
    private Sheet peakResultSht;

    public PeakIntervalTree(String resultExcel) {
        File peakCallingExcel = new File(resultExcel);
        try {
            this.workbook = Workbook.getWorkbook(peakCallingExcel);
            this.peakResultSht = workbook.getSheet("peak");
        } catch (IOException | BiffException xlse) {
            xlse.printStackTrace();
        }
    }

    public HashMap<String, HashMap<String, IntervalTree>> peakTrees() {
        HashMap<String, HashMap<String, IntervalTree>> treeMap = this.buildPeakIntervalTree();
        cleanup();

        return treeMap;
    }

    /**
     * get field name index of the peak calling result file
     */
    private HashMap<String, Integer> getFieldIndex() {
        HashMap<String, Integer> fieldIdx = new HashMap<String, Integer>();
        int totalCols = this.peakResultSht.getColumns();

        // the first row is field names
        for (int i = 0; i < totalCols; i++) {
            Cell cell = this.peakResultSht.getCell(0, i);
            String field = cell.getContents();
            if (field.equalsIgnoreCase("chr"))
                fieldIdx.put("chr", i);
            else if (field.equalsIgnoreCase("chromstart"))
                fieldIdx.put("chromstart", i);
            else if (field.equalsIgnoreCase("chromend"))
                fieldIdx.put("chromend", i);
            else if (field.equalsIgnoreCase("strand"))
                fieldIdx.put("strand", i);
            else if (field.equalsIgnoreCase("blockcount"))
                fieldIdx.put("blockcount", i);
            else if (field.equalsIgnoreCase("blocksizes"))
                fieldIdx.put("blocksizes", i);
            else if (field.equalsIgnoreCase("blockstarts"))
                fieldIdx.put("blockstarts", i);
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

        String chrNum, start, end, strand, blockCount, blockSize, blockStart;
        int totalRows = this.peakResultSht.getRows();
        for (int i = 1; i < totalRows; i++) {
            chrNum = this.peakResultSht.getCell(fieldIdx.get("chr"), i).getContents();
            if (!chrNum.equals(preChr)) {
                IntervalTree posStrandTree = new IntervalTree();
                IntervalTree negStrandTree = new IntervalTree();
                HashMap<String, IntervalTree> peakTree = new HashMap<>();
                peakTree.put("+", posStrandTree);
                peakTree.put("-", negStrandTree);

                peakTreeMap.put(chrNum, peakTree);
            }
            strand = this.peakResultSht.getCell(fieldIdx.get("strand"), i).getContents();
            start = this.peakResultSht.getCell(fieldIdx.get("chromstart"), i).getContents();
            end = this.peakResultSht.getCell(fieldIdx.get("chromend"), i).getContents();
            blockCount = this.peakResultSht.getCell(fieldIdx.get("blockcount"), i).getContents();
            blockSize = this.peakResultSht.getCell(fieldIdx.get("blocksizes"), i).getContents();
            blockStart = this.peakResultSht.getCell(fieldIdx.get("blockstarts"), i).getContents();
            if (Integer.parseInt(blockCount) > 1) {
                IntervalTreeNode[] nodes = this.multiBlock(start, blockSize, blockStart);
                for (IntervalTreeNode node : nodes) {
                    IntervalTree it = peakTreeMap.get(chrNum).get(strand);
                    it = it.insertNode(it, node);
                    peakTreeMap.get(chrNum).put(strand, it);
                }
            } else {
                IntervalTreeNode newNode = new IntervalTreeNode(Integer.parseInt(start), Integer.parseInt(end));
                IntervalTree it = peakTreeMap.get(chrNum).get(strand);
                it = it.insertNode(it, newNode);
                peakTreeMap.get(chrNum).put(strand, it);
            }

            preChr = chrNum;
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
    private IntervalTreeNode[] multiBlock(String start, String blockSize, String blockStarts) {
        String[] sizes = blockSize.split(",");
        String[] starts = blockStarts.split(",");
        int chrStart = Integer.parseInt(start);

        IntervalTreeNode[] nodes = new IntervalTreeNode[starts.length];
        for (int i = 0; i < starts.length; i++) {
            int intervalStart = chrStart + Integer.parseInt(starts[i]);
            int intervalEnd = intervalStart + Integer.parseInt(sizes[i]);
            IntervalTreeNode newNode = new IntervalTreeNode(intervalStart, intervalEnd);
            nodes[i] = newNode;
        }

        return nodes;
    }

    /**
     * release resource, garbage collection
     */
    private void cleanup() {
        if (this.workbook != null)
            workbook.close();
    }
}
