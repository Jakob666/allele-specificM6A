package GTFComponent;

import heterozygoteSiteAnalysis.IntervalTree;
import heterozygoteSiteAnalysis.IntervalTreeNode;

import java.io.*;
import java.util.HashMap;

/**
 * 将VCF文件中的record对应到各个基因中
 */
public class VcfSnpMatchGene {
    protected String vcfFile;
    // 每条染色体的GTF区间树
    protected HashMap<String, IntervalTree> gtfIntervalTree;
    // geneAlleleReads = {"geneId->geneName": {pos1: "major":count1, "minor": count1}, ...}
    private HashMap<String, HashMap<Integer, HashMap<String, Integer>>> geneAlleleReads = new HashMap<>();
    private HashMap<String, int[]> geneMajorAlleleStrand = new HashMap<>();

    public VcfSnpMatchGene(String vcfFile, String gtfFile) {
        this.vcfFile = vcfFile;
        GTFIntervalTree git = new GTFIntervalTree(gtfFile);
        git.parseGTFFile();
        this.gtfIntervalTree = git.getGtfIntervalTrees();
    }

    public void parseVcfFile() {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.vcfFile))));
            String line = "", chrNum, refNc, altNc;
            int position, ref, refReverse, alt, altReverse, majorCount, minorCount, majorStrand;
            int[] readsCount;
            IntervalTree it;
            IntervalTreeNode itn;
            String[] info;
            GTFIntervalTreeNode gitn;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    refNc = info[3];
                    altNc = info[4];
                    if (refNc.length() > 1 || altNc.length() > 1)
                        continue;
                    chrNum = info[0];
                    it = this.gtfIntervalTree.getOrDefault(chrNum, null);
                    if (it == null)
                        continue;
                    position = Integer.parseInt(info[1]);
                    readsCount = this.parseDp4(info[7]);
                    ref = readsCount[0];
                    refReverse = readsCount[1];
                    alt = readsCount[2];
                    altReverse = readsCount[3];
                    if (alt+altReverse==0)  // ref+refReverse==0 ||
                        continue;
                    if (ref + refReverse >= alt + altReverse) {
                        majorCount = ref + refReverse;
                        minorCount = alt + altReverse;
                        majorStrand = 0;
                    } else {
                        majorCount = alt + altReverse;
                        minorCount = ref + refReverse;
                        majorStrand = 1;
                    }

                    itn = it.search(it.root, position);
                    if (itn == null)
                        continue;
                    gitn = (GTFIntervalTreeNode) itn;
                    this.recursiveSearch(gitn, majorCount, minorCount, majorStrand, position);
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

    protected int[] parseDp4(String dp4String) {
        String[] info = dp4String.split(";"), readsCount = new String[4];
        for (String s: info) {
            if (s.startsWith("DP4")) {
                readsCount = s.split("=")[1].split(",");
            }
        }
        int[] reads = new int[4];
        for (int i=0; i<readsCount.length; i++) {
            reads[i] = Integer.parseInt(readsCount[i]);
        }

        return reads;
    }

    private void renewGeneAlleleReads(String geneId, String geneName, int mutPosition, int majorCount, int minorCount) {
        String label = String.join("->", new String[]{geneId, geneName});
        HashMap<Integer, HashMap<String, Integer>> geneSnp = this.geneAlleleReads.getOrDefault(label, new HashMap<>());
        HashMap<String, Integer> alleleReads = new HashMap<>();
        alleleReads.put("major", majorCount);
        alleleReads.put("minor", minorCount);

        geneSnp.put(mutPosition, alleleReads);
        this.geneAlleleReads.put(label, geneSnp);
    }

    public HashMap<String, HashMap<Integer, HashMap<String, Integer>>> getGeneAlleleReads() {
        return this.geneAlleleReads;
    }

    public HashMap<String, int[]> getGeneMajorAlleleStrand() {
        return this.geneMajorAlleleStrand;
    }

    public void recursiveSearch(GTFIntervalTreeNode gitn, int majorCount, int minorCount, int majorStrand, int position) {
        String geneName = gitn.geneName;
        String geneId = gitn.geneId;
        this.renewGeneAlleleReads(geneId, geneName, position, majorCount, minorCount);
        int[] majorAlleleRecord;
        majorAlleleRecord = this.geneMajorAlleleStrand.getOrDefault(geneId, new int[2]);
        majorAlleleRecord[majorStrand] = majorAlleleRecord[majorStrand] + 1;
        this.geneMajorAlleleStrand.put(geneId, majorAlleleRecord);

        if (gitn.rightChild != null) {
            GTFIntervalTreeNode rc = (GTFIntervalTreeNode) gitn.rightChild;
            if (rc.intervalStart <= position && position <= rc.intervalEnd)
                this.recursiveSearch(rc, majorCount, minorCount, majorStrand, position);
        }
        if (gitn.leftChild != null) {
            GTFIntervalTreeNode lc = (GTFIntervalTreeNode) gitn.leftChild;
            if (lc.intervalStart <= position && position <= lc.intervalEnd)
                this.recursiveSearch(lc, majorCount, minorCount, majorStrand, position);
        }
    }

}
