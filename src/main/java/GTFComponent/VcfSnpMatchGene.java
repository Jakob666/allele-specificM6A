package GTFComponent;

import heterozygoteSiteAnalysis.IntervalTree;
import heterozygoteSiteAnalysis.IntervalTreeNode;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * 将VCF文件中的record对应到各个基因中
 */
public class VcfSnpMatchGene {
    private String vcfFile;
    private HashMap<String, IntervalTree> gtfIntervalTree;
    // geneAlleleReads = {"geneId->geneName": {"major":[count1, count2,...], "minor": [count1, count2,...]}, ...}
    private HashMap<String, HashMap<String, ArrayList<Integer>>> geneAlleleReads = new HashMap<>();
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
            String line = "", chrNum, geneName, geneId;
            int position, ref, refReverse, alt, altReverse, majorCount, minorCount, majorStrand;
            int[] readsCount, majorAlleleRecord;
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
                    if (ref+refReverse==0 || alt+altReverse==0)
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
                    else {
                        gitn = (GTFIntervalTreeNode) itn;
                        geneName = gitn.geneName;
                        geneId = gitn.geneId;
                        this.renewGeneAlleleReads(geneId, geneName, majorCount, minorCount);
                        majorAlleleRecord = this.geneMajorAlleleStrand.getOrDefault(geneId, new int[2]);
                        majorAlleleRecord[majorStrand] = majorAlleleRecord[majorStrand] + 1;
                        this.geneMajorAlleleStrand.put(geneId, majorAlleleRecord);
                    }
                    // 基因之间可能存在重叠的情况
                    while (itn.leftChild != null && itn.leftChild.intervalEnd >= position ) {
                        gitn = (GTFIntervalTreeNode) itn.leftChild;
                        geneName = gitn.geneName;
                        geneId = gitn.geneId;
                        this.renewGeneAlleleReads(geneId, geneName, majorCount, minorCount);
                        majorAlleleRecord = this.geneMajorAlleleStrand.getOrDefault(geneId, new int[2]);
                        majorAlleleRecord[majorStrand] = majorAlleleRecord[majorStrand] + 1;
                        this.geneMajorAlleleStrand.put(geneId, majorAlleleRecord);
                    }
                    while (itn.rightChild != null && itn.rightChild.intervalStart <= position) {
                        gitn = (GTFIntervalTreeNode) itn.rightChild;
                        geneName = gitn.geneName;
                        geneId = gitn.geneId;
                        this.renewGeneAlleleReads(geneId, geneName, majorCount, minorCount);
                        majorAlleleRecord = this.geneMajorAlleleStrand.getOrDefault(geneId, new int[2]);
                        majorAlleleRecord[majorStrand] = majorAlleleRecord[majorStrand] + 1;
                        this.geneMajorAlleleStrand.put(geneId, majorAlleleRecord);
                    }
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

    private int[] parseDp4(String dp4String) {
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

    private void renewGeneAlleleReads(String geneId, String geneName, int majorCount, int minorCount) {
        String label = String.join("->", new String[]{geneId, geneName});
        HashMap<String, ArrayList<Integer>> geneSnp = this.geneAlleleReads.getOrDefault(label, new HashMap<>());
        ArrayList<Integer> majorAllele = geneSnp.getOrDefault("major", new ArrayList<>());
        ArrayList<Integer> minorAllele = geneSnp.getOrDefault("minor", new ArrayList<>());
        majorAllele.add(majorCount);
        minorAllele.add(minorCount);
        geneSnp.put("major", majorAllele);
        geneSnp.put("minor", minorAllele);
        this.geneAlleleReads.put(label, geneSnp);
    }

    public HashMap<String, HashMap<String, ArrayList<Integer>>> getGeneAlleleReads() {
        return this.geneAlleleReads;
    }

    public HashMap<String, int[]> getGeneMajorAlleleStrand() {
        return this.geneMajorAlleleStrand;
    }

}
