package GTFComponent;

import heterozygoteSiteAnalysis.DbsnpAnnotation;
import heterozygoteSiteAnalysis.IntervalTree;
import heterozygoteSiteAnalysis.IntervalTreeNode;

import java.io.*;
import java.util.HashMap;
import java.util.LinkedList;

/**
 * match VCF file mutation record to corresponding gene
 */
public class VcfSnpMatchGene {
    protected String vcfFile;
    private int readsCoverageThreshold;
    // GTF interval tree for each chromosome
    public HashMap<String, IntervalTree> gtfIntervalTree;
    // geneAlleleReads = {"geneId->geneName": {pos1: [refAllele:count, altAllele: count1]}, ...}
    private HashMap<String, HashMap<Integer, String[]>> geneAlleleReads = new HashMap<>();
    private HashMap<String, HashMap<Integer, String>> geneMajorAlleleNucleotide = new HashMap<>();

    public VcfSnpMatchGene(String vcfFile, String gtfFile, int readsCoverageThreshold) {
        this.vcfFile = vcfFile;
        GTFIntervalTree git = new GTFIntervalTree(gtfFile);
        git.parseGTFFile();
        this.gtfIntervalTree = git.getGtfIntervalTrees();
        this.readsCoverageThreshold = readsCoverageThreshold;
    }

    public HashMap<String, HashMap<Integer, String[]>> getGeneAlleleReads() {
        return this.geneAlleleReads;
    }

    public HashMap<String, HashMap<Integer, String>> getGeneMajorAlleleNucleotide() {
        return this.geneMajorAlleleNucleotide;
    }

    /**
     * 解析SNP calling输出的VCF文件，使用其中的信息更新 geneAlleleReads和geneMajorAlleleStrand
     */
    public void parseVcfFile(HashMap<String, LinkedList<DbsnpAnnotation.DIYNode>> dbsnpRecord) {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.vcfFile))));
            String line = "", chrNum, refNc, altNc;
            int position, referenceCount, alternativeCount, refOrAlt;
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
                    // only take single nucleotide mutations into consideration
                    if (refNc.length() > 1 || altNc.length() > 1)
                        continue;
                    chrNum = info[0];
                    it = this.gtfIntervalTree.getOrDefault(chrNum, null);
                    if (it == null)
                        continue;
                    // the mutation not in dbsnp record
                    if (dbsnpRecord != null && !DbsnpAnnotation.getSearchRes(dbsnpRecord, chrNum, info[1]))
                        continue;
                    position = Integer.parseInt(info[1]);
                    readsCount = this.parseDp4(info[7]);
                    referenceCount = readsCount[0] + readsCount[1];
                    alternativeCount = readsCount[2] + readsCount[3];
                    // filter homo SNP sites
                    if (referenceCount == 0 || alternativeCount==0)
                        continue;
                    if (Math.max(referenceCount, alternativeCount) < this.readsCoverageThreshold)
                        continue;
                    // make sure the major allele equals to reference or alternative
                    if (referenceCount > alternativeCount)
                        refOrAlt = 0;
                    else if (referenceCount < alternativeCount)
                        refOrAlt = 1;
                    else
                        refOrAlt = Integer.MIN_VALUE;

                    // locate the SNV site on particular gene using GTF interval tree
                    itn = it.search(it.root, position);
                    if (itn == null)
                        continue;
                    gitn = (GTFIntervalTreeNode) itn;
                    this.recursiveSearch(gitn, refNc, altNc, referenceCount, alternativeCount, refOrAlt, position);
                    info = null;
                }
            }
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
     * 从VCF文件INFO列解析得到reference和alternative nucleotide的reads count数目
     * @param dp4String INFO列信息
     * @return int[] {ref, ref_reverse, alt, alt_reverse}
     */
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

    /**
     * 搜索到的结果更新 geneAlleleReads
     * @param geneId gene id
     * @param geneName gene name
     * @param referenceAllele reference allele对应的nucleotide
     * @param alternativeAllele alternative allele对应的nucleotide
     * @param mutPosition 突变位点
     * @param referenceCount reference allele覆盖的reads
     * @param alternativeCount alternative allele覆盖的reads
     */
    private void renewGeneAlleleReads(String geneId, String geneName, String referenceAllele, String alternativeAllele,
                                      int mutPosition, int referenceCount, int alternativeCount) {
        String label = String.join("->", new String[]{geneId, geneName});
        HashMap<Integer, String[]> geneSnp = this.geneAlleleReads.getOrDefault(label, new HashMap<>());
        String reference = String.join(":", new String[] {referenceAllele, Integer.toString(referenceCount)});
        String alternative = String.join(":", new String[] {alternativeAllele, Integer.toString(alternativeCount)});

        geneSnp.put(mutPosition, new String[] {reference, alternative});
        this.geneAlleleReads.put(label, geneSnp);
    }

    /**
     * 递归搜索
     * @param gitn GTFIntervalTree搜索返回的GTFIntervalTreeNode对象
     * @param referenceAllele reference allele对应的nucleotide
     * @param alternativeAllele alternative allele对应的nucleotide
     * @param referenceCount reference allele覆盖的reads
     * @param alternativeCount alternative allele覆盖的reads
     * @param refOrAlt major allele是reference还是alternative，分别使用0和1表示
     * @param position 突变位点
     */
    public void recursiveSearch(GTFIntervalTreeNode gitn, String referenceAllele, String alternativeAllele,
                                int referenceCount, int alternativeCount, int refOrAlt, int position) {
        String geneName = gitn.geneName;
        String geneId = gitn.geneId;
        this.renewGeneAlleleReads(geneId, geneName, referenceAllele, alternativeAllele, position, referenceCount, alternativeCount);
        HashMap<Integer, String> majorAlleleRecord;

        majorAlleleRecord = this.geneMajorAlleleNucleotide.getOrDefault(geneId, new HashMap<>());
        if (refOrAlt == 1)
            majorAlleleRecord.put(position, "alt");
        else
            majorAlleleRecord.put(position, "ref");
        this.geneMajorAlleleNucleotide.put(geneId, majorAlleleRecord);

        if (gitn.rightChild != null) {
            GTFIntervalTreeNode rc = (GTFIntervalTreeNode) gitn.rightChild;
            if (rc.intervalStart <= position && position <= rc.intervalEnd)
                this.recursiveSearch(rc, referenceAllele, alternativeAllele, referenceCount, alternativeCount, refOrAlt, position);
        }
        if (gitn.leftChild != null) {
            GTFIntervalTreeNode lc = (GTFIntervalTreeNode) gitn.leftChild;
            if (lc.intervalStart <= position && position <= lc.intervalEnd)
                this.recursiveSearch(lc, referenceAllele, alternativeAllele, referenceCount, alternativeCount, refOrAlt, position);
        }
    }

}
