package GTFComponent;

import heterozygoteSiteAnalysis.DbsnpAnnotation;
import heterozygoteSiteAnalysis.IntervalTree;
import heterozygoteSiteAnalysis.IntervalTreeNode;

import java.io.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;

/**
 * match VCF file mutation record to corresponding gene
 */
public class VcfSnpMatchGene {
    protected String vcfFile;
    private int readsCoverageThreshold;
    // GTF interval tree for each chromosome
    private HashMap<String, IntervalTree> gtfIntervalTree;
    private HashMap<String, HashMap<String, IntervalTree>> geneExonIntervalTree;
    // geneAlleleReads = {"geneId->geneName": {pos1: [majorAllele:count, minorAllele: count1]}, ...}
    private HashMap<String, HashMap<Integer, String[]>> geneAlleleReads = new HashMap<>();

    public VcfSnpMatchGene(String vcfFile, String gtfFile, int readsCoverageThreshold) {
        this.vcfFile = vcfFile;
        GTFIntervalTree git = new GTFIntervalTree(gtfFile);
        git.parseGTFFile();
        GeneExonIntervalTree geit = new GeneExonIntervalTree(gtfFile);
        geit.generateExonTree();
        this.gtfIntervalTree = git.getGtfIntervalTrees();
        this.geneExonIntervalTree = geit.getGeneExonIntervalTree();
        this.readsCoverageThreshold = readsCoverageThreshold;
    }

    public HashMap<String, HashMap<Integer, String[]>> getGeneAlleleReads() {
        return this.geneAlleleReads;
    }

    /**
     * parse SNP calling output VCF format file, match each SNV site to corresponding gene
     */
    public void parseVcfFile(HashMap<String, LinkedList<DbsnpAnnotation.DIYNode>> dbsnpRecord) {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.vcfFile))));
            String line = "", chrNum, refNc, altNc, majorNc, minorNc;
            int position, allele1Count, allele2Count, majorAlleleCount, minorAlleleCount;
            boolean specialMutation;
            int[] readsCount;
            IntervalTree it;
            IntervalTreeNode itn;
            HashMap<String, IntervalTree> eit;
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
                    // whether SNV site contains multiple possible mutations
                    specialMutation = altNc.length() > 1 && altNc.contains(",");

                    chrNum = info[0];
                    it = this.gtfIntervalTree.getOrDefault(chrNum, null);
                    if (it == null)
                        continue;
                    eit = this.geneExonIntervalTree.getOrDefault(chrNum, null);
                    // the mutation not in dbsnp record
                    if (dbsnpRecord != null && !DbsnpAnnotation.getSearchRes(dbsnpRecord, chrNum, info[1]))
                        continue;
                    position = Integer.parseInt(info[1]);

                    if (info.length >= 10) {    // may contains FORMAT and sample
                        readsCount = this.parseAd(info[8], info[9]);   // parse AD(allele depth) from VCF record
                        if (readsCount != null) {
                            if (!specialMutation) {
                                allele1Count = readsCount[0];
                                allele2Count = readsCount[1];
                                if (allele1Count >= allele2Count) {
                                    majorAlleleCount = allele1Count;
                                    minorAlleleCount = allele2Count;
                                    majorNc = refNc;
                                    minorNc = altNc;
                                } else {
                                    majorAlleleCount = allele2Count;
                                    minorAlleleCount = allele1Count;
                                    majorNc = altNc;
                                    minorNc = refNc;
                                }
                            } else {
                                String[] possibleMutations = altNc.split(",");
                                int[] originReads = Arrays.copyOf(readsCount, readsCount.length);
                                Arrays.sort(readsCount);
                                majorAlleleCount = readsCount[readsCount.length - 1];
                                minorAlleleCount = readsCount[readsCount.length - 2];
                                int[] indexRes = findIndex(originReads, majorAlleleCount, minorAlleleCount);
                                int majorAlleleIdx = indexRes[0];
                                int minorAlleleIdx = indexRes[1];
                                majorNc = (majorAlleleIdx == 0)? refNc : possibleMutations[majorAlleleIdx - 1];
                                minorNc = (minorAlleleIdx == 0)? refNc : possibleMutations[minorAlleleIdx - 1];
                            }
                        } else {    // if FORMAT contains no AD information, try to parse DP4 of INFO column
                            readsCount = this.parseDp4(info[7]);
                            if (readsCount == null) {
                                System.out.println("invalid VCF record: " + line );
                                System.exit(2);
                            }
                            allele1Count = readsCount[0] + readsCount[1];
                            allele2Count = readsCount[2] + readsCount[3];
                            majorAlleleCount = (allele1Count >= allele2Count)? allele1Count: allele2Count;
                            minorAlleleCount = (allele1Count < allele2Count)? allele1Count: allele2Count;
                            majorNc = (allele1Count >= allele2Count)? refNc : altNc.split(",")[0];
                            minorNc = (allele1Count < allele2Count)? refNc: altNc.split(",")[0];
                        }
                    } else {  // contains INFO but no INFO and Sample
                        readsCount = this.parseDp4(info[7]);
                        if (readsCount == null) {
                            System.out.println("invalid VCF record: " + line );
                            System.exit(2);
                        }
                        allele1Count = readsCount[0] + readsCount[1];
                        allele2Count = readsCount[2] + readsCount[3];
                        majorAlleleCount = (allele1Count >= allele2Count)? allele1Count: allele2Count;
                        minorAlleleCount = (allele1Count < allele2Count)? allele1Count: allele2Count;
                        majorNc = (allele1Count >= allele2Count)? refNc : altNc.split(",")[0];
                        minorNc = (allele1Count < allele2Count)? refNc: altNc.split(",")[0];
                    }

                    if (majorAlleleCount < this.readsCoverageThreshold)
                        continue;

                    // locate the SNV site on particular gene using GTF interval tree
                    itn = it.search(it.root, position);
                    if (itn == null)
                        continue;
                    gitn = (GTFIntervalTreeNode) itn;
                    this.recursiveSearch(gitn, eit, majorNc, minorNc, majorAlleleCount, minorAlleleCount, position);
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
            this.gtfIntervalTree = null;
            this.geneExonIntervalTree = null;
        }
    }

    /**
     * get reference and alternative nucleotide from VCF format file
     * @param dp4String INFO column in VCF format file
     * @return int[] {ref, ref_reverse, alt, alt_reverse}
     */
    protected int[] parseDp4(String dp4String) {
        String[] info = dp4String.split(";"), readsCount = null;
        for (String s: info) {
            if (s.startsWith("DP4")) {
                readsCount = s.split("=")[1].split(",");
            }
        }
        if (readsCount == null)
            return null;
        int[] reads = new int[4];
        for (int i=0; i<readsCount.length; i++) {
            reads[i] = Integer.parseInt(readsCount[i]);
        }

        return reads;
    }

    /**
     * get reference and alternative nucleotide from VCF format file
     * @param format FORMAT information in VCF file
     * @param sampleInfo sample information in VCF file
     * @return int[] {ref, ref_reverse, alt, alt_reverse}
     */
    protected int[] parseAd(String format, String sampleInfo) {
        String[] formatContent = format.split(":");
        String[] sampleData = sampleInfo.split(":");
        int idx = -1;
        for (int i=0; i<formatContent.length; i++) {
            if (formatContent[i].equals("AD")) {
                idx = i;
                break;
            }
        }
        if (idx < 0)
            return null;
        String[] readsRecord = sampleData[idx].split(",");
        int[] readsCount = new int[readsRecord.length];
        for (int i=0; i<readsRecord.length; i++) {
            readsCount[i] = Integer.valueOf(readsRecord[i]);
        }
        formatContent = null;
        sampleData = null;
        readsRecord = null;
        return readsCount;
    }

    public int[] findIndex(int[] readsCount, int majorCount, int minorCount) {
        int majorIdx = -1, minorIdx = -1;
        for (int i=0; i<readsCount.length; i++) {
            if (readsCount[i] == majorCount && majorIdx < 0)
                majorIdx = i;
            if (readsCount[i] == minorCount && minorIdx < 0 && majorIdx != i)
                minorIdx = i;
        }
        return new int[] {majorIdx, minorIdx};
    }

    /**
     * renew geneAlleleReads
     * @param geneId gene id
     * @param geneName gene name
     * @param majorAllele major allele nucleotide
     * @param minorAllele minor allele nucleotide
     * @param mutPosition mutation position
     * @param majorAlleleCount major allele reads count
     * @param minorAlleleCount minor allele reads count
     */
    private void renewGeneAlleleReads(String geneId, String geneName, String majorAllele, String minorAllele,
                                      int mutPosition, int majorAlleleCount, int minorAlleleCount) {
        String label = String.join("->", new String[]{geneId, geneName});
        HashMap<Integer, String[]> geneSnp = this.geneAlleleReads.getOrDefault(label, new HashMap<>());
        String reference = String.join(":", new String[] {majorAllele, String.valueOf(majorAlleleCount)});
        String alternative = String.join(":", new String[] {minorAllele, String.valueOf(minorAlleleCount)});

        geneSnp.put(mutPosition, new String[] {reference, alternative});
        this.geneAlleleReads.put(label, geneSnp);
    }

    /**
     * recursive search
     * @param gitn GTFIntervalTreeNode instance on GTFIntervalTree
     * @param exonIntervalTree IntervalTree records genes' exon region
     * @param majorAllele major allele nucleotide
     * @param minorAllele minor allele nucleotide
     * @param majorAlleleCount major allele reads count
     * @param minorAlleleCount minor allele reads count
     * @param position mutation position
     */
    public void recursiveSearch(GTFIntervalTreeNode gitn, HashMap<String, IntervalTree> exonIntervalTree,
                                String majorAllele, String minorAllele,
                                int majorAlleleCount, int minorAlleleCount, int position) {
        String geneName = gitn.geneName;
        String geneId = gitn.geneId;
        IntervalTree intervalTree = exonIntervalTree.getOrDefault(geneId, null);

        if (intervalTree != null) {
            IntervalTreeNode itn = intervalTree.search(intervalTree.root, position);
            if (itn != null) {
                this.renewGeneAlleleReads(geneId, geneName, majorAllele, minorAllele,
                        position, majorAlleleCount, minorAlleleCount);
            }
        }

        if (gitn.rightChild != null) {
            GTFIntervalTreeNode rc = (GTFIntervalTreeNode) gitn.rightChild;
            if (rc.intervalStart <= position && position <= rc.intervalEnd)
                this.recursiveSearch(rc, exonIntervalTree, majorAllele, minorAllele, majorAlleleCount,
                                     minorAlleleCount, position);
        }
        if (gitn.leftChild != null) {
            GTFIntervalTreeNode lc = (GTFIntervalTreeNode) gitn.leftChild;
            if (lc.intervalStart <= position && position <= lc.intervalEnd)
                this.recursiveSearch(lc, exonIntervalTree, majorAllele, minorAllele, majorAlleleCount,
                                     minorAlleleCount, position);
        }
    }

}
