package AseSeqSimulator;

import GTFComponent.ElementRecord;
import heterozygoteSiteAnalysis.IntervalTree;
import heterozygoteSiteAnalysis.IntervalTreeNode;
import org.apache.commons.math3.distribution.UniformRealDistribution;

import java.util.*;

/**
 * generate SNP on simulated sequencing reads
 */
public class SNPGenerator {
    private ArrayList<String> bases = new ArrayList<String>(Arrays.asList("A", "C", "T", "G"));
    private ArrayList<Integer> mutSiteNum;
    private HashMap<String, HashSet<Integer>> mutGenePosition;
    private HashMap<String, String> originExonSequence = new HashMap<>();
    private HashMap<String, String> mutatedExonSeqence = new HashMap<>();
    private HashMap<String, LinkedList<Gene>> selectedGenes;
    private int mutateGeneNum;
    private String vcfFile;

    /**
     * Constructor
     */
    public SNPGenerator(HashMap<String, LinkedList<Gene>> selectedGenes, int mutateGeneNum, String vcfFile, int maxMutNum) {
        this.selectedGenes = selectedGenes;
        this.mutateGeneNum = mutateGeneNum;
        this.setMutSiteNum(maxMutNum);
        this.vcfFile = vcfFile;
        if (this.vcfFile == null) {
            this.randomMutateGene();
            this.randomMutatedExonSequence();
        } else
            this.vcfFileMutateGene();

    }

    /**
     * set the maximum mutation site number on fragments
     * @param num maximum mutation site number
     */
    private void setMutSiteNum(int num) {
        ArrayList<Integer> mutNum = new ArrayList<>();
        for (int i = 1; i <= num; i++) {
            mutNum.add(i);
        }

        this.mutSiteNum = mutNum;
    }

    /**
     * randomly select mutated genes from all selected genes, and randomly select several SNP site on its exon sequence
     */
    private void randomMutateGene() {
        ArrayList<Gene> totalGene = new ArrayList<>();
        HashMap<String, HashSet<Integer>> mutatedGene = new HashMap<>();
        for (LinkedList<Gene> chrGenes: selectedGenes.values()) {
            totalGene.addAll(chrGenes);
        }
        int i = 0;
        while (i < this.mutateGeneNum) {
            Collections.shuffle(totalGene);
            Gene mutGene = totalGene.get(0);
            if (!mutatedGene.keySet().contains(mutGene.getGeneId())) {
                String exonSeq = mutGene.getExonSeq();
                this.originExonSequence.put(mutGene.getGeneId(), exonSeq);

                UniformRealDistribution urd = new UniformRealDistribution(0, exonSeq.length()-1);
                Collections.shuffle(this.mutSiteNum);
                int mutNum = this.mutSiteNum.get(0);
                int order = 0;
                HashSet<Integer> geneMutPosition = new HashSet<>();
                while (order < mutNum){
                    Integer mutPosition = (int) urd.sample();
                    if (geneMutPosition.contains(mutPosition))
                        continue;
                    geneMutPosition.add(mutPosition);
                    order ++;
                }
                mutatedGene.put(mutGene.getGeneId(), geneMutPosition);
                i ++;
            }
        }

        this.mutGenePosition = mutatedGene;
    }

    /**
     * get SNP site from user VCF file
     */
    private void vcfFileMutateGene() {
        HashMap<String, HashSet<Integer>> mutatedGenePos = new HashMap<>();
        HashMap<String, IntervalTree> chrExonTrees = this.generateChrExonTrees();
        HashMap<String, LinkedList<VcfRecord>> chrVcfs = this.getVcfFromFile();

        // get common chromosome of select genes and vcfFiles
        Set<String> commonChrs = chrExonTrees.keySet();
        commonChrs.retainAll(chrVcfs.keySet());
        if (commonChrs.size() == 0) {
            System.out.println("no common chromosomes");
            System.exit(2);
        }

        for (String chr: commonChrs) {
            LinkedList<VcfRecord> vcfRecords = chrVcfs.get(chr);
            IntervalTree it = chrExonTrees.get(chr);

            for (VcfRecord record: vcfRecords) {
                ExonTreeNode etn = (ExonTreeNode) it.search(it.root, record.vcfSite);
                if (etn != null) {
                    String geneId = etn.geneName;
                    HashSet<Integer> geneMutSites = mutatedGenePos.getOrDefault(geneId, new HashSet<>());
                    geneMutSites.add(record.vcfSite-etn.intervalStart);
                    mutatedGenePos.put(geneId, geneMutSites);
                    // renew the mutated exon sequence
                    this.mutateViaVcfRecord(chr, geneId, record.vcfSite, record.alt);
                }
            }
        }

        this.mutGenePosition = mutatedGenePos;
    }

    /**
     * use gene exons to establish interval trees
     * @return Interval trees
     */
    private HashMap<String, IntervalTree> generateChrExonTrees() {
        HashMap<String, IntervalTree> chrExonTrees = new HashMap<>();
        for (String chrNum: this.selectedGenes.keySet()) {
            IntervalTree exonTree = new IntervalTree();
            LinkedList<Gene> chrGenes = this.selectedGenes.get(chrNum);
            for (Gene gene: chrGenes) {
                LinkedList<ElementRecord> exons = gene.getExonList();
                for (ElementRecord exon: exons) {
                    int exonStart = exon.getElementStart();
                    int exonEnd = exon.getElementEnd();
                    ExonTreeNode etn = new ExonTreeNode(exonStart, exonEnd, gene.getGeneId(), chrNum);
                    exonTree = exonTree.insertNode(exonTree, etn);
                }
            }
            chrExonTrees.put(chrNum, exonTree);
        }
        return chrExonTrees;
    }

    /**
     * get vcf sites of each chromosome from VCF file
     * @return HashMap
     */
    private HashMap<String, LinkedList<VcfRecord>> getVcfFromFile() {
        VcfReader vcfReader = new VcfReader(this.vcfFile);
        vcfReader.getChrVcfSites();
        HashMap<String, LinkedList<VcfRecord>> chrVcfs = vcfReader.getChrVcfs();

        return chrVcfs;
    }

    /**
     * randomly mutate nucleotide sequence of each mutate gene
     */
    private void randomMutatedExonSequence() {
        for (String geneId: this.mutGenePosition.keySet()) {
            HashSet<Integer> mutPositions = this.mutGenePosition.get(geneId);
            this.mutatedExonSeqence.put(geneId, this.originExonSequence.get(geneId));
            for (Integer position : mutPositions) {
                this.singleSiteMutation(geneId, position);
            }
        }
    }

    /**
     * mutate the exonic sequence according to the records in VCF file
     */
    private void mutateViaVcfRecord(String chrNum, String geneId, int mutatePosition, String alt) {
        LinkedList<Gene> chrGenes = this.selectedGenes.get(chrNum);

        for (Gene gene: chrGenes) {
            int idx = 0;
            if (!gene.getGeneId().equals(geneId))
                continue;
            LinkedList<ElementRecord> exons = gene.getExonList();
            String exonSeq = this.mutatedExonSeqence.getOrDefault(geneId, gene.getExonSeq());
            for (ElementRecord exon: exons) {
                int start = exon.getElementStart();
                int end = exon.getElementEnd();
                if (start > mutatePosition)
                    continue;
                else if (mutatePosition > end)
                    idx = idx + exon.getElementEnd() - exon.getElementStart() + 1;
                else
                    idx = idx + mutatePosition - start + 1;
            }
            String mutExonSeq = exonSeq.substring(0, idx) + alt + exonSeq.substring(idx+1);
            this.mutatedExonSeqence.put(geneId, mutExonSeq);
        }
    }

    /**
     * site mutation on sequencing read
     */
    private void singleSiteMutation(String geneId, int mutPosition) {
        String exonSeq = this.mutatedExonSeqence.get(geneId);
        String mutateSeq;
        String ref = Character.toString(exonSeq.charAt(mutPosition-1));
        String mut = ref;
        while (mut.equals(ref)) {
            Collections.shuffle(this.bases);
            mut = this.bases.get(0);
        }

        mutateSeq = exonSeq.substring(0, mutPosition) + mut + exonSeq.substring(mutPosition+1);
        this.mutatedExonSeqence.put(geneId, mutateSeq);
    }
    
    public HashMap<String, HashSet<Integer>> getMutGenePosition() {
        return this.mutGenePosition;
    }
    
    public HashMap<String, String> getOriginExonSequence() {
        return this.originExonSequence;
    }
    
    public HashMap<String, String> getMutatedExonSeqence() {
        return this.mutatedExonSeqence;
    }

}

class ExonTreeNode extends IntervalTreeNode {
    public String chrNum, geneName;
    public ExonTreeNode parent, leftChild, rightChild;

    public ExonTreeNode(int start, int end, String geneName, String chrNum) {
        super(start, end, start, end);
        this.geneName = geneName;
        this.chrNum = chrNum;
    }
}
