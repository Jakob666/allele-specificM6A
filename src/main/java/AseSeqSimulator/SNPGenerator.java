package AseSeqSimulator;

import GTFComponent.ElementRecord;
import heterozygoteSiteAnalysis.IntervalTree;
import org.apache.commons.math3.distribution.UniformIntegerDistribution;

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
    private HashMap<String, HashMap<Integer, String[]>> mutationRefAlt = new HashMap<>();
    private HashMap<String, LinkedList<Gene>> selectedGenes;
    private HashMap<String, HashMap<Integer, Integer>> mutGenomePosition = new HashMap<>();
    private double mutateProp;
    private String vcfFile;

    /**
     * Constructor
     */
    public SNPGenerator(HashMap<String, LinkedList<Gene>> selectedGenes, double mutateProportion, String vcfFile, int minMutNum, int maxMutNum) {
        this.selectedGenes = selectedGenes;
        this.mutateProp = mutateProportion;
        this.setMutSiteNum(minMutNum, maxMutNum);
        this.vcfFile = vcfFile;
        // if no VCF file support, randomly generate SNP on exon sequence
        if (this.vcfFile == null) {
            this.randomMutateGene();
            this.randomMutatedExonSequence();
        } else // vcf-based SNP generate
            this.vcfFileMutateGene();
    }

    /**
     * set the maximum mutation site number on fragments
     */
    private void setMutSiteNum(int minNum, int maxNum) {
        ArrayList<Integer> mutNum = new ArrayList<>();
        for (int i = minNum; i <= maxNum; i++) {
            mutNum.add(i);
        }

        this.mutSiteNum = mutNum;
    }

    /**
     * randomly select mutated genes from selected genes on each chromosome, and randomly select several
     * SNP site on its exon sequence
     */
    private void randomMutateGene() {
        HashMap<String, HashSet<Integer>> mutatedGene = new HashMap<>();
        for (LinkedList<Gene> chrGenes: selectedGenes.values()) {
            // measure the number of mutate gene on current chromosome and randomly select from list
            int mutGeneNum = (int) (this.mutateProp * chrGenes.size());
            Collections.shuffle(chrGenes);
            List<Gene> mutGenes = chrGenes.subList(0, mutGeneNum);

            for (Gene mutGene: mutGenes) {
                String exonSeq = mutGene.getExonSeq();
                this.originExonSequence.put(mutGene.getGeneId(), exonSeq);

                UniformIntegerDistribution uid = new UniformIntegerDistribution(1, exonSeq.length()-1);
                Collections.shuffle(this.mutSiteNum);
                int mutNum = this.mutSiteNum.get(0);
                int order = 0;
                HashSet<Integer> geneMutPosition = new HashSet<>();
                while (order < mutNum){
                    Integer mutPosition = uid.sample();
                    if (geneMutPosition.contains(mutPosition))
                        continue;
                    geneMutPosition.add(mutPosition);
                    HashMap<Integer, Integer> genomePos  = this.mutGenomePosition.getOrDefault(mutGene.getGeneId(), new HashMap<>());
                    genomePos.put(mutPosition, this.getGenomePosition(mutGene, mutPosition));
                    this.mutGenomePosition.put(mutGene.getGeneId(), genomePos);
                    order ++;
                }
                mutatedGene.put(mutGene.getGeneId(), geneMutPosition);
            }
        }

        this.mutGenePosition = mutatedGene;
    }

    /**
     * trans the exon sequence position to genome position
     * @param gene Gene instance
     * @param exonMutPos mutate exon sequence position
     * @return genome position
     */
    private int getGenomePosition(Gene gene, int exonMutPos) {
        ElementRecord exon = gene.getExonList();
        int genomePosition;
        int length = 0, distance = exonMutPos;
        while (exon != null) {
            length = length + exon.getElementEnd() - exon.getElementStart() + 1;
            if (length >= exonMutPos)
                break;
            exon = exon.getNextElement();
            distance = exonMutPos - length;
        }
        if (gene.getStrand().equals("+")) {
            int start = exon.getElementStart();
            genomePosition = start + distance;
        } else {
            int end = exon.getElementEnd();
            genomePosition = end - distance;
        }

        return genomePosition;
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
                ElementRecord exon = gene.getExonList();
                while (exon != null) {
                    int exonStart = exon.getElementStart();
                    int exonEnd = exon.getElementEnd();
                    ExonTreeNode etn = new ExonTreeNode(exonStart, exonEnd, gene.getGeneId(), chrNum);
                    exonTree = exonTree.insertNode(exonTree, etn);
                    exon = exon.getNextElement();
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
            ElementRecord exon = gene.getExonList();
            String exonSeq = this.mutatedExonSeqence.getOrDefault(geneId, gene.getExonSeq());
            while (exon != null) {
                int start = exon.getElementStart();
                int end = exon.getElementEnd();
                if (start > mutatePosition)
                    continue;
                else if (mutatePosition > end)
                    idx = idx + exon.getElementEnd() - exon.getElementStart() + 1;
                else
                    idx = idx + mutatePosition - start + 1;
                exon = exon.getNextElement();
            }
            String mutExonSeq = exonSeq.substring(0, idx) + alt + exonSeq.substring(idx+1);
            this.mutatedExonSeqence.put(geneId, mutExonSeq);
        }
    }

    /**
     * site mutation on sequencing read
     */
    private void singleSiteMutation(String geneId, int mutPosition) {
        HashMap<Integer, String[]> mutationType = this.mutationRefAlt.getOrDefault(geneId, new HashMap<>());
        String exonSeq = this.mutatedExonSeqence.get(geneId);
        String mutateSeq;
        String ref = Character.toString(exonSeq.charAt(mutPosition));
        String mut = ref;
        while (mut.equals(ref)) {
            Collections.shuffle(this.bases);
            mut = this.bases.get(0);
        }

        mutateSeq = exonSeq.substring(0, mutPosition) + mut + exonSeq.substring(mutPosition+1);
        this.mutatedExonSeqence.put(geneId, mutateSeq);
        mutationType.put(mutPosition, new String[]{ref, mut});
        this.mutationRefAlt.put(geneId, mutationType);
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

    public HashMap<String, HashMap<Integer, String[]>> getMutationRefAlt() {
        return this.mutationRefAlt;
    }

    public HashMap<String, HashMap<Integer, Integer>> getMutGenomePosition() {
        return this.mutGenomePosition;
    }

}
