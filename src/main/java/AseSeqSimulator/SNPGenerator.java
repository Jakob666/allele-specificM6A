package AseSeqSimulator;

import org.apache.commons.math3.distribution.UniformRealDistribution;

import java.util.*;

/**
 * generate SNP on simulated sequencing reads
 */
public class SNPGenerator {
    private ArrayList<String> bases = new ArrayList<String>(Arrays.asList("A", "C", "T", "G"));
    private ArrayList<Integer> mutSiteNum = new ArrayList<>(Arrays.asList(1, 2, 3, 4, 5));
    private HashMap<String, HashSet<Integer>> mutGenePosition;
    private HashMap<String, String> originExonSequence = new HashMap<>();
    private HashMap<String, Double> refProportion = new HashMap<>();
    private HashMap<String, String> mutatedExonSeqence = new HashMap<>();
    private HashMap<String, LinkedList<Gene>> selectedGenes;
    private int mutateGeneNum;

    /**
     * Constructor
     */
    public SNPGenerator(HashMap<String, LinkedList<Gene>> selectedGenes, int mutateGeneNum) {
        this.selectedGenes = selectedGenes;
        this.mutateGeneNum = mutateGeneNum;
        this.randomMutateGene();
//        this.refAltProportion();
        this.mutatedExonSeqenceuence();
    }

    /**
     * random select mutated genes from all selected genes, and select serveral SNP site on its exon sequence
     */
    private void randomMutateGene() {
        ArrayList<Gene> totalGene = new ArrayList<>();
        HashMap<String, HashSet<Integer>> mutatedGene = new HashMap<>();
        for (LinkedList<Gene> chrGenes: selectedGenes.values()) {
            totalGene.addAll(chrGenes);
        }
        for (int i = 0; i < this.mutateGeneNum; i++) {
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
            }
        }

        this.mutGenePosition = mutatedGene;
    }

    /**
     * proportion of reference and alternative reads for each gene
     */
    private void refAltProportion() {
        Set<String> mutGeneIds = this.mutGenePosition.keySet();
        UniformRealDistribution urd = new UniformRealDistribution(0, 1);
        double refProportion;
        for (String geneId: mutGeneIds) {
            refProportion = urd.sample();
            this.refProportion.put(geneId, refProportion);
        }
    }

    /**
     * measure the mutated nucleotide sequence of each mutate gene
     */
    private void mutatedExonSeqenceuence() {
        for (String geneId: this.mutGenePosition.keySet()) {
            HashSet<Integer> mutPositions = this.mutGenePosition.get(geneId);
            this.mutatedExonSeqence.put(geneId, this.originExonSequence.get(geneId));
            for (Integer position : mutPositions) {
                this.singleSiteMutation(geneId, position);
            }
        }    
    }

    /**
     * site mutation on sequencing read
     */
    private void singleSiteMutation(String geneId, int mutPosition) {
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
    }
    
    public HashMap<String, HashSet<Integer>> getMutGenePosition() {
        return this.mutGenePosition;
    }
    
    public HashMap<String, String> getOriginExonSequence() {
        return this.originExonSequence;
    }
    
    public HashMap<String, String> getmutatedExonSeqence() {
        return this.mutatedExonSeqence;
    }

    public HashMap<String, Double> getRefProportion() {
        return this.refProportion;
    }

}
