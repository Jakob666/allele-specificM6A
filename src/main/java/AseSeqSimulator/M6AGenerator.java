package AseSeqSimulator;

import GTFComponent.ElementRecord;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.UniformIntegerDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;

import java.io.*;
import java.util.*;

/**
 * 给每个基因的mRNA上生成随机数目个m6A修饰位点
 */
public class M6AGenerator {
    private HashMap<String, HashSet<Integer>> mutGenePosition;
    private NormalDistribution m6aFrequency = new NormalDistribution(6, 2);

    /**
     * 构造方法
     * @param mutGenePosition 突变基因中SNP位点的信息
     */
    public M6AGenerator(HashMap<String, HashSet<Integer>> mutGenePosition) {
        this.mutGenePosition = mutGenePosition;
    }

    /**
     * 随机生成基因m6A修饰位点，该方法为公有方法，供外界调用，返回m6A修饰位点外显子的位置及其基因组的位置
     * @param m6aModifyGene 需要生成m6A修饰位点的Gene对象
     * @return m6A修饰位点
     */
    public HashMap<Integer, Integer> generateM6aSites(Gene m6aModifyGene) {
        String modifiedGeneIds = m6aModifyGene.getGeneId();
        int geneExonSeqLength = m6aModifyGene.getExonSeq().length();
        UniformIntegerDistribution uid;
        HashSet<Integer> exonMutations, geneM6aModifySites = new HashSet<>();
        int head, tail, modifyNum, modifyPosition;

        exonMutations = this.mutGenePosition.getOrDefault(modifiedGeneIds, null);

        // 获取突变位点的范围，用于生成均匀分布的抽样范围
        head = (exonMutations != null)? Collections.min(exonMutations) : 1;
        head = (head - 100 > 0)? (head - 100): 1;
        tail = (exonMutations != null)? Collections.max(exonMutations) : geneExonSeqLength;
        tail = (tail + 100 <= geneExonSeqLength)? (tail + 100): geneExonSeqLength - 1;
        uid = new UniformIntegerDistribution(head, tail);
        // 随机选取mRNA上m6a修饰位点数目并生成m6a修饰位点(外显子序列的位点，并非基因组位点)
        modifyNum = Math.abs((int) this.m6aFrequency.sample());
        for (int i = 0; i < modifyNum; ) {
            modifyPosition = uid.sample();
            if (exonMutations != null &&!exonMutations.contains(modifyPosition) && !geneM6aModifySites.contains(modifyPosition)) {
                geneM6aModifySites.add(modifyPosition);
                i++;
            } else if (exonMutations == null && !geneM6aModifySites.contains(modifyPosition)) {
                geneM6aModifySites.add(modifyPosition);
                i++;
            }
        }

        int genomePosition;
        HashMap<Integer, Integer> m6aModificationSites = new HashMap<>();
        for (Integer m6aSite: geneM6aModifySites) {
            genomePosition = this.m6aGenomePosition(m6aModifyGene, m6aSite);
            m6aModificationSites.put(m6aSite, genomePosition);
        }
        // 释放内存
        geneM6aModifySites = null;
        exonMutations = null;
        uid = null;

        return m6aModificationSites;
    }

    /**
     * 将随机生成的基因的m6A位点写入到文件中
     * @param simulatedM6aSites 记录每个基因修饰位点的哈希表，HashMap<String, HashMap<Integer, Integer>> geneId: exon position: genome position
     * @param outputFile 输出文件
     */
    public HashMap<String, HashMap<Integer, Double>> storeGeneM6aSites(HashMap<String, HashMap<Integer, Integer>> simulatedM6aSites, File outputFile) {
        HashMap<String, HashMap<Integer, Double>> asmRatio = new HashMap<>();
        UniformRealDistribution urd = new UniformRealDistribution(0.6, 0.8);
        BufferedWriter bfw = null;
        try {
            bfw = new BufferedWriter(
                    new OutputStreamWriter(new FileOutputStream(outputFile))
            );
            bfw.write("geneId\texonPosition\tgenomePosition\tmajorASMRatio\n");
            HashMap<Integer, Integer> m6aSites;
            Integer genomePosition;
            for (String geneId: simulatedM6aSites.keySet()) {
                m6aSites = simulatedM6aSites.get(geneId);
                for (Integer exonPosition: m6aSites.keySet()) {
                    double randNum = Math.random(), ase = 0.5;
                    if (randNum > 0.5)
                        ase = urd.sample();
                    genomePosition = m6aSites.get(exonPosition);
                    bfw.write(geneId + "\t" + exonPosition + "\t" + genomePosition + "\t" + ase);
                    bfw.newLine();
                    HashMap<Integer, Double> geneAsm = asmRatio.getOrDefault(geneId, new HashMap<>());
                    geneAsm.put(exonPosition, ase);
                    asmRatio.put(geneId, geneAsm);
                }
            }
            bfw.close();
        } catch (IOException ie) {
            ie.printStackTrace();
        } finally {
            if (bfw != null) {
                try {
                    bfw.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }

        return asmRatio;
    }

    /**
     * 通过m6A在基因外显子区域的位置，得到其基因组上的位置
     * @param gene Gene对象
     * @param exonM6aPos 外显子序列m6A突变位置
     * @return 基因组上的对应位置
     */
    private int m6aGenomePosition(Gene gene, int exonM6aPos) {
        ElementRecord exon = gene.getExonList();
        int genomePosition;
        int length = 0, distance = exonM6aPos;
        while (exon != null) {
            length = length + exon.getElementEnd() - exon.getElementStart() + 1;
            if (length >= exonM6aPos)
                break;
            exon = exon.getNextElement();
            distance = exonM6aPos - length;
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
}
