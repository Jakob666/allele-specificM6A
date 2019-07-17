package AseSeqSimulator;

import GTFComponent.*;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;

import java.io.*;
import java.util.*;

public class ReadsGenerator {
    private HashMap<String, LinkedList<Gene>> ChrGeneMap = new HashMap<String, LinkedList<Gene>>();
    private GTFReader readGTF = new GTFReader();
    private int readLength, sequencingDepth;
    private double geneProp;
    private TwoBitParser twoBit;
    private long librarySize;
    private String gtfFile, outputDir;
    private HashMap<String, String> genesOriginExonSeq, genesMutatedExonSeq;
    private HashMap<String, HashSet<Integer>> geneMutatedPosition;
    private HashMap<String, HashMap<Integer, String[]>> geneMutationRefAlt;
    private HashMap<String, HashMap<Integer, Integer>> geneMutGenomePosition, geneM6aGenomePosition;
    private HashMap<String, HashMap<Integer, Double>> geneM6aAsmRatio;
    private HashMap<String, HashSet<String>> selectedGeneChr = new HashMap<>();
    private HashMap<String, Double> aseGeneMajorAlleleRatio = new HashMap<>();
    private SequencingError seqErrorModel;

    public ReadsGenerator(String gtfFile, double geneProp, String twoBitFile) {
        this.readGTF.readFromFile(gtfFile);
        System.out.println("read GTF complete");
        this.geneProp = geneProp;
        this.gtfFile = gtfFile;
        try {
            this.twoBit = new TwoBitParser(new File(twoBitFile));
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void setLibrarySize(long librarySize) {
        this.librarySize = librarySize;
    }

    private void setSequencingDepth(int sequencingDepth) {
        this.sequencingDepth = sequencingDepth;
    }

    private void setReadLength(int readLength) {
        this.readLength = readLength;
    }

    public HashMap<String, LinkedList<Gene>> getChrGeneMap() {
        return this.ChrGeneMap;
    }

    /**
     * 每条染色体上选取一定数目的基因，这些基因将被存放在 chrGeneMap 变量中。
     * @param overlapGene 基因间是否有允许重叠。true为允许有重叠；false表示没有重叠。
     */
    private void selectGene(boolean overlapGene) {
        HashMap<String, ChromosomeRecord> chrMap = this.readGTF.getChromosomeMap();
        try {
            String randomGeneId;
            HashMap<String, GeneRecord> genesOnChr;
            ArrayList<String> geneIds;
            LinkedList<Gene> geneList;
            LinkedList<int[]> transcriptRegion;

            for (String chr : chrMap.keySet()) {
                HashSet<String> chrGenes = new HashSet<>();
                transcriptRegion = new LinkedList<int[]>();
                // 2bit文件中不存在chrMT，线粒体
                if (chr.equals("MT"))
                    continue;
                ChromosomeRecord chromosomeRecord = chrMap.get(chr);
                genesOnChr = chromosomeRecord.getChrGenes();

                // 将染色体上的基因随机打乱，之后从中随机抽取
                geneIds = new ArrayList<String>(genesOnChr.keySet());
                int geneNum = (int) (this.geneProp * genesOnChr.size());
                System.out.println(chr+":"+geneNum);
                geneList = new LinkedList<>();
                int order = 0, exonLength = 0, transStart, transEnd;
                boolean overlap = false;
                Collections.shuffle(geneIds);
                for (int rank = 0; order < geneNum && rank < geneIds.size(); rank++) {
                    randomGeneId = geneIds.get(rank);
                    GeneRecord gr = genesOnChr.get(randomGeneId);
                    TranscriptRecord longestTranscript = CommonMethod.findLongestTranscript(gr);
                    HashMap<String, ElementRecord> elements = longestTranscript.getElementList();
                    ElementRecord exons = elements.getOrDefault("exon", null);
                    while (exons != null) {
                        int start = exons.getElementStart();
                        int end = exons.getElementEnd();
                        exonLength += end - start + 1;
                        exons = exons.getNextElement();
                    }

                    // 如果exon区域过短则直接跳过当前基因
                    if (exonLength <= 1000)
                        continue;
                    transStart = longestTranscript.getTranscriptStart();
                    transEnd = longestTranscript.getTranscriptEnd();

                    if (!overlapGene) { // 如果要求基因之间没有重叠，则判断当前基因与已选的基因之间是否存在重叠区域
                        for (int[] region: transcriptRegion) {
                            int start_max = Math.max(region[0], transStart);
                            int end_min = Math.min(region[1], transEnd);
                            if (start_max <= end_min) {
                                overlap = true;
                                break;
                            }
                        }
                        if (overlap)
                            continue;
                    }

                    transcriptRegion.add(new int[]{transStart, transEnd});
                    Gene gene = new Gene(gr.getGeneId(), gr.getGeneName(), gr.getStrand(), chr);
                    gene.setLongestTranscriptRecord(longestTranscript);
                    // 获取某个基因外显子最长的转录本，通过 Gene 对象的splicing方法
                    this.twoBit.setCurrentSequence("chr"+chr);
                    gene.splicing(this.twoBit);
                    geneList.add(gene);
                    twoBit.close();
                    chrGenes.add(gene.getGeneId());
//                    System.out.println(order + ":" + gene.getGeneId()+"->"+longestTranscript.getTranscriptId());
                    order ++;
                }
                this.selectedGeneChr.put(chr, chrGenes);
                this.ChrGeneMap.put(chr, geneList);
            }
            // 释放内存
            transcriptRegion = null;
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        // 释放内存
        chrMap = null;
        this.readGTF = null;
    }

    /**
     * 对选中的基因随机生成SNP位点
     * @param vcfFile 可以指定某个VCF文件生成SNP位点
     * @param minimumMut 外显子序列最小突变数目
     * @param maximumMut 外显子序列最大突变数目
     * @param mutProp 突变基因占基因总数的比例
     */
    private void randomMutateGene(String vcfFile, int minimumMut, int maximumMut, double mutProp) {
        SNPGenerator snpGenerator = new SNPGenerator(this.ChrGeneMap, mutProp, vcfFile, minimumMut, maximumMut);
        snpGenerator.generateSNP();
        this.genesOriginExonSeq = snpGenerator.getOriginExonSequence();     // 基因未突变的外显子序列
        this.genesMutatedExonSeq = snpGenerator.getMutatedExonSeqence();    // 突变基因的外显子序列
        this.geneMutatedPosition = snpGenerator.getMutGenePosition();       // 突变基因的突变位点在外显子上的位置
        this.geneMutationRefAlt = snpGenerator.getMutationRefAlt();
        this.geneMutGenomePosition = snpGenerator.getMutGenomePosition();   // 突变基因每个外显子突变位点对应的基因组位置
        snpGenerator = null;
    }

    /**
     * 设置具有SNP位点的基因 major allele和 minor allele的比例
     */
    private void mutatedGeneAseRatio() {
        UniformRealDistribution urd = new UniformRealDistribution(0.60, 0.80);
        // 获取到突变基因的 geneID集合
        Set<String> mutGeneIds = this.geneMutatedPosition.keySet();
        double ref, randNum;
        for (String mutGeneId: mutGeneIds) {
            randNum = Math.random();
            if (randNum > 0.5)
                this.aseGeneMajorAlleleRatio.put(mutGeneId, 0.5);
            else {
                ref = urd.sample();
                this.aseGeneMajorAlleleRatio.put(mutGeneId, ref);
            }
        }
    }

    /**
     * 依据随机选取的基因生成一个简化的GTF文件.
     * @param outputFile 简化的GTF文件输出路径
     */
    private void createSimpleGTF(String outputFile) {
        BufferedReader bfr = null;
        BufferedWriter bfw = null;
        try {
            bfr = new BufferedReader(
                    new InputStreamReader(new FileInputStream(new File(this.gtfFile)))
            );
            bfw = new BufferedWriter(
                    new OutputStreamWriter(new FileOutputStream(new File(outputFile)))
            );

            Set<String> chrs = this.ChrGeneMap.keySet();
            Set<String> selectGeneList = new HashSet<>(), selectedTranscriptList = new HashSet<>();
            for (String chr: chrs) {
                LinkedList<Gene> chrSelectGene = this.ChrGeneMap.get(chr);
                for (Gene gene: chrSelectGene) {
                    String geneId = gene.getGeneId();
                    selectGeneList.add(geneId);
                    TranscriptRecord longestTranscript = gene.getLongestTranscriptRecord();
                    String transcriptId = longestTranscript.getTranscriptId();
                    selectedTranscriptList.add(transcriptId);
                }
            }
            String line = "";
            String[] lineInfo;
            GTFReader gtfReader = new GTFReader();
            GeneAttribute geneAttribute;
            String geneId, transcriptId;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#")) {
                        bfw.write(line);
                        bfw.newLine();
                    } else {
                        lineInfo = line.split("\t");
                        geneAttribute = gtfReader.parseAttributes(lineInfo[8]);
                        geneId = geneAttribute.getGeneID();
                        transcriptId = geneAttribute.getTranscriptID();
                        if (selectGeneList.contains(geneId) && lineInfo[2].equals("gene")){
                            bfw.write(line);
                            bfw.newLine();
                        } else if (selectGeneList.contains(geneId) && selectedTranscriptList.contains(transcriptId)) {
                            bfw.write(line);
                            bfw.newLine();
                        }
                    }
                }
            }
            bfr.close();
            bfw.close();
        }catch (IOException ie) {
            ie.printStackTrace();
        } finally {
            if (bfw != null) {
                try {
                    bfw.close();
                } catch (IOException ie) {
                    ie.printStackTrace();
                }
            }
            if (bfr != null) {
                try {
                    bfr.close();
                } catch (IOException ie) {
                    ie.printStackTrace();
                }
            }
        }
    }

    /**
     * 为每个挑选出的基因设置RPKM值和PM参数，通过设置的参数计算表达值和reads数目，写入到文件
     * @param outputfile 输出文件路径
     */
    private void transcriptionParameter(String outputfile, String geneExpFile) {
        // 用于随机抽取RPKM值
        UniformRealDistribution unidiform = new UniformRealDistribution(100, 400);
        GeneExpDistribution geneExp = GeneExpDistribution.getInstance();
        HashMap<String, double[]> geneExpValue = null;
        Set<String> mutatedGeneId = this.geneMutatedPosition.keySet();
        // 如果有基因表达的参考数据，则读取数据
        if (geneExpFile != null)
            geneExpValue = geneExp.experimentalGeneExp(geneExpFile);
        try {
            FileWriter fw = new FileWriter(outputfile);
            fw.write("chr\tgeneId\tgeneName\tRPKM\texonLength\treadsCount\tmajorAlleleRatio\tminorAlleleRatio\n");
            for (Map.Entry<String, LinkedList<Gene>> entry : this.ChrGeneMap.entrySet()) {
                LinkedList<Gene> geneList = entry.getValue();
                String chr = entry.getKey();
                String geneId;
                double majorAlleleRatio, minorAlleleRatio;
                for (Gene gene: geneList) {
                    geneId = gene.getGeneId();
                    double RPKM;
                    // 如果没有设置测序深度
                    if (this.sequencingDepth == 0) {
                        if (geneExpValue == null)   // 没有参考的基因表达文件
                            RPKM = unidiform.sample();
                        else {      // 有参考的基因表达文件
                            double[] expData = geneExpValue.getOrDefault(gene.getGeneName(), null);
                            if (expData == null)    // 如果参考文件中没有记录该基因，则还是随机抽取RPKM值
                                RPKM = unidiform.sample();
                            else {                  // 如果参考文件中有记录该基因，则拟合正态分布，从分布中抽取RPKM值
                                NormalDistribution dist = new NormalDistribution(expData[0], expData[1] + 0.0001);
                                RPKM = Math.abs(dist.sample());
                                dist = null;
                            }
                        }
                        gene.setRPKM(RPKM); // 为基因设置RPKM值，在没有设置深度时通过文库大小和RPKM计算基因在测序时的reads数
                        gene.calculateReadsCountViaLibrarySize(this.librarySize);
                    } else { // 如果设定了测序深度，则直接通过深度求取reads数目和RPKM
                        gene.calculateReadsCountViaSequencingDepth(this.sequencingDepth, this.readLength, this.librarySize);
                        RPKM = gene.getRPKM();
                    }
                    if (mutatedGeneId.contains(geneId)) {
                        majorAlleleRatio = this.aseGeneMajorAlleleRatio.get(geneId);
                        minorAlleleRatio = 1.0 - majorAlleleRatio;
                    } else {
                        majorAlleleRatio = 1.0;
                        minorAlleleRatio = 0.0;
                    }
                    fw.write(chr + "\t" + geneId + "\t" + gene.getGeneName() + "\t" + RPKM + "\t"
                                  + gene.getExonSeq().length() + "\t" + gene.getReadsCount() + "\t"
                                  + majorAlleleRatio+ "\t" + minorAlleleRatio+"\n");
                }
                geneList = null;
            }
            fw.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        // 释放内存
        geneExpValue = null;
    }

    /**
     * 对每个选中的基因进行m6A修饰，在其外显子序列上随机生成一定数目的m6A修饰位点并将其记录到文件中，同时得到每个甲基化位点的ASM ratio
     */
    private void geneM6aModification() {
        this.geneM6aGenomePosition = new HashMap<>();

        M6AGenerator m6AGenerator = new M6AGenerator(this.geneMutatedPosition);
        HashMap<Integer, Integer> m6aSites;
        for (String chrNum: this.ChrGeneMap.keySet()) {
            LinkedList<Gene> chrGenes = this.ChrGeneMap.get(chrNum);
            for (Gene gene: chrGenes) {
                m6aSites = m6AGenerator.generateM6aSites(gene);
                this.geneM6aGenomePosition.put(gene.getGeneId(), m6aSites);
            }
        }
        File recordFile = new File(this.outputDir, "simulateM6aSites.txt");
        this.geneM6aAsmRatio = m6AGenerator.storeGeneM6aSites(this.geneM6aGenomePosition, recordFile);
    }

    /**
     * 将模拟的基因上SNP位点信息写入到文件
     * @param outputFile 输出文件路径
     */
    private void mutateGeneInfomation(String outputFile) {
        Set<String> geneIds = this.geneMutatedPosition.keySet();
        BufferedWriter bfw = null;
        try {
            bfw = new BufferedWriter(
                    new OutputStreamWriter(new FileOutputStream(new File(outputFile)))
            );
            bfw.write("chr\tgeneId\tgenomePos\texonSeqPos\tref\talt\n");
            String locateChr, writeOut;
            List<Integer> positions;
            for (String geneId: geneIds) {
                locateChr = this.getGeneChr(geneId);
                positions = new ArrayList<>(this.geneMutatedPosition.get(geneId));
                Collections.sort(positions);
                for (Integer position: positions) {
                    int genomePos = this.geneMutGenomePosition.get(geneId).get(position);
                    String[] refAndAlt = this.geneMutationRefAlt.get(geneId).get(position);
                    String ref = refAndAlt[0];
                    String alt = refAndAlt[1];
                    writeOut = locateChr + "\t" + geneId + "\t" + genomePos + "\t" + position + "\t" + ref.toUpperCase() + "\t" + alt.toUpperCase();
                    bfw.write(writeOut);
                    bfw.newLine();
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
    }

    /**
     * 获取某个基因所在的染色体号
     * @param geneId Gene 对象的geneId
     * @return 染色体号
     */
    private String getGeneChr(String geneId) {
        String locateChr = null;
        for (String chr: this.selectedGeneChr.keySet()) {
            HashSet<String> chrGenes = this.selectedGeneChr.get(chr);
            if (chrGenes.contains(geneId))
                locateChr = chr;
        }
        return locateChr;
    }

    /**
     * 选取的基因生成模拟测序数据的reads。将生成的reads写入文件
     * @param fragmentMean fragment平均长度
     * @param fragmentTheta fragment长度标准差
     * @param inputOutputFile1 input样本输出文件
     * @param inputOutputFile2 input样本输出文件(生成pair-end时传入，否则设置为null)
     * @param ipOutputFile1 ip样本输出文件
     * @param ipOutputFile2 ip样本输出文件(生成pair-end时传入，否则设置为null)
     * @param inputMultiple 生物学重复
     */
    private void generateReads(int fragmentMean, int fragmentTheta, String inputOutputFile1, String inputOutputFile2,
                               String ipOutputFile1, String ipOutputFile2, int inputMultiple) {
        BufferedWriter inputMate1File = null, inputMate2File = null, ipMate1File = null, ipMate2File = null;
        String direct = "SE";
        try {
            inputMate1File = new BufferedWriter(
                    new OutputStreamWriter(new FileOutputStream(new File(inputOutputFile1)))
            );
            if (inputOutputFile2 != null) {
                inputMate2File = new BufferedWriter(
                        new OutputStreamWriter(new FileOutputStream(new File(inputOutputFile2)))
                );
                direct = "PE";
            }
            ipMate1File = new BufferedWriter(
                    new OutputStreamWriter(new FileOutputStream(new File(ipOutputFile1)))
            );
            if (inputOutputFile2 != null) {
                ipMate2File = new BufferedWriter(
                        new OutputStreamWriter(new FileOutputStream(new File(ipOutputFile2)))
                );
                direct = "PE";
            }
            // 获取到突变基因的 geneID集合
            Set<String> mutGeneIds = this.geneMutatedPosition.keySet();
            for (Map.Entry<String, LinkedList<Gene>> entry : this.ChrGeneMap.entrySet()) {
                LinkedList<Gene> GeneList = entry.getValue();
                Set<Integer> geneM6aSites;  // 基因外显子区域的m6A修饰位点集合
                Set<Integer> geneMutateSites; // 基因外显子区域突变位点集合
                for (Gene gene: GeneList) {
                    geneMutateSites = null;
                    String geneId = gene.getGeneId();
                    // 如果当前基因在突变基因的列表中，则获取它突变后的外显子序列
                    String mutExonSeq = null;
                    double ref = 1.0;
                    if (mutGeneIds.contains(geneId)) {
                        ref = this.aseGeneMajorAlleleRatio.get(geneId);
                        mutExonSeq = this.genesMutatedExonSeq.get(geneId);
                        geneMutateSites = this.geneMutatedPosition.get(geneId);
                    }
                    geneM6aSites = this.geneM6aGenomePosition.get(geneId).keySet();
                    // 生成基因mRNA的fragment，并将其富集到INPUT和IP两个类群
                    gene.enrichInputFragment(fragmentMean, fragmentTheta, this.readLength, geneMutateSites);
                    gene.enrichIpFragment(fragmentMean, fragmentTheta, this.readLength, geneM6aSites, geneMutateSites);
                    gene.generateReads(inputMate1File, inputMate2File, this.readLength, fragmentMean, fragmentTheta, ref,
                                       mutExonSeq, direct, this.seqErrorModel, "input", null);
                    HashMap<Integer, Double> geneM6aAsm = this.geneM6aAsmRatio.getOrDefault(geneId, null);
                    gene.generateReads(ipMate1File, ipMate2File, this.readLength, fragmentMean, fragmentTheta, ref, mutExonSeq,
                                       direct, this.seqErrorModel, "ip", geneM6aAsm);
                    gene.peakFragmentFromBackground(ipMate1File, ipMate2File, this.readLength, fragmentMean, fragmentTheta, ref,
                                                    mutExonSeq, direct, this.seqErrorModel, geneMutateSites, geneM6aAsm);
                    gene.release();
                }
            }
            inputMate1File.close();
            ipMate1File.close();
            if (inputMate2File != null)
                inputMate1File.close();
            if (ipMate2File != null)
                ipMate2File.close();
        } catch (Exception io) {
            io.printStackTrace();
        } finally {
            if (inputMate1File != null) {
                try {
                    inputMate1File.close();
                } catch (IOException ie) {
                    ie.printStackTrace();
                }
            }
            if (inputMate2File != null) {
                try {
                    inputMate2File.close();
                } catch (IOException ie) {
                    ie.printStackTrace();
                }
            }
            if (ipMate2File != null) {
                try {
                    ipMate2File.close();
                } catch (IOException ie) {
                    ie.printStackTrace();
                }
            }
            if (ipMate1File != null) {
                try {
                    ipMate1File.close();
                } catch (IOException ie) {
                    ie.printStackTrace();
                }
            }
        }
    }


    /**
     * 生成模拟数据的方法，公有方法，供外界调用
     */
    public void simulateSequencing(String dataPath, String vcfFile, int librarySize, int sequencingDepth, int readLength,
                                   int minimumMut, int maximumMut, int fragmentMean, int fragmentTheta, double mutProp,
                                   int Multiple, int repeat, double pcrErrorProb, boolean overlap, boolean singleEnd,
                                   String geneExpFile) {

        this.outputDir = new File(dataPath).getAbsolutePath();
        this.seqErrorModel = new SequencingError(pcrErrorProb);
        this.setLibrarySize(librarySize);
        this.setSequencingDepth(sequencingDepth);
        this.setReadLength(readLength);
        // 每条染色体上选取基因用于生成数据
        this.selectGene(overlap);
        System.out.println("complete select");
        // 随机从选取的基因中抽取一部分进行突变，并设置ASE表达的major allele的比例
        this.randomMutateGene(vcfFile, minimumMut, maximumMut, mutProp);
        this.mutatedGeneAseRatio();
        System.out.println("complete variant");
        // 将随机生成的基因SNP信息写入文件
        this.mutateGeneInfomation(new File(this.outputDir, "mutations.txt").getAbsolutePath());
        // m6A修饰在mRNA上很常见，对每个选中的基因随机生成m6A修饰位点并写入文件
        this.geneM6aModification();
        System.out.println("complete m6A modification");
        // 生成简化的 GTF文件
        this.createSimpleGTF(new File(outputDir, "sim.chr.gtf").getAbsolutePath());
        // 计算基因的RPKM值并写入到文件
        this.transcriptionParameter(new File(outputDir, "rpkm.txt").getAbsolutePath(), geneExpFile);

        // 依据设定的实验重复次数生成IP和INPUT文件
        for (int i = 0; i < repeat; i++) {
            String InputOutputfile = new File(outputDir, "Input"+i+".fasta").getAbsolutePath();
            String IpOutputfile = new File(outputDir, "Ip"+i+".fasta").getAbsolutePath();
            String InputMatefile = null;
            String IpMatefile = null;
            // 单端测序
            if (singleEnd) {
                this.generateReads(fragmentMean,fragmentTheta, InputOutputfile, InputMatefile, IpOutputfile, IpMatefile, Multiple);
            } else { // 双端测序
                InputMatefile = new File(outputDir, "Input_mate"+i+".fasta").getAbsolutePath();
                IpMatefile = new File(outputDir, "Ip_mate"+i+".fasta").getAbsolutePath();
                this.generateReads(fragmentMean,fragmentTheta, InputOutputfile, InputMatefile, IpOutputfile, IpMatefile, Multiple);
            }
        }
        this.ChrGeneMap = null;
    }

    public HashMap<String, HashSet<Integer>> getGeneMutatedPosition() {
        return this.geneMutatedPosition;
    }

    public HashMap<String, String> getGenesMutatedExonSeq() {
        return this.genesMutatedExonSeq;
    }
}
