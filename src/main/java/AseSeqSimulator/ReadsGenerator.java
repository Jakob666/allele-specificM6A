package AseSeqSimulator;

import GTFComponent.*;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;

import java.io.*;
import java.util.*;

public class ReadsGenerator {
    private HashMap<String, LinkedList<Gene>> ChrGeneMap = new HashMap<String, LinkedList<Gene>>();
    private GTFReader readGTF = new GTFReader();
    private int readLength, sequencingDepth, m6aPeakLength;
    private double geneProp, aseInfimum, aseSupremum;
    private TwoBitParser twoBit;
    private long librarySize;
    private String gtfFile, outputDir;
    private HashMap<String, String> genesOriginExonSeq, genesMutatedExonSeq;
    private HashMap<String, HashSet<Integer>> geneMutatedPosition;
    private HashMap<String, HashMap<Integer, String[]>> geneMutationRefAlt;
    private HashMap<String, HashMap<Integer, Integer>> geneMutGenomePosition;
    private HashMap<String, HashMap<String, HashMap<Integer, Integer>>> geneM6aGenomePosition;
    private HashMap<String, HashMap<Integer, Double>> geneM6aAsmRatio;
    private HashMap<String, HashMap<Integer, Boolean>> geneM6aAsmBias;
    private HashMap<String, HashSet<String>> selectedGeneChr = new HashMap<>();
    private HashMap<String, Double> aseGeneMajorAlleleRatio = new HashMap<>();
    private HashMap<String, ArrayList<String>> geneM6aPeakRanges = new HashMap<>();
    private HashMap<String, ArrayList<Boolean>> geneM6aPeakAsm = new HashMap<>();
    private SequencingError seqErrorModel;
    private HashMap<String, HashMap<Integer, int[]>> geneRnaMutationCoverage = new HashMap<>(), geneDnaMutationCoverage = new HashMap<>();

    public ReadsGenerator(String gtfFile, double geneProp, double aseInfimum, double aseSupremum, String twoBitFile) {
        this.readGTF.readFromFile(gtfFile);
        System.out.println("read GTF complete");
        this.geneProp = geneProp;
        this.aseInfimum = aseInfimum;
        this.aseSupremum = aseSupremum;
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

    private void setM6aPeakLength(int peakLength) {
        this.m6aPeakLength = peakLength;
    }

    public HashMap<String, LinkedList<Gene>> getChrGeneMap() {
        return this.ChrGeneMap;
    }

    /**
     * select genes from each chromosome
     * @param overlapGene if overlap between genes
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
                // 2bit file contains no chrMT
                if (chr.equals("MT") || chr.equals("X") || chr.equals("Y"))
                    continue;
                ChromosomeRecord chromosomeRecord = chrMap.get(chr);
                genesOnChr = chromosomeRecord.getChrGenes();

                geneIds = new ArrayList<>(genesOnChr.keySet());
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

                    // skip gene with short exonic regions
                    if (exonLength <= 1000)
                        continue;
                    transStart = longestTranscript.getTranscriptStart();
                    transEnd = longestTranscript.getTranscriptEnd();

                    if (!overlapGene) {
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
                    // get longest transcriptome
                    this.twoBit.setCurrentSequence("chr"+chr);
                    gene.splicing(this.twoBit);
                    geneList.add(gene);
                    twoBit.close();
                    chrGenes.add(gene.getGeneId());
                    order ++;
                }
                this.selectedGeneChr.put(chr, chrGenes);
                this.ChrGeneMap.put(chr, geneList);
            }
            // release
            transcriptRegion = null;
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        // release
        chrMap.clear();
        this.readGTF = null;
    }

    /**
     * generate SNP site
     * @param vcfFile form SNV sites via existed VCF file, default null
     * @param minimumMut minimum mutation number on exon
     * @param maximumMut maximum mutation number on exon
     * @param mutProp proportion of mutation gene in total gene
     */
    private void randomMutateGene(String vcfFile, int minimumMut, int maximumMut, double mutProp) {
        SNPGenerator snpGenerator = new SNPGenerator(this.ChrGeneMap, this.geneM6aPeakRanges, mutProp, vcfFile,
                                                     minimumMut, maximumMut);
        snpGenerator.generateSNP();
        this.genesOriginExonSeq = snpGenerator.getOriginExonSequence();     // gene exon sequence without SNV site
        this.genesMutatedExonSeq = snpGenerator.getMutatedExonSeqence();    // gene mutated exon sequence
        this.geneMutatedPosition = snpGenerator.getMutGenePosition();       // gene mutation position on exonic region
        this.geneMutationRefAlt = snpGenerator.getMutationRefAlt();
        this.geneMutGenomePosition = snpGenerator.getMutGenomePosition();   // gene mutation site genome position
        snpGenerator = null;
    }

    /**
     * mutated gene major allele frequency
     */
    private void mutatedGeneAseRatio() {
        UniformRealDistribution urd = null;
        if (!(Math.abs(this.aseInfimum-this.aseSupremum) < 0.00001))
            urd = new UniformRealDistribution(this.aseInfimum, this.aseSupremum);
        // mutation geneID
        Set<String> mutGeneIds = this.geneMutatedPosition.keySet();
        double ref, randNum;
        for (String mutGeneId: mutGeneIds) {
            randNum = Math.random();
            if (randNum > 0.5)
                this.aseGeneMajorAlleleRatio.put(mutGeneId, 0.5);
            else {
                ref = (urd != null)? urd.sample(): this.aseInfimum;
                this.aseGeneMajorAlleleRatio.put(mutGeneId, ref);
            }
        }
    }

    /**
     * generate simplified GTF file
     * @param outputFile GTF file path
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
     * set gene RPKM value
     * @param rnaRPKMFile rna RPKM output file path
     * @param dnaRPKMFile dna depth output file path
     * @param geneExpFile gene expression file, default null
     * @param minReadsCount minimum reads coverage
     * @param maxReadsCount maximum reads coverage
     */
    private void transcriptionParameter(String rnaRPKMFile, String dnaRPKMFile, String geneExpFile,
                                        int minReadsCount, int maxReadsCount) {
        UniformRealDistribution unidiform = new UniformRealDistribution(minReadsCount, maxReadsCount);
        GeneExpDistribution geneExp = GeneExpDistribution.getInstance();
        HashMap<String, double[]> geneExpValue = null;
        Set<String> mutatedGeneId = this.geneMutatedPosition.keySet();
        if (geneExpFile != null)
            geneExpValue = geneExp.experimentalGeneExp(geneExpFile);

        BufferedWriter rfw = null, dfw = null;
        try {
            rfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(rnaRPKMFile))));
            dfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(dnaRPKMFile))));
            rfw.write("chr\tgeneId\tgeneName\tRPKM\texonLength\treadsCount\tmajorAlleleRatio\tminorAlleleRatio\n");
            dfw.write("chr\tgeneId\tgeneName\tRPKM\texonLength\treadsCount\tmajorAlleleRatio\tminorAlleleRatio\n");
            for (Map.Entry<String, LinkedList<Gene>> entry : this.ChrGeneMap.entrySet()) {
                LinkedList<Gene> geneList = entry.getValue();
                String chr = entry.getKey();
                String geneId;
                double majorAlleleRatio, minorAlleleRatio;
                for (Gene gene: geneList) {
                    geneId = gene.getGeneId();
                    double rnaRPKM, dnaRPKM = 0;
                    // RNA-seq表达值
                    if (geneExpValue == null)
                        rnaRPKM = unidiform.sample();
                    else {      // 有参考的基因表达文件
                        double[] expData = geneExpValue.getOrDefault(gene.getGeneName(), null);
                        if (expData == null)
                            rnaRPKM = unidiform.sample();
                        else {
                            NormalDistribution dist = new NormalDistribution(expData[0], expData[1] + 0.0001);
                            rnaRPKM = Math.abs(dist.sample());
                            dist = null;
                        }
                    }
                    gene.setRnaRPKM(rnaRPKM);
                    gene.calculateReadsCountViaLibrarySize(this.librarySize);
                    // DNA sequencing depth
                    if (this.sequencingDepth != 0) {
                        gene.calculateReadsCountViaSequencingDepth(this.sequencingDepth, this.readLength);
                        dnaRPKM = gene.getDnaRPKM();
                        gene.setDnaRPKM(dnaRPKM);
                    }

                    if (mutatedGeneId.contains(geneId)) {
                        majorAlleleRatio = this.aseGeneMajorAlleleRatio.get(geneId);
                        minorAlleleRatio = 1.0 - majorAlleleRatio;
                    } else {
                        majorAlleleRatio = 1.0;
                        minorAlleleRatio = 0.0;
                    }
                    rfw.write(chr + "\t" + geneId + "\t" + gene.getGeneName() + "\t" + rnaRPKM + "\t"
                            + gene.getExonSeq().length() + "\t" + gene.getRnaReadsCount() + "\t"
                            + majorAlleleRatio+ "\t" + minorAlleleRatio+"\n");

                    dfw.write(chr + "\t" + geneId + "\t" + gene.getGeneName() + "\t" + dnaRPKM + "\t"
                            + gene.getExonSeq().length() + "\t" + gene.getDnaReadsCount() + "\t"
                            + majorAlleleRatio+ "\t" + minorAlleleRatio+"\n");

                }
                geneList = null;
            }
        } catch (IOException ex) {
            ex.printStackTrace();
        } finally {
            if (rfw != null) {
                try {
                    rfw.close();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
            if (dfw != null) {
                try {
                    dfw.close();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
        // release
        geneExpValue = null;
    }

    /**
     * generate m6A modification position
     */
    private void geneM6aPeakRange(int maxPeakNum) {
        M6AGenerator m6AGenerator = new M6AGenerator();
        M6APeaks mp = M6APeaks.getInstance(m6AGenerator);
        for (String chrNum: this.ChrGeneMap.keySet()) {
            List<Gene> chrGenes = this.ChrGeneMap.get(chrNum);
            for (Gene gene: chrGenes) {
                mp.geneM6aPeakRange(gene, this.m6aPeakLength, this.readLength, maxPeakNum);
            }
        }
        this.geneM6aGenomePosition = mp.getGeneM6aSites();
        for (String label: this.geneM6aGenomePosition.keySet()) {
            String[] info = label.split(":");
            String geneId = info[1];
            HashMap<String, HashMap<Integer, Integer>> genePeakRecords = this.geneM6aGenomePosition.get(label);
            ArrayList<String> genePeaks = new ArrayList<>(genePeakRecords.keySet().size());
            ArrayList<Boolean> genePeakAse = new ArrayList<>(genePeakRecords.keySet().size());
            for (String peakRange: genePeakRecords.keySet()) {
                boolean ase;
                genePeaks.add(peakRange);
                double random = Math.random();
                ase = random > 0.5;
                genePeakAse.add(ase);
            }
            this.geneM6aPeakRanges.put(geneId, genePeaks);
            this.geneM6aPeakAsm.put(geneId, genePeakAse);
        }

        File recordFile = new File(this.outputDir, "simulateM6aSites.txt");
        m6AGenerator.storeGeneM6aSites(this.geneM6aGenomePosition, recordFile, this.geneM6aPeakRanges, this.geneM6aPeakAsm);
        this.geneM6aAsmRatio = m6AGenerator.getAseRatio();
        this.geneM6aAsmBias = m6AGenerator.getAseBias();
    }

    /**
     * record simulated m6A signal peak into BED format file
     * @param outputBedFile BED format file path
     */
    private void m6aPeakBedFile(String outputBedFile) {
        BufferedWriter bfw = null;
        try {
            bfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(outputBedFile))));
            bfw.write("# chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tASM\n");
            for (String label: this.geneM6aGenomePosition.keySet()) {
                String[] info = label.split(":");
                String chrNum = info[0], geneId = info[1], strand = info[2];
                ArrayList<String> genePeaks = this.geneM6aPeakRanges.get(geneId);
                ArrayList<Boolean> peakAsm = this.geneM6aPeakAsm.get(geneId);
                assert peakAsm.size() == genePeaks.size();
                for (int i =0; i < peakAsm.size(); i++) {
                    String peakRange = genePeaks.get(i);
                    boolean asm = peakAsm.get(i);
                    String[] range = peakRange.split(":");
                    String genomeStart = range[0], genomeEnd = range[1];
                    String[] lineInfo = new String[] {chrNum, genomeStart, genomeEnd, geneId, "0.001", strand,
                                        genomeStart, genomeEnd, "0", "1", Integer.toString(this.m6aPeakLength), "0",
                                        Boolean.toString(asm)};
                    String line = String.join("\t", lineInfo);
                    bfw.write(line);
                    bfw.newLine();
                }
            }
        } catch (IOException ie) {
            ie.printStackTrace();
        } finally {
            if (bfw != null) {
                try {
                    bfw.close();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
    }

    /**
     * record simulated SNV information
     * @param outputFile output file path
     */
    private void mutateGeneInfomation(String outputFile) {
        Set<String> geneIds = this.geneMutatedPosition.keySet();
        BufferedWriter bfw = null;
        try {
            bfw = new BufferedWriter(
                    new OutputStreamWriter(new FileOutputStream(new File(outputFile)))
            );
            bfw.write("#chr\tgeneId\tgenomePos\texonSeqPos\tref\talt\trnaMajorCoverage\trnaMinorCoverage\tdnaMajorCoverage\tdnaMinorCoverage\n");
            String locateChr, writeOut;
            List<Integer> positions;
            for (String geneId: geneIds) {
                locateChr = this.getGeneChr(geneId);
                HashMap<Integer, int[]> rnaMutationCoverage = this.geneRnaMutationCoverage.get(geneId),
                                        dnaMutationCoverage = this.geneDnaMutationCoverage.get(geneId);

                positions = new ArrayList<>(this.geneMutatedPosition.get(geneId));
                Collections.sort(positions);
                for (Integer position: positions) {
                    int genomePos = this.geneMutGenomePosition.get(geneId).get(position);
                    String[] refAndAlt = this.geneMutationRefAlt.get(geneId).get(position);
                    String ref = refAndAlt[0];
                    String alt = refAndAlt[1];
                    // sites covered by 0 reads caused by low expression
                    String rnaMajorCoverage = Integer.toString(rnaMutationCoverage.getOrDefault(position, new int[]{0, 0})[0]),
                           rnaMinorCoverage = Integer.toString(rnaMutationCoverage.getOrDefault(position, new int[]{0, 0})[1]),
                           dnaMajorCoverage = Integer.toString(dnaMutationCoverage.getOrDefault(position, new int[]{0, 0})[0]),
                           dnaMinorCoverage = Integer.toString(dnaMutationCoverage.getOrDefault(position, new int[]{0, 0})[1]);
                    writeOut = String.join("\t", new String[] {locateChr, geneId, Integer.toString(genomePos),
                            Integer.toString(position), ref.toUpperCase(), alt.toUpperCase(), rnaMajorCoverage,
                            rnaMinorCoverage, dnaMajorCoverage, dnaMinorCoverage});
                    bfw.write(writeOut);
                    bfw.newLine();
                }
            }
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
     * get gene chromosome
     * @param geneId geneId
     * @return chromosome
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
     * simulate sequencing reads
     * @param fragmentMean average fragment length
     * @param fragmentTheta std of fragment length
     * @param inputOutputFile1 INPUT sample file
     * @param inputOutputFile2 INPUT sample file(only use when generate pair-end data; otherwise, null)
     * @param ipOutputFile1 IP sample file
     * @param ipOutputFile2 IP sample file(only use when generate pair-end data; otherwise, null)
     * @param type dna or rna
     */
    private void generateReads(int fragmentMean, int fragmentTheta, String inputOutputFile1, String inputOutputFile2,
                               String ipOutputFile1, String ipOutputFile2, String type) {
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

            // mutation geneIDs
            Set<String> mutGeneIds = this.geneMutatedPosition.keySet();
            for (Map.Entry<String, LinkedList<Gene>> entry : this.ChrGeneMap.entrySet()) {
                LinkedList<Gene> GeneList = entry.getValue();
                Set<Integer> geneM6aSites = new HashSet<>();  // m6A sites on exon
                Set<Integer> geneMutateSites; // mutations on exon
                for (Gene gene: GeneList) {
                    geneMutateSites = null;
                    String geneId = gene.getGeneId(), chrNum = gene.getChr(), strand = gene.getStrand();
                    String label = String.join(":", new String[]{chrNum, geneId, strand});
                    String mutExonSeq = null;
                    double ref = 1.0;
                    if (mutGeneIds.contains(geneId)) {
                        if (type.equals("rna")) // generate RNA-seq data
                            ref = this.aseGeneMajorAlleleRatio.get(geneId);
                        else // do not need to take expression value when generate DNA sequencing data
                            ref = 0.5;
                        mutExonSeq = this.genesMutatedExonSeq.get(geneId);
                        geneMutateSites = this.geneMutatedPosition.get(geneId);
                    }
                    HashMap<String, HashMap<Integer, Integer>> peakCoveredSites = this.geneM6aGenomePosition.get(label);
                    for (String peak: peakCoveredSites.keySet()) {
                        HashMap<Integer, Integer> m6aSites = peakCoveredSites.get(peak);
                        for (Integer exonPosition: m6aSites.keySet())
                            geneM6aSites.add(exonPosition);
                    }
                    // INPUT and IP fragments
                    gene.enrichInputFragment(fragmentMean, fragmentTheta, this.readLength, geneMutateSites, type);
                    gene.generateReads(inputMate1File, inputMate2File, this.readLength, fragmentMean, fragmentTheta, ref,
                                       mutExonSeq, direct, this.seqErrorModel, "input", type, null, null);
                    if (type.equals("rna"))
                        this.geneRnaMutationCoverage.put(geneId, gene.getRnaMutationCoverage());
                    else
                        this.geneDnaMutationCoverage.put(geneId, gene.getDnaMutationCoverage());
                    HashMap<Integer, Double> geneM6aAsm = this.geneM6aAsmRatio.getOrDefault(geneId, null);
                    HashMap<Integer, Boolean> m6aBias = this.geneM6aAsmBias.getOrDefault(geneId, null);
                    gene.enrichIpFragment(fragmentMean, fragmentTheta, this.readLength, geneM6aSites, geneMutateSites, type);
                    gene.generateReads(ipMate1File, ipMate2File, this.readLength, fragmentMean, fragmentTheta, ref, mutExonSeq,
                                       direct, this.seqErrorModel, "ip", type, geneM6aAsm, m6aBias);
                    gene.peakFragmentFromBackground(ipMate1File, ipMate2File, this.readLength, fragmentMean, fragmentTheta, ref,
                                                    mutExonSeq, direct, this.seqErrorModel, geneMutateSites, geneM6aAsm, m6aBias);
                    gene.release();
                }
            }
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
     * generate simulated data
     */
    public void simulateSequencing(String dataPath, String vcfFile, int librarySize, int sequencingDepth, int readLength,
                                   int minimumMut, int maximumMut, int fragmentMean, int fragmentTheta, int minReadsCount,
                                   int maxReadsCount, double mutProp, int peakLength, int maxPeakNum, double pcrErrorProb,
                                   boolean overlap, boolean singleEnd, String geneExpFile) {

        this.outputDir = new File(dataPath).getAbsolutePath();
        this.seqErrorModel = new SequencingError(pcrErrorProb);
        this.setLibrarySize(librarySize);
        this.setSequencingDepth(sequencingDepth);
        this.setReadLength(readLength);
        this.setM6aPeakLength(peakLength);
        this.selectGene(overlap);
        System.out.println("complete select");
        // simulate m6A peaks
        this.geneM6aPeakRange(maxPeakNum);
        String outputBedFile = new File(this.outputDir, "simulatePeak.bed").getAbsolutePath();
        this.m6aPeakBedFile(outputBedFile);
        System.out.println("complete m6A modification");
        // simulate mutated genes
        this.randomMutateGene(vcfFile, minimumMut, maximumMut, mutProp);
        this.mutatedGeneAseRatio();
        System.out.println("complete variant");

        // generate simplified GTF file
        this.createSimpleGTF(new File(outputDir, "sim.chr.gtf").getAbsolutePath());
        // set genes' RPKM and record into file
        String rnaseqRPKMFile = new File(outputDir, "rnaseq_rpkm.txt").getAbsolutePath();
        String dnaseqRPKMFile = new File(outputDir, "dnaseq_rpkm.txt").getAbsolutePath();
        this.transcriptionParameter(rnaseqRPKMFile, dnaseqRPKMFile, geneExpFile, minReadsCount, maxReadsCount);

        // generate MeRIP-seq、MeDIP-seq IP and INPUT file
        File meRipOutputDir = new File(outputDir, "merip_seq");
        meRipOutputDir.mkdir();
        File meDipOutputDir = new File(outputDir, "medip_seq");
        meDipOutputDir.mkdir();
        // generate MeRIP-seq data
        String InputOutputFile = new File(meRipOutputDir, "Input.fasta").getAbsolutePath();
        String IpOutputFile = new File(meRipOutputDir, "Ip.fasta").getAbsolutePath();
        String InputMateFile, IpMateFile;
        if (singleEnd) {  // single-end
            this.generateReads(fragmentMean,fragmentTheta, InputOutputFile, null, IpOutputFile, null, "rna");
        } else { // pair-end
            InputMateFile = new File(meRipOutputDir, "Input_mate.fasta").getAbsolutePath();
            IpMateFile = new File(meRipOutputDir, "Ip_mate.fasta").getAbsolutePath();
            this.generateReads(fragmentMean,fragmentTheta, InputOutputFile, InputMateFile, IpOutputFile, IpMateFile, "rna");
        }
        // generate MeDIP-seq data
        InputOutputFile = new File(meDipOutputDir, "Input.fasta").getAbsolutePath();
        IpOutputFile = new File(meDipOutputDir, "Ip.fasta").getAbsolutePath();
        if (singleEnd) {  // single-end
            this.generateReads(fragmentMean,fragmentTheta, InputOutputFile, null, IpOutputFile, null, "dna");
        } else { // pair-end
            InputMateFile = new File(meDipOutputDir, "Input_mate.fasta").getAbsolutePath();
            IpMateFile = new File(meDipOutputDir, "Ip_mate.fasta").getAbsolutePath();
            this.generateReads(fragmentMean,fragmentTheta, InputOutputFile, InputMateFile, IpOutputFile, IpMateFile, "dna");
        }

        this.mutateGeneInfomation(new File(this.outputDir, "mutations.txt").getAbsolutePath());

        this.geneM6aGenomePosition.clear();
        this.geneM6aPeakAsm.clear();
        this.geneM6aPeakRanges.clear();
        this.geneMutatedPosition.clear();
        this.geneMutGenomePosition.clear();
        this.geneM6aAsmBias.clear();
        this.ChrGeneMap.clear();
        this.geneRnaMutationCoverage.clear();
        this.geneDnaMutationCoverage.clear();
    }

    public HashMap<String, HashSet<Integer>> getGeneMutatedPosition() {
        return this.geneMutatedPosition;
    }

    public HashMap<String, String> getGenesMutatedExonSeq() {
        return this.genesMutatedExonSeq;
    }
}
