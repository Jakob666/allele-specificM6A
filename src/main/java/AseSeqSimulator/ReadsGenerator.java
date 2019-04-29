package AseSeqSimulator;

import GTFComponent.*;
import PeakSimulator.PeakSimulator;
import org.apache.commons.math3.distribution.UniformRealDistribution;

import java.io.*;
import java.util.*;

public class ReadsGenerator {
    private HashMap<String, LinkedList<Gene>> ChrGeneMap = new HashMap<String, LinkedList<Gene>>();
    private GTFReader readGTF = new GTFReader();
    private int readLength, peakLength;
    private double geneProp;
    private TwoBitParser twoBit;
    private long librarySize;
    private String gtfFile;
    private HashMap<String, String> genesOriginExonSeq, genesMutatedExonSeq;
    private HashMap<String, HashSet<Integer>> geneMutatedPosition;
    private HashMap<String, HashMap<Integer, String[]>> geneMutationRefAlt;
    private HashMap<String, HashMap<Integer, Integer>> geneMutGenomePosition;
    private HashMap<String, HashSet<String>> selectedGeneChr = new HashMap<>();
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

    public void setLibrarySize(long librarySize) {
        this.librarySize = librarySize;
    }

    public void setReadLength(int readLength) {
        this.readLength = readLength;
    }

    public void setPeakLength(int peakLength) {
        this.peakLength = peakLength;
    }

    public HashMap<String, LinkedList<Gene>> getChrGeneMap() {
        return this.ChrGeneMap;
    }

    public GTFReader getReadGTF() {
        return this.readGTF;
    }

    /**
     * select a group of genes from each chromosome, the group size equals to geneNum. These selected genes are stored
     * in chrGeneMap
     */
    public void selectGene(boolean overlapGene) {
        HashMap<String, ChromosomeRecord> chrMap = readGTF.getChromosomeMap();
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
                if (chr.equals("MT"))
                    continue;
                ChromosomeRecord chromosomeRecord = chrMap.get(chr);
                genesOnChr = chromosomeRecord.getChrGenes();

                // get all geneId and then shuffle for random select
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

                    // ignore the short exonic region gene
                    if (exonLength <= 1000)
                        continue;
                    transStart = longestTranscript.getTranscriptStart();
                    transEnd = longestTranscript.getTranscriptEnd();

                    if (overlapGene) { // there's no overlap regions between each two selected genes
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
                    Gene gene = new Gene(gr.getGeneId(), gr.getGeneStart(), gr.getGeneEnd(), gr.getStrand(), chr);
                    gene.setLongestTranscriptRecord(longestTranscript);
                    // get the longest transcript exon sequence with Gene splicing method
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
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    /**
     * generate a simplified GTF file of selected genes.
     * @param outputFile output simplified GTF file
     */
    public void createSimpleGTF(String outputFile) {
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
     * set RPKM and PM parameters for each selected gene
     * @param outputfile output file path
     */
    public void transcriptionParameter(String outputfile) {
        double chr_pm_range = 0.95 / ChrGeneMap.keySet().size();
        double lower;
        double upper = 0;
        UniformRealDistribution unidiform = new UniformRealDistribution(10, 1000);
        try {
            FileWriter fw = new FileWriter(outputfile);
            for (Map.Entry<String, LinkedList<Gene>> entry : ChrGeneMap.entrySet()) {
                LinkedList<Gene> GeneList = entry.getValue();
                String chr = entry.getKey();
                double gene_pm_range = chr_pm_range / GeneList.size();
                for (Gene gene: GeneList) {
                    lower = upper;
                    upper = upper + gene_pm_range;
                    gene.setPmRange(lower, upper);
                    double RPKM = unidiform.sample();
                    gene.setRPKM(RPKM);
                    fw.write(chr + "\t" + gene.getGeneId() + "\t" + RPKM + "\n");
                }
            }
            fw.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    /**
     * randomly generate m6A peak for selected genes
     */
    public void simulateM6aPeaks() {
        Set<String> mutatedGenesId = this.geneMutatedPosition.keySet();
        LinkedList peaks;
        for (String chrNum: this.ChrGeneMap.keySet()) {
            LinkedList<Gene> genes = this.ChrGeneMap.get(chrNum);
            for (Gene gene: genes) {
                String geneId = gene.getGeneId();
                if (mutatedGenesId.contains(geneId)) {
                    HashSet<Integer> mutPositions = this.geneMutatedPosition.get(geneId);
                    peaks = PeakSimulator.altGenePeakSimulation(this.peakLength, gene, mutPositions);
                    gene.setPeakList(peaks);
                } else {
                    peaks = PeakSimulator.refGenePeakSimulation(this.peakLength, gene);
                    gene.setPeakList(peaks);
                }
            }
        }
    }

    /**
     * write the simulated mutation informations into file system
     */
    public void mutateGeneInfomation(String outputFile) {
        Set<String> geneIds = this.geneMutatedPosition.keySet();
        BufferedWriter bfw = null;
        try {
            bfw = new BufferedWriter(
                    new OutputStreamWriter(new FileOutputStream(new File(outputFile)))
            );
            bfw.write("chr\tgeneId\tgenomePos\texonSeqPos\tref\talt\n");
            int refCount, altCount;
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
                    writeOut = locateChr + "\t" + geneId + "\t" + genomePos + "\t" + position + "\t" + ref + "\t" + alt;
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
     * generate INPUT sequencing reads for each selected gene. The simulated sequencing reads will be write out into the
     * output file
     * @param fragmentMean fragment mean
     * @param fragmentTheta fragment std
     * @param inputOutputFile output file path
     * @param inputMultiple replicates for each experiment
     */
    public void generateInputReads(int fragmentMean, int fragmentTheta, String inputOutputFile, int inputMultiple) {
        BufferedWriter fw = null;
        UniformRealDistribution urd = new UniformRealDistribution(0.0, 0.95);
        try {
            fw = new BufferedWriter(
                    new OutputStreamWriter(new FileOutputStream(new File(inputOutputFile)))
            );
            // get mutated genes' ID
            Set<String> mutGeneIds = this.geneMutatedPosition.keySet();
            for (Map.Entry<String, LinkedList<Gene>> entry : this.ChrGeneMap.entrySet()) {
                LinkedList<Gene> GeneList = entry.getValue();
                for (Gene gene: GeneList) {
                    String geneId = gene.getGeneId();
                    // if gene in mutate gene list, get its mutated exon sequence
                    String mutExonSeq = null;
                    double ref = 1.0;
                    if (mutGeneIds.contains(geneId)) {
                        ref = urd.sample();
                        mutExonSeq = this.genesMutatedExonSeq.get(geneId);
                    }
                    gene.calculateReadsCount(librarySize);

                    gene.generateInputReads(fragmentMean, fragmentTheta, fw, this.readLength, inputMultiple, ref, mutExonSeq,
                                       this.geneMutatedPosition.getOrDefault(geneId, null), this.seqErrorModel);
//                    HashMap<Integer, int[]> refAltReads = gene.getReadsCountRecord();
//                    this.geneMutReadsCount.put(geneId, refAltReads);
                }
            }
            fw.close();
        } catch (Exception io) {
            io.printStackTrace();
        } finally {
            if (fw != null) {
                try {
                    fw.close();
                } catch (IOException ie) {
                    ie.printStackTrace();
                }
            }
        }
    }

    /**
     * generate IP sequencing reads for each selected gene. The simulated sequencing reads will be write out into the
     * output file
     * @param fragmentMean fragment mean
     * @param fragmentTheta fragment std
     * @param ipOutputFile output file path
     * @param ipMultiple replicates for each experiment
     */
    public void generateIpReads(int fragmentMean, int fragmentTheta, String ipOutputFile, int ipMultiple) {
        BufferedWriter fw = null;
        UniformRealDistribution urd = new UniformRealDistribution(0.0, 0.95);
        try {
            fw = new BufferedWriter(
                    new OutputStreamWriter(new FileOutputStream(new File(ipOutputFile)))
            );
            // get mutated genes' ID
            Set<String> mutGeneIds = this.geneMutatedPosition.keySet();

            for (Map.Entry<String, LinkedList<Gene>> entry : ChrGeneMap.entrySet()) {
                LinkedList<Gene> GeneList = entry.getValue();
                for (Gene gene: GeneList) {
                    String geneId = gene.getGeneId();
                    String mutExonSeq = null;
                    double ref = 1.0;
                    if (mutGeneIds.contains(geneId)) {
                        ref = urd.sample();
                        mutExonSeq = this.genesMutatedExonSeq.get(geneId);
                    }

                    gene.calculateReadsCount(librarySize);

                    gene.generateIpReads(fragmentMean, fragmentTheta, fw, this.readLength, ipMultiple, ref, mutExonSeq,
                            this.geneMutatedPosition.getOrDefault(geneId, null), this.seqErrorModel);

//                    if (mutGeneIds.contains(geneId)) {
//                        int refCount = gene.getRefReadsCount();
//                        int altCount = gene.getAltReadsCount();
//                        int[] refAndAlt = new int[]{refCount, altCount};
//                        this.mutGeneRefAltCount.put(geneId, refAndAlt);
//                    }
                }
            }
            fw.close();
        } catch (Exception io) {
            io.printStackTrace();
        } finally {
            if (fw != null) {
                try {
                    fw.close();
                } catch (IOException ie) {
                    ie.printStackTrace();
                }
            }
        }
    }

    /**
     * sequencing reads simulating
     */
    public void simulateSequencing(String dataPath, String vcfFile, int librarySize, int peakLength, int readLength,
                                   int minimumMut, int maximumMut, int fragmentMean, int fragmentTheta, double mutProp,
                                   int Multiple, int repeat, double pcrErrorProb, boolean overlap) {

        this.seqErrorModel = new SequencingError(pcrErrorProb);
        this.setLibrarySize(librarySize);
        this.setPeakLength(peakLength);
        this.setReadLength(readLength);
        this.selectGene(overlap);
        System.out.println("complete select");

        // randomly generate mutations on selected genes' exon sequence
        SNPGenerator snpGenerator = new SNPGenerator(this.ChrGeneMap, mutProp, vcfFile, minimumMut, maximumMut);
        System.out.println("complete variant");
        this.genesOriginExonSeq = snpGenerator.getOriginExonSequence();
        this.genesMutatedExonSeq = snpGenerator.getMutatedExonSeqence();
        this.geneMutatedPosition = snpGenerator.getMutGenePosition();
        this.geneMutationRefAlt = snpGenerator.getMutationRefAlt();
        this.geneMutGenomePosition = snpGenerator.getMutGenomePosition();
        String outputDir = new File(dataPath).getAbsolutePath();

        // randomly generate m6A peak, write out the simulated peak record
//        this.simulateM6aPeaks();
//        File simulatedPeakFile = new File(outputDir, "simulatedPeak.txt");
//        PeakSimulator.storeSimulatedPeaksRecord(this.ChrGeneMap, simulatedPeakFile);
//        File simulatedPeakBedFile = new File(outputDir, "simulatedPeak.bed");
//        PeakSimulator.Trans2bed(simulatedPeakFile, simulatedPeakBedFile);
//        File simulatedPeakCenterFile = new File(outputDir, "simulatedPeakCenter.bed");
//        PeakSimulator.Trans2bed_peaksite(simulatedPeakFile, simulatedPeakCenterFile, readLength / 2);

        // generate a simplified GTF file
        this.createSimpleGTF(new File(outputDir, "sim.chr.gtf").getAbsolutePath());
        // record the RPKM data
        this.transcriptionParameter(new File(outputDir, "rpkm.txt").getAbsolutePath());

        for (int i = 0; i < repeat; i++) {
            String InputOutputfile = new File(outputDir, "Input"+i+".fasta").getAbsolutePath();
            String IpOutputfile = new File(outputDir, "Ip"+i+".fasta").getAbsolutePath();
            this.generateInputReads(fragmentMean,fragmentTheta, InputOutputfile,Multiple);
            this.generateIpReads(fragmentMean,fragmentTheta, IpOutputfile,Multiple);
            // write out gene mutation records of each experiment
            mutateGeneInfomation(new File(outputDir, "mutations"+i+".txt").getAbsolutePath());

        }

    }

    public HashMap<String, HashSet<Integer>> getGeneMutatedPosition() {
        return this.geneMutatedPosition;
    }

    public HashMap<String, String> getGenesMutatedExonSeq() {
        return this.genesMutatedExonSeq;
    }
}
