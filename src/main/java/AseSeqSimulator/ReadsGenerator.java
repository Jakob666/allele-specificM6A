package AseSeqSimulator;

import GTFComponent.*;
import org.apache.commons.math3.distribution.UniformRealDistribution;

import java.io.*;
import java.util.*;

public class ReadsGenerator {
    private HashMap<String, LinkedList<Gene>> ChrGeneMap = new HashMap<String, LinkedList<Gene>>();
    private GTFReader readGTF = new GTFReader();
    private int geneNum, readLength;
    private TwoBitParser twoBit;
    private long librarySize;
    private String gtfFile;
    private HashMap<String, String> genesOriginExonSeq, genesMutatedExonSeq;
    private HashMap<String, HashSet<Integer>> geneMutatedPosition;
    private HashMap<String, int[]> mutGeneRefAltCount = new HashMap<>();
    private SequencingError seqErrorModel;

    public ReadsGenerator(String gtfFile, int geneNum, String twoBitFile) {
        this.readGTF.readFromFile(gtfFile);
        this.geneNum = geneNum;
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
    public void selectGene() {
        HashMap<String, ChromosomeRecord> chrMap = readGTF.getChromosomeMap();
        try {
            String chr, randomGeneId;
            HashMap<String, GeneRecord> genesOnChr;
            ArrayList<String> geneIds;
            LinkedList<Gene> geneList;
            LinkedList<int[]> transcriptRegion;

            for (Map.Entry<String, ChromosomeRecord> entry : chrMap.entrySet()) {
                transcriptRegion = new LinkedList<int[]>();
                chr = entry.getKey();
                ChromosomeRecord chromosomeRecord = entry.getValue();
                genesOnChr = chromosomeRecord.getChrGenes();

                // get all geneId and then shuffle for random select
                geneIds = new ArrayList<String>(genesOnChr.keySet());
                Collections.shuffle(geneIds);
                geneList = new LinkedList<>();
                int order = 0, exonLength = 0, transStart, transEnd;
                boolean overlap = false;
                while (geneList.size() < this.geneNum && order < geneIds.size()) {
                    randomGeneId = geneIds.get(order);
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
                    // there's no overlap regions between each two selected genes
                    for (Iterator<int[]> region_iter = transcriptRegion.iterator(); region_iter.hasNext(); ) {
                        int[] region = region_iter.next();
                        int start_max = Math.max(region[0], transStart);
                        int end_min = Math.min(region[1], transEnd);
                        if (start_max <= end_min) {
                            overlap = true;
                            break;
                        }
                    }
                    if (overlap)
                        continue;
                    transcriptRegion.add(new int[]{transStart, transEnd});
                    Gene gene = new Gene(gr.getGeneId(), gr.getGeneStart(), gr.getGeneEnd(), gr.getStrand(), chr);
                    gene.setLongestTranscriptRecord(longestTranscript);
                    this.twoBit.setCurrentSequence("chr"+chr);
                    gene.splicing(this.twoBit);
                    geneList.add(gene);
                    twoBit.close();

                    order ++;
                }
                ChrGeneMap.put(chr, geneList);
            }
        } catch (Exception ex) {
            System.out.println("select two bit Exception");
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
                for (Iterator<Gene> geneIterator = chrSelectGene.iterator(); geneIterator.hasNext();) {
                    Gene gene = geneIterator.next();
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
                for (Iterator<Gene> gene_iterator = GeneList.iterator(); gene_iterator.hasNext(); ) {
                    lower = upper;
                    upper = upper + gene_pm_range;
                    Gene gene = gene_iterator.next();
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
     * write the simulated mutation informations into file system
     */
    public void mutateGeneInfomation(String outputFile) {
        Set<String> geneIds = this.geneMutatedPosition.keySet();
        BufferedWriter bfw = null;
        try {
            bfw = new BufferedWriter(
                    new OutputStreamWriter(new FileOutputStream(new File(outputFile)))
            );
            int refCount, altCount;
            String refExonSeq, altExonSeq, writeOut;
            HashSet<Integer> mutPositions;
            Object[] positions;
            for (String geneId: geneIds) {
                refExonSeq = this.genesOriginExonSeq.get(geneId);
                altExonSeq = this.genesMutatedExonSeq.get(geneId);
                mutPositions = this.geneMutatedPosition.get(geneId);
                positions = mutPositions.toArray();
                Arrays.sort(positions);
                refCount = this.mutGeneRefAltCount.get(geneId)[0];
                altCount = this.mutGeneRefAltCount.get(geneId)[1];
                writeOut = geneId + ": " + Arrays.toString(positions) + "\n" + "ref: " + refExonSeq + "\n" +
                           "alt: " + altExonSeq + "\n" + "ref&alt reads count: " + refCount + "\t" + altCount;
                bfw.write(writeOut);
                bfw.newLine();
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
     * generate sequencing reads for each selected gene. The simulated sequencing reads will be write out into the
     * output file
     * @param fragmentMean fragment mean
     * @param fragmentTheta fragment std
     * @param inputOutputFile output file path
     * @param InputMultiple replicates for each experiment
     */
    public void generateReads(int fragmentMean, int fragmentTheta, String inputOutputFile, int InputMultiple) {
        BufferedWriter fw = null;
        try {
            fw = new BufferedWriter(
                    new OutputStreamWriter(new FileOutputStream(new File(inputOutputFile)))
            );
            // get mutated genes' ID
            Set<String> mutGeneIds = this.geneMutatedPosition.keySet();

            for (Map.Entry<String, LinkedList<Gene>> entry : ChrGeneMap.entrySet()) {
                LinkedList<Gene> GeneList = entry.getValue();
                Gene gene;

                for (Iterator<Gene> iterator = GeneList.iterator(); iterator.hasNext(); ) {
                    gene = iterator.next();
                    String geneId = gene.getGeneId();

                    // if gene in mutate gene list, change its exon sequence to variant sequence
                    if (mutGeneIds.contains(geneId)) {
                        gene.setExonSeq(this.genesMutatedExonSeq.get(geneId));
                    }

                    gene.calculateReadsCount(librarySize);

                    gene.generateReads(fragmentMean, fragmentTheta, fw, this.readLength, InputMultiple,
                                       this.geneMutatedPosition.getOrDefault(geneId, null), this.seqErrorModel);

                    if (mutGeneIds.contains(geneId)) {
                        int refCount = gene.getRefReadsCount();
                        int altCount = gene.getAltReadsCount();
                        int[] refAndAlt = new int[]{refCount, altCount};
                        this.mutGeneRefAltCount.put(geneId, refAndAlt);
                    }
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
    public void simulateSequencing(String dataPath, String vcfFile, int librarySize, int readLength, int maximumMut,
                                   int fragmentMean, int fragmentTheta, int mutGeneNum, int inputMultiple, int repeat, double pcrErrorProb) {

        this.seqErrorModel = new SequencingError(pcrErrorProb, readLength);
        this.setLibrarySize(librarySize);
        this.setReadLength(readLength);
        this.selectGene();

        // random select mutated genes
        SNPGenerator snpGenerator = new SNPGenerator(this.ChrGeneMap, mutGeneNum, vcfFile, maximumMut);
        this.genesOriginExonSeq = snpGenerator.getOriginExonSequence();
        this.genesMutatedExonSeq = snpGenerator.getMutatedExonSeqence();
        this.geneMutatedPosition = snpGenerator.getMutGenePosition();

        // generate a simplified GTF file
        this.createSimpleGTF(dataPath + "\\Homo_sapiens.GRCh38.sim.chr.gtf");
        // record the RPKM data
        this.transcriptionParameter(dataPath + "\\rpkm.txt");

        for (int i = 0; i < repeat; i++) {
            String InputOutputfile = dataPath + "\\Input"+i+".fasta";
            this.generateReads(fragmentMean,fragmentTheta, InputOutputfile,inputMultiple);
            mutateGeneInfomation(dataPath + "\\mutations"+i+".txt");
        }
    }

    public HashMap<String, HashSet<Integer>> getGeneMutatedPosition() {
        return this.geneMutatedPosition;
    }

    public HashMap<String, String> getGenesMutatedExonSeq() {
        return this.genesMutatedExonSeq;
    }
}
