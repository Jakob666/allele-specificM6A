package AseSeqSimulator;

import GTFComponent.ElementRecord;
import GTFComponent.TranscriptRecord;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.UniformIntegerDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;

/**
 * random fragment a gene transcript sequence, the fragment length obeys normal distribution
 */
public class Gene {
    private String geneId, strand, geneName;
    private TranscriptRecord longestTranscriptRecord;
    private double rnaRPKM, dnaRPKM;
    private String chr, exonSeq;
    private int rnaReadsCount, dnaReadsCount;
    private double[] pmRange = new double[]{0.8, 0.9};
    private ElementRecord exonList = null;
    private HashMap<Integer, ArrayList<int[]>> m6aSiteFragments = new HashMap<>(), m6aSiteMutateFragments = new HashMap<>();
    // fragment which contains the reads cover SNV sites
    private HashMap<Integer, ArrayList<int[]>> inputMutateFragments = new HashMap<>(), ipMutateFragments = new HashMap<>();
    private HashMap<Integer, int[]> rnaMutationCoverage = new HashMap<>(), dnaMutationCoverage = new HashMap<>();
    // reads with same location on exons
    private ArrayList<int[]> inputFragment = new ArrayList<>(), ipFragment = new ArrayList<>(), m6aSiteNormalFragments = new ArrayList<>();

    public Gene(String geneId, String geneName, String strand, String chr) {
        this.geneId = geneId;
        this.strand = strand;
        this.geneName = geneName;
        this.chr = chr;
    }

    /**
     * get the longest transcript of a certain gene, and set it a property for a Gene instance
     * @param longestTranscriptRecord TranscriptRecord instance with the longest sequence length
     */
    public void setLongestTranscriptRecord(TranscriptRecord longestTranscriptRecord) {
        this.longestTranscriptRecord = longestTranscriptRecord;
    }

    public void setPmRange(double lower, double upper) {
        this.pmRange = new double[]{lower, upper};
    }

    public void setRnaRPKM(double RPKM) {
        this.rnaRPKM = RPKM;
    }

    public void setDnaRPKM(double RPKM) {
        this.dnaRPKM = RPKM;
    }

    public String getGeneId() {
        return this.geneId;
    }

    public String getGeneName() {
        return this.geneName;
    }

    public String getStrand() {
        return this.strand;
    }

    public String getChr() {
        return this.chr;
    }

    public int getRnaReadsCount() {
        return this.rnaReadsCount;
    }

    public int getDnaReadsCount() {
        return this.dnaReadsCount;
    }

    public double getDnaRPKM() {
        return this.dnaRPKM;
    }

    public ElementRecord getExonList() {
        return this.exonList;
    }

    public String getExonSeq() {
        return this.exonSeq;
    }

    public TranscriptRecord getLongestTranscriptRecord() {
        return this.longestTranscriptRecord;
    }

    public HashMap<Integer, int[]> getRnaMutationCoverage() {
        return this.rnaMutationCoverage;
    }

    public HashMap<Integer, int[]> getDnaMutationCoverage() {
        return this.dnaMutationCoverage;
    }

    /**
     * merge gene exonic region to obtain completed exon sequence
     * @param twoBit TwoBitParser instance
     */
    public void splicing(TwoBitParser twoBit) {
        HashMap<String, ElementRecord> elementsOnTranscript = longestTranscriptRecord.getElementList();
        ElementRecord transcriptExon = elementsOnTranscript.getOrDefault("exon", null);
        StringBuilder sb = new StringBuilder();
        if (transcriptExon != null)
            this.exonList = transcriptExon;
        try {
            while (transcriptExon != null) {
                String strand = transcriptExon.getStrand();
                int exonStart = transcriptExon.getElementStart();
                int exonEnd = transcriptExon.getElementEnd();
                twoBit.reset();
                String exonSeq = twoBit.loadFragment(exonStart-1, (exonEnd - exonStart + 1));
                if (strand.equals("-")) {
                    sb.append(CommonMethod.AntiChain(exonSeq));
                } else {
                    sb.append(exonSeq);
                }
                transcriptExon = transcriptExon.getNextElement();
            }
        } catch (IOException ie) {
            ie.printStackTrace();
        }
        this.exonSeq = sb.toString();
        sb = null;
    }

    /**
     * calculate reads count via RPKM and library size
     * @param librarySize library size
     */
    public void calculateReadsCountViaLibrarySize(long librarySize) {
        this.rnaReadsCount = (int) ((this.exonSeq.length() * librarySize * this.rnaRPKM) / Math.pow(10, 9));
    }

    /**
     * calculate reads count via sequence depth and library size
     * @param depth sequencing depth
     * @param readLength read length
     */
    public void calculateReadsCountViaSequencingDepth(int depth, int readLength) {
        this.dnaReadsCount = this.exonSeq.length() * depth / readLength;
    }

    /**
     * calculate reads count via m6A IP sample enrichment ratio
     * @param pm m6A enrichment ratio
     * @param backgroundSize background size
     */
    private int m6aReadsCount(double pm, int backgroundSize) {
        return  (int) (backgroundSize * (pm / (1 - pm)));
    }

    /**
     * randomly choose fragments to simulate INPUT data
     * @param fragmentMean average length of fragment
     * @param fragmentTheta std of fragment length
     */
    public void enrichInputFragment(int fragmentMean, int fragmentTheta, int readLength, Set<Integer> mutatePositions,
                                    String type) {
        int readsCount = (type.equals("rna"))? this.rnaReadsCount: this.dnaReadsCount;
        List<Integer> mutations = null;
        if (mutatePositions != null) {
            mutations = new ArrayList<>(mutatePositions);
            Collections.sort(mutations);
        }
        // randomly form fragment length
        NormalDistribution nordi = new NormalDistribution(fragmentMean, fragmentTheta);
        int curReadsCount = 0, fragmentLength, breakPoint, endPoint;
        boolean cover;
        while (curReadsCount < readsCount) {
            cover = false;
            // get the start and end point of the fragment on exonic region
            fragmentLength = Math.abs((int) nordi.sample());
            int[] breakAndEndPoint = this.getBreakEndPoint(fragmentLength);
            breakPoint = breakAndEndPoint[0];
            endPoint = breakAndEndPoint[1];

            // if the fragment cover SNV sites, record it
            if (mutations != null)
                cover = this.ifCoverMutateStie(breakPoint, endPoint, readLength, mutations, this.inputMutateFragments);
            if (!cover)
                this.inputFragment.add(breakAndEndPoint);
            curReadsCount++;
        }
        if (mutations != null)
            mutations.clear();
        nordi = null;
    }

    /**
     * randomly choose fragments to simulate IP data
     * @param fragmentMean average length of fragment
     * @param fragmentTheta std of fragment length
     * @param geneM6aSites m6A modification sites on gene
     */
    public void enrichIpFragment(int fragmentMean, int fragmentTheta, int readLength, Set<Integer> geneM6aSites,
                                 Set<Integer> mutatePositions, String type) {
        int readsCount = (type.equals("rna"))? rnaReadsCount: dnaReadsCount;
        List<Integer> mutations = null;
        if (mutatePositions != null) {
            mutations = new ArrayList<>(mutatePositions);
            Collections.sort(mutations);
        }
        List<Integer> m6aSites = new ArrayList<>(geneM6aSites);
        Collections.sort(m6aSites);
        NormalDistribution nordi = new NormalDistribution(fragmentMean, fragmentTheta);
        int curReadsCount = 0, fragmentLength, break_point, end_point;
        boolean cover;
        while (curReadsCount < readsCount) {
            cover = false;
            // get the start and end point of the fragment on exonic region
            fragmentLength = Math.abs((int) nordi.sample());
            int[] breakAndEndPoint = this.getBreakEndPoint(fragmentLength);
            break_point = breakAndEndPoint[0];
            end_point = breakAndEndPoint[1];

            // if the fragment cover SNV sites, record it
            if (mutations != null)
                cover = this.ifCoverMutateStie(break_point, end_point, readLength, mutations, this.ipMutateFragments);
            if (!cover)
                this.ipFragment.add(breakAndEndPoint);

            // if the fragment cover m6A sites, put it into IP data
            if (geneM6aSites.size() != 0)
                this.ifCoverM6aSite(break_point, end_point, m6aSites);
            curReadsCount++;
        }
        m6aSites.clear();
        if (mutations != null)
            mutations.clear();
        nordi = null;
    }

    /**
     * form simulated sequencing file
     * @param mateFile1 output file 1
     * @param mateFile2 output file 2(only use when form pair-end data; otherwise, null)
     * @param readLength sequencing read length
     * @param fragmentMean average length of fragment
     * @param fragmentTheta std of fragment length
     * @param direct single-end or pair-end
     * @param mutExonSeq exon mutation sequence if the gene contains SNV sites; otherwise null
     * @param refProp major allele frequency
     * @param seqErrorModel sequencing error
     * @param sample ip or input
     * @param geneM6aAsm if gene contains m6A sites, the ASM ratio for each site
     */
    public void generateReads(BufferedWriter mateFile1, BufferedWriter mateFile2, int readLength, int fragmentMean,
                              int fragmentTheta, double refProp, String mutExonSeq, String direct, SequencingError seqErrorModel,
                              String sample, String type, HashMap<Integer, Double> geneM6aAsm, HashMap<Integer, Boolean> minorBias) {
        ArrayList<int[]> fragmentRanges;
        ArrayList<Integer> m6aSites = null;
        if (geneM6aAsm != null) {
            m6aSites = new ArrayList<>(geneM6aAsm.keySet());
            Collections.sort(m6aSites);
        }
        HashMap<Integer, ArrayList<int[]>> mutateFragmentRanges;
        fragmentRanges = (sample.equals("ip"))? this.ipFragment : this.inputFragment;
        mutateFragmentRanges = (sample.equals("ip"))? this.ipMutateFragments : this.inputMutateFragments;
        HashMap<Integer, int[]> coverageRecords = null;
        if (!sample.equals("ip"))
            coverageRecords = (type.equals("rna"))? rnaMutationCoverage: dnaMutationCoverage;

        try {
            // fragment reads without SNP sites
            this.writeIn(fragmentRanges, direct, readLength, seqErrorModel, mateFile1, mateFile2, this.exonSeq, "ref");

            // fragment reads with SNP sites
            ArrayList<int[]> mutateFragments;
            int count, majorAlleleCount;
            List<int[]> majorAlleleFragmentRanges, minorAlleleFragmentRanges;
            for (Integer mutateSite: mutateFragmentRanges.keySet()) {
                // SNV-covered fragment start and end position
                mutateFragments = mutateFragmentRanges.get(mutateSite);
                if (sample.equals("ip") && m6aSites != null) {
                     int site = this.ifInM6aPeakRange(mutateSite, m6aSites, fragmentTheta+2*fragmentTheta);
                     if (site != -1)
                         refProp = geneM6aAsm.get(site);
                }
                count = mutateFragments.size();
                if (Math.abs(refProp - 0.5) < 0.00001) {
                    if (count % 2 == 0)
                        majorAlleleCount = count / 2;
                    else
                        majorAlleleCount = (count + 1) / 2;
                    if (coverageRecords != null)
                        coverageRecords.put(mutateSite, new int[] {majorAlleleCount, majorAlleleCount});
                } else {
                    majorAlleleCount = (int) Math.ceil(count * refProp);
                    if (coverageRecords != null)
                        coverageRecords.put(mutateSite, new int[] {majorAlleleCount, count-majorAlleleCount});
                }
                boolean bias = false;
                if (sample.equals("ip") && m6aSites != null) {
                    int site = this.ifInM6aPeakRange(mutateSite, m6aSites, fragmentMean + 2 * fragmentTheta);
                    if (site != -1) {
                        refProp = geneM6aAsm.get(site);
                        bias = minorBias.get(site);
                        if (!bias)
                            majorAlleleCount = count - majorAlleleCount;
                    }
                }
                Collections.shuffle(mutateFragments);
                majorAlleleFragmentRanges = mutateFragments.subList(0, majorAlleleCount);
                minorAlleleFragmentRanges = mutateFragments.subList(majorAlleleCount, mutateFragments.size());
                this.writeIn(majorAlleleFragmentRanges, direct, readLength, seqErrorModel, mateFile1, mateFile2, this.exonSeq, "maj");
                majorAlleleFragmentRanges = null;
                this.writeIn(minorAlleleFragmentRanges, direct, readLength, seqErrorModel, mateFile1, mateFile2, mutExonSeq, "min");
                minorAlleleFragmentRanges = null;
            }

        } catch (Exception io) {
            io.printStackTrace();
        }
    }

    /**
     * enrichment reads under m6a peak
     * @param mateFile1 output file 1
     * @param mateFile2 output file 2(only use when form pair-end data; otherwise, null)
     * @param readLength sequencing read length
     * @param fragmentMean average length of fragment
     * @param fragmentTheta std of fragment length
     * @param direct single-end or pair-end
     * @param mutExonSeq exon mutation sequence if the gene contains SNV sites; otherwise null
     * @param refProp major allele frequency
     * @param seqErrorModel sequencing error
     */
    public void peakFragmentFromBackground(BufferedWriter mateFile1, BufferedWriter mateFile2, int readLength, int fragmentMean,
                                           int fragmentTheta, double refProp, String mutExonSeq, String direct,
                                           SequencingError seqErrorModel, Set<Integer> geneMutateSites,
                                           HashMap<Integer, Double> geneM6aAsm, HashMap<Integer, Boolean> minorBias) {
        List<Integer> mutations = null;
        if (geneMutateSites != null) {
            mutations = new ArrayList<>(geneMutateSites);
            Collections.sort(mutations);
        }
        ArrayList<Integer> m6aSites = null;
        if (geneM6aAsm != null) {
            m6aSites = new ArrayList<>(geneM6aAsm.keySet());
            Collections.sort(m6aSites);
        }
        int peakReadsCount, peakFragmentCount;
        UniformRealDistribution urd = new UniformRealDistribution(this.pmRange[0], this.pmRange[1]);
        ArrayList<int[]> fragmentRanges;
        List<int[]> randomFragmentRanges = new ArrayList<>();
        try {
            for (Integer m6aSite: this.m6aSiteFragments.keySet()) {
                fragmentRanges = this.m6aSiteFragments.getOrDefault(m6aSite, null);
                if (fragmentRanges == null)
                    continue;
                peakFragmentCount = fragmentRanges.size();
                // calculate enrichment reads number on each m6A sites and form reads from m6A covered fragments
                peakReadsCount = this.m6aReadsCount(urd.sample(), peakFragmentCount);
                for (int i = 0; i < peakReadsCount; i++) {
                    int[] range = fragmentRanges.get((int) (Math.random() * fragmentRanges.size()));
                    randomFragmentRanges.add(range);
                }

                // if the reads covered SNV sites, record it
                boolean cover;
                for (int[] fragmentRange: randomFragmentRanges) {
                    cover = false;
                    if (mutations != null)
                        cover = this.ifCoverMutateStie(fragmentRange[0], fragmentRange[1], readLength, mutations,
                                                       this.m6aSiteMutateFragments);
                    if (!cover)
                        this.m6aSiteNormalFragments.add(fragmentRange);
                }
                randomFragmentRanges.clear();
                // form reads without SNV sites
                this.writeIn(this.m6aSiteNormalFragments, direct, readLength, seqErrorModel, mateFile1, mateFile2, this.exonSeq, "ref");
                // form reads with SNV sites
                ArrayList<int[]> mutateFragments;
                int count, majorAlleleCount;
                List<int[]> majorAlleleFragmentRanges, minorAlleleFragmentRanges;
                for (Integer mutateSite: this.m6aSiteMutateFragments.keySet()) {
                    mutateFragments = this.m6aSiteMutateFragments.get(mutateSite);
                    boolean bias = false;
                    // judge if the SNV site is covered by m6A signal range [site-fragmentMean-2*theta, site+fragmentMean+2*theta]
                    if (m6aSites != null) {
                        int site = this.ifInM6aPeakRange(mutateSite, m6aSites, fragmentMean+2*fragmentTheta);
                        // if SNV covered by m6A signal, get simulated ASM ratio
                        if (site != -1) {
                            refProp = geneM6aAsm.get(site);
                            bias = minorBias.get(site);
                        }
                    }
                    count = mutateFragments.size();
                    if (bias)
                        majorAlleleCount = (int) (count * refProp);
                    else
                        majorAlleleCount = count - (int) (count * refProp);
                    Collections.shuffle(mutateFragments);
                    majorAlleleFragmentRanges = mutateFragments.subList(0, majorAlleleCount);
                    minorAlleleFragmentRanges = mutateFragments.subList(majorAlleleCount, mutateFragments.size());
                    this.writeIn(majorAlleleFragmentRanges, direct, readLength, seqErrorModel, mateFile1, mateFile2, this.exonSeq, "maj");
                    majorAlleleFragmentRanges = null;
                    this.writeIn(minorAlleleFragmentRanges, direct, readLength, seqErrorModel, mateFile1, mateFile2, mutExonSeq, "min");
                    minorAlleleFragmentRanges = null;
                }
                this.m6aSiteNormalFragments.clear();
                this.m6aSiteMutateFragments.clear();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        this.m6aSiteFragments.clear();
    }

    /**
     * if the reads cover m6A modification site
     * @param start fragment start position
     * @param end fragment end position
     * @param sites m6A modification sites
     */
    private void ifCoverM6aSite(int start, int end, List<Integer> sites) {
        ArrayList<int[]> fragmentRange;
        for (Integer site: sites) {
            if (start <= site && site <= end) {
                fragmentRange = this.m6aSiteFragments.getOrDefault(site, new ArrayList<>());
                fragmentRange.add(new int[]{start, end});
                this.m6aSiteFragments.put(site, fragmentRange);
                break;
            }
        }
    }

    /**
     * if the reads cover mutation site
     * @param start fragment start position
     * @param end fragment end position
     * @param sites fragment end position
     * @return true, if covered; otherwise, false
     */
    private boolean ifCoverMutateStie(int start, int end, int readLength, List<Integer> sites,
                                      HashMap<Integer, ArrayList<int[]>> fragments) {
        ArrayList<int[]> fragmentRange;
        int readEnd = ((start + readLength -1) > this.exonSeq.length())? end: (start + readLength);
        for (Integer site: sites) {
            if (start <= site && site <= readEnd) {
                fragmentRange = fragments.getOrDefault(site, new ArrayList<>());
                fragmentRange.add(new int[]{start, end});
                fragments.put(site, fragmentRange);
                return true;
            }
        }
        return false;
    }

    /**
     * if the mutated read is covered by ASM m6A signal
     * @param mutatePosition mutation site position
     * @param m6aSites m6A modification sites
     * @param fragmentLength fragment length
     * @return true, if covered; otherwise, false
     */
    private int ifInM6aPeakRange(int mutatePosition, ArrayList<Integer> m6aSites, int fragmentLength) {
        for (Integer site: m6aSites) {
            if (Math.max(0, site-fragmentLength) <= mutatePosition && mutatePosition <= site + fragmentLength)
                return site;
        }
        return -1;
    }

    /**
     * get break and end point of fragment on exon sequence, according to the random generated position
     * @param fragmentLength fragment length
     * @return break and end point
     */
    private int[] getBreakEndPoint(int fragmentLength) {
        int break_point, end_point;
        int exonSeqLength = this.exonSeq.length();
        if (fragmentLength > exonSeqLength) {
            break_point = 0;
            end_point = exonSeqLength;
        } else {
            UniformIntegerDistribution uniformdi = new UniformIntegerDistribution(0, exonSeqLength - fragmentLength);
            break_point = uniformdi.sample();
            end_point = break_point + fragmentLength;
            uniformdi = null;
        }

        return new int[]{break_point, end_point};
    }

    private void writeIn(List<int[]> fragmentRanges, String direct, int readLength, SequencingError seqErrorModel,
                         BufferedWriter mateFile1, BufferedWriter mateFile2, String exonSequence, String type) {
        String fragmentString;
        Fragmentation fragment;
        int break_point, endPoint;
        for (int[] fragmentRange: fragmentRanges) {
            // fragment start and end position
            break_point = fragmentRange[0];
            endPoint = break_point + readLength - 1;

            if (endPoint > exonSequence.length())
                endPoint = exonSequence.length() - 1;
            fragmentString = exonSequence.substring(break_point, endPoint+1);
            fragment = this.getSequencingReads(fragmentString, break_point, endPoint, readLength, direct);
            if (direct.equals("SE"))
                this.writeReadInFile(fragment, seqErrorModel, mateFile1, readLength, break_point, endPoint, type);
            else
                this.pairReadToFile(fragment, seqErrorModel, mateFile1, mateFile2, readLength, break_point, endPoint);
            fragment = null;
            fragmentString = null;
        }
    }

    /**
     * single-end read into file
     * @param fragment Fragmentation instance
     * @param seqErrorModel sequencing error model instance
     * @param fw BufferedWriter instance
     * @param readLength read length
     * @param break_point fragment start position on exon
     * @param end_point fragment end position on exon
     */
    private void writeReadInFile(Fragmentation fragment, SequencingError seqErrorModel, BufferedWriter fw,
                                 int readLength, int break_point, int end_point, String type) {
        try {
            String sequencingRead = fragment.getSingleEndRead();
            sequencingRead = seqErrorModel.pcrErrorReads(sequencingRead);
            if (sequencingRead.length() != readLength) {
                int baseNNum = readLength - sequencingRead.length();
                String baseNSeq = this.fillBaseN(baseNNum);
                sequencingRead = sequencingRead + baseNSeq;
            }

            fw.write(">chr" + this.chr + "_" + this.geneId + "_" +break_point + "_" + end_point + " " + type);
            fw.newLine();
            fw.write(sequencingRead);
            fw.newLine();
        } catch (IOException | StringIndexOutOfBoundsException e) {
            e.printStackTrace();
        }
    }

    /**
     * pair-end read into file
     * @param fragment Fragmentation instance
     * @param seqErrorModel sequencing error instance
     * @param mateFile1 mate file 1
     * @param mateFile2 mate file 2
     * @param readLength read length
     * @param break_point fragment start position on exon
     * @param end_point fragment end position on exon
     */
    private void pairReadToFile(Fragmentation fragment, SequencingError seqErrorModel, BufferedWriter mateFile1,
                                BufferedWriter mateFile2, int readLength, int break_point, int end_point) {
        try {
            String[] sequencingRead = fragment.getPairEndRead();
            String mate1 = sequencingRead[0];
            String mate2 = sequencingRead[1];
            mate1 = seqErrorModel.pcrErrorReads(mate1);
            mate2 = seqErrorModel.pcrErrorReads(mate2);
            if (mate1.length() != readLength) {
                int baseNNum = readLength - mate1.length();
                String baseNSeq = this.fillBaseN(baseNNum);
                mate1 = mate1 + baseNSeq;
                mate2 = baseNSeq + mate2;
            }

            mateFile1.write(">chr" + this.chr + "_" + this.geneId + "_" +break_point + "_" + end_point);
            mateFile1.newLine();
            mateFile1.write(mate1);
            mateFile1.newLine();

            mateFile2.write(">chr" + this.chr + "_" + this.geneId + "_"+break_point + "_" + end_point);
            mateFile2.newLine();
            mateFile2.write(mate2);
            mateFile2.newLine();
        } catch (IOException | StringIndexOutOfBoundsException e) {
            e.printStackTrace();
        }
    }

    /**
     * use base N to fill the reads, when read length less than 50
     * @param baseNNum number of base N to fill
     * @return base N sequence
     */
    private String fillBaseN(int baseNNum) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < baseNNum; i++) {
            sb.append("N");
        }
        return sb.toString();
    }

    /**
     * generate fragment and form sequencing reads
     * @param break_point start position of the fragment on exon sequence
     * @param end_point end position of the fragment on exon sequence
     * @param readLength sequencing read length
     * @param direction "SE" single-end, "PE" pair-end
     * @return fragment instance
     */
    private Fragmentation getSequencingReads(String fragmentString, int break_point, int end_point, int readLength, String direction) {
        Fragmentation fragment = new Fragmentation(fragmentString, break_point, end_point);
        double randNum = Math.random();
        String strand = (randNum < 0.5)? "+" : "-";
        if (direction.equals("SE")) {   // single-end sequencing
            fragment.singleEndSequencingRead(strand, readLength);
        } else {    // pair-end sequencing
            fragment.pairEndSequencingRead(strand, readLength);
        }

        return fragment;
    }

    /**
     * release
     */
    public void release() {
        this.ipFragment.clear();
        this.inputFragment.clear();
        this.inputMutateFragments.clear();
        this.ipMutateFragments.clear();
    }
}

