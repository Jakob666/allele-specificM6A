package PeakSimulator;

import java.util.LinkedList;

public class PeakRecord {
    private String chr, bamFileWord, strand;
    private LinkedList<int[]> peakRegionList = new LinkedList<>();
    private int peakStart, peakEnd;
    private double reduceProportion;

    public PeakRecord(String chr, int peakStart, int peakEnd, String bamFileWord, String strand) {
        this.chr = chr;
        this.peakStart = peakStart;
        this.peakEnd = peakEnd;
        this.bamFileWord = bamFileWord;
        this.strand = strand;
    }

    public void AddPeakRegion(int regionStart, int regionEnd){
        int[] region = new int[]{regionStart, regionEnd};
        this.peakRegionList.add(region);
    }

    public void setReduce_proportion(double reduceProportion){
        this.reduceProportion = reduceProportion;
    }

    public String getChr() {
        return this.chr;
    }

    public String getStrand() {
        return this.strand;
    }

    public String getBamFileWord() {
        return this.bamFileWord;
    }

    public LinkedList<int[]> getPeakRegionList() {
        return this.peakRegionList;
    }

    public int getPeakStart() {
        return this.peakStart;
    }

    public int getPeakEnd() {
        return this.peakEnd;
    }

    public double getReduceProportion() {
        return this.reduceProportion;
    }
}
