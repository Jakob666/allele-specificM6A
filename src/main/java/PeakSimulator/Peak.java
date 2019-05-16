package PeakSimulator;

import AseSeqSimulator.Fragmentation;

import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.Map;

/**
 * model to describe m6A peak
 */
public class Peak {
    private String chr;
    private LinkedHashMap<Integer, Integer> RegionMap = new LinkedHashMap<Integer, Integer>();
    private String PeakString;
    private double PM;
    private int m6A_reads_count;
    private String controlString;
    private String GeneID;
    private String strand;
    private int peak_start;
    private int peak_end;

    public Peak() {

    }

    public void setChr(String chr) {
        this.chr = chr;
    }

    public void addRegion(int start, int end) {
        RegionMap.put(start, end);
    }

    public String getChr() {
        return chr;
    }

    public String getRegion() {
        String region_string = "";
        for (Map.Entry<Integer, Integer> entry : RegionMap.entrySet()) {
            int start = entry.getKey();
            int end = entry.getValue();
            region_string = region_string + start + "-" + end + ";";
        }
        return region_string;
    }

    public int getPeakLength() {
        return PeakString.length();
    }

    public void setPeakString(String peakString) {
        this.PeakString = peakString;
    }

    public String getPeakString() {
        return PeakString;
    }

//    public LinkedList<Fragmentation> getFragmentList() {
//        return FragmentList;
//    }

//    public void releaseFragmentList(){
//        this.FragmentList.clear();
//    }

    public void setPM(double PM) {
        this.PM = PM;
    }

    public double getPM() {
        return PM;
    }

    public void setGeneID(String geneID) {
        this.GeneID = geneID;
    }

    public String getGeneID() {
        return GeneID;
    }

    public String getControlString() {
        return controlString;
    }

    public void setStrand(String strand) {
        this.strand = strand;
    }

    public String getStrand() {
        return strand;
    }

    public void setPeak_start(int peak_start) {
        this.peak_start = peak_start;
    }

    public void setPeak_end(int peak_end){
        this.peak_end = peak_end;
    }

    public int getPeak_start(){
        return peak_start;
    }

    public int getPeak_end(){
        return peak_end;
    }

//    public void addFragments(Fragmentation fragment){
//        FragmentList.add(fragment);
//    }

    public void setM6A_reads_count(double RPKM, long library_size) {
        int Input_reads = (int) ((PeakString.length() * library_size * RPKM) / Math.pow(10, 9));
        this.m6A_reads_count = (int) (Input_reads * (PM / (1 - PM)));
        this.controlString = RPKM + "\t" + Input_reads + "\t" + m6A_reads_count;
    }

    public void setM6A_reads_count(double RPKM, long library_size, int background){
        this.m6A_reads_count = (int) (background * (PM / (1 - PM)));
        this.controlString = RPKM + "\t" + background + "\t" + m6A_reads_count;
    }

    public int getM6A_reads_count(){
        return m6A_reads_count;
    }


}
