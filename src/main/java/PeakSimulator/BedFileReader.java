package PeakSimulator;

import java.io.*;
import java.util.HashMap;
import java.util.LinkedList;

public class BedFileReader {
    private String bedFile;
    // if proportion is true, will generate methylation data
    private boolean proportion = false;
    private HashMap<String, LinkedList<PeakRecord>> chrPeaks = new HashMap<String, LinkedList<PeakRecord>>();

    public BedFileReader(String bedFile) {
        this.bedFile = bedFile;
        this.readBedFile();
    }

    public BedFileReader(String bedFile, boolean proportion){
        this.bedFile = bedFile;
        this.proportion = proportion;
        this.readBedFile();
    }

    /**
     * get bed file content, establish HashMap for peaks on each chromosome
     */
    private void readBedFile() {
        BufferedReader br = null;
        try {
            br = new BufferedReader(
                    new InputStreamReader(new FileInputStream(new File(this.bedFile)))
            );
            String line = "";
            String[] dataArr;
            String chr, strand, bamfileword;
            int chr_start, chr_end;
            double reduce_proportion;
            LinkedList<PeakRecord> PeakList;
            while (line != null) {
                line = br.readLine();
                if (line == null || line.startsWith("#"))
                    continue;
                dataArr = line.split("\t");
                chr = dataArr[0];
                PeakList = this.chrPeaks.containsKey(chr)? this.chrPeaks.get(chr):new LinkedList<PeakRecord>();
                chr_start = Integer.parseInt(dataArr[1]);
                chr_end = Integer.parseInt(dataArr[2]);
                strand = dataArr[5];
                bamfileword = chr + "\t" + chr_start + "\t" + chr_end + "\t" + dataArr[4];

                PeakRecord peak = new PeakRecord(chr, chr_start, chr_end, bamfileword, strand);
                if (this.proportion) {
                    reduce_proportion = Double.parseDouble(dataArr[dataArr.length-1]);
                    peak.setReduce_proportion(reduce_proportion);
                }

                if (dataArr.length < 10) {
                    peak.AddPeakRegion(chr_start,chr_end);
                } else if (Integer.parseInt(dataArr[9]) < 2) {
                    peak.AddPeakRegion(chr_start,chr_end);
                } else {
                    String[] peakstarts = dataArr[11].split(",");
                    String[] peaklengths = dataArr[10].split(",");
                    for (int i = 0; i < peakstarts.length; i++) {
                        int start = chr_start + Integer.parseInt(peakstarts[i]) - 1;
                        int end = start + Integer.parseInt(peaklengths[i]) - 1;
                        peak.AddPeakRegion(start,end);
                    }
                }
                PeakList.add(peak);
                this.chrPeaks.put(chr, PeakList);
            }
            br.close();
        } catch (IOException ie) {
            ie.printStackTrace();
        } finally {
            if (br != null) {
                try {
                    br.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    public HashMap<String, LinkedList<PeakRecord>> getChrPeaks() {
        return this.chrPeaks;
    }

}
