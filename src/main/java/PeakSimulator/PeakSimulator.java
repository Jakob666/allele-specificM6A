package PeakSimulator;

import AseSeqSimulator.Fragmentation;
import AseSeqSimulator.Gene;
import GTFComponent.ElementRecord;
import org.apache.commons.math3.distribution.UniformRealDistribution;

import java.io.*;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

public class PeakSimulator {

    //设置Peak区域
    public static LinkedList<Peak> peakSimulation(int peak_num, int peakLength, String exonSeq, String chrNum,
                                                  String geneId, LinkedList<ElementRecord> exonList) {
        LinkedList<Peak> peakList;
        LinkedList<Integer> peakStartList = new LinkedList<Integer>();
        if (exonSeq.length() > peakLength) {
            // 起始位点的选取区间
            int start_range = exonSeq.length() - peakLength;
            int count = 0;
            while (peakStartList.size() < peak_num && count <= 100) {
                count++;
                int peak_start = (int) (Math.random() * start_range);
                int temp_length = peakLength;
                if(peak_start + peakLength > exonSeq.length()){
                    temp_length = exonSeq.length() - peak_start;
                }
                // 选取的peak不能相互重叠
                boolean coverage = false;
                for (Iterator<Integer> iterator = peakStartList.iterator(); iterator.hasNext(); ) {
                    int prev_start = iterator.next();
                    int prev_length = peakLength;
                    if(prev_start + peakLength > exonSeq.length()){
                        prev_length = exonSeq.length() - prev_start;
                    }
                    int start_max = Math.max(prev_start, peak_start);
                    int end_min = Math.min((prev_start + prev_length), (peak_start + temp_length));
                    if (start_max < end_min) {
                        coverage = true;
                        break;
                    }
                }
                // 如果当前生成的peak与其他peak不重叠，则将起始位点加入到peakStartList中
                if (!coverage) {
                    peakStartList.add(peak_start);
                }
            }
        }else{  // 如果设定的peak length比exon区域长度长，则将0加入到peakStartList中
            peakStartList.add(0);
        }

        peakList = PeakSimulator.peakLocation(exonList, peakStartList, peakLength, chrNum, exonSeq, geneId);

        return peakList;
    }

    //确定peak区域在染色体上的位置，将结果赋予Gene类对象的peakList属性
    private static LinkedList<Peak> peakLocation(LinkedList<ElementRecord> exonList, LinkedList<Integer> peakStartList,
                                                 int peakLength, String chr, String exonSeq, String geneID) {
        LinkedList<Peak> PeakList = new LinkedList<Peak>();
        try {
            for (Iterator<Integer> iterator = peakStartList.iterator(); iterator.hasNext(); ) {
                int accumulation = 0;
                int peak_start = iterator.next();
                int peak_end;
                // 由起始位点和peak length确定peak的终止位点
                if ((exonSeq.length() - peak_start) > peakLength) {
                    peak_end = peak_start + peakLength;
                } else {
                    peak_end = exonSeq.length() - 1;
                }
                // 截取peak的序列
                String PeakString = exonSeq.substring(peak_start, peak_end);
                Peak peak = new Peak();
                peak.setChr(chr);
                peak.setPeakString(PeakString);
                peak.setPeak_start(peak_start);
                peak.setPeak_end(peak_end);
                peak.setGeneID(geneID);
                boolean first_region = true;
                String region_string = "";
                for (Iterator<ElementRecord> exon_iterator = exonList.iterator(); exon_iterator.hasNext(); ) {
                    ElementRecord exon = exon_iterator.next();
                    int exon_start = exon.getElementStart();
                    int exon_end = exon.getElementEnd();
                    String strand = exon.getStrand();
                    int exon_length = exon_end - exon_start +1;
                    int prev_accumulation = accumulation;
                    accumulation = accumulation + exon_length;
                    int region_start;
                    int region_end;
                    String temp = "";
                    if (strand.equals("+")) {
                        peak.setStrand("+");
                        if (accumulation > peak_end) {
                            if (first_region) {
                                region_start = exon_start + peak_start - prev_accumulation;
                                region_end = exon_start + peak_end - prev_accumulation -1;
                            } else {
                                region_start = exon_start;
                                region_end = exon_start + peak_end - prev_accumulation -1;
                            }
                            peak.addRegion(region_start, region_end);
                            break;
                        } else if (accumulation > peak_start) {
                            if (first_region) {
                                region_start = exon_start + peak_start - prev_accumulation;
                                region_end = exon_end;
                            } else {
                                region_start = exon_start;
                                region_end = exon_end;
                            }
                            peak.addRegion(region_start, region_end);
                            first_region = false;
                        }
                    } else {
                        peak.setStrand("-");
                        if (accumulation > peak_end) {
                            if (first_region) {
                                region_start = exon_end - peak_end + prev_accumulation +1;
                                region_end = exon_end - peak_start + prev_accumulation;
                            } else {
                                region_start = exon_end - peak_end + prev_accumulation + 1;
                                region_end = exon_end;

                            }
                            peak.addRegion(region_start, region_end);
                            break;

                        } else if (accumulation > peak_start) {
                            if (first_region) {
                                region_start = exon_start;
                                region_end = exon_end - peak_start + prev_accumulation;
                            } else {
                                region_start = exon_start;
                                region_end = exon_end;
                            }
                            peak.addRegion(region_start, region_end);

                            first_region = false;
                        }
                    }
                }
                PeakList.add(peak);
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        return PeakList;
    }

    /**
     * write simulated peaks into bed file
     */
    public static void storeSimulatedPeaksRecord(HashMap<String, LinkedList<Gene>> selectedGenes, File simulatedPeakBedFile) {
        BufferedWriter bfw = null;
        try {
            bfw = new BufferedWriter(
                    new OutputStreamWriter(new FileOutputStream(simulatedPeakBedFile))
            );
            bfw.write("# chr\tgene ID\tPeak Region\tPeak Length\tPM\tRPKM\tInput Reads\tm6A reads");
            bfw.newLine();
            for (String chr: selectedGenes.keySet()) {
                LinkedList<Gene> genes = selectedGenes.get(chr);
                for (Gene gene: genes) {
                    String geneID = gene.getGeneId();
                    LinkedList<Peak> peakList = gene.getPeakList();
                    for (Peak peak: peakList) {
                        String peak_region = peak.getRegion();
                        int peakLength = peak.getPeakLength();
                        double pm = peak.getPM();
                        bfw.write(chr + "\t" + geneID + "\t" + peak_region + "\t" + peakLength + "\t" + pm + "\t" + peak.getControlString());
                        bfw.newLine();
                    }
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
     * trans the peak records into bed format record
     * @param peakfile output file of storeSimulatedPeaksRecord method
     * @param bedfile output BED file
     */
    public static void Trans2bed(File peakfile, File bedfile) {
        try {
            BufferedReader br = new BufferedReader(new FileReader(peakfile));
            FileWriter fw = new FileWriter(bedfile);
            br.readLine();//跳过表头
            fw.write("# chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickstart\tthickend\tItemRgb\tblockCount\tblockSizes\tblockStarts\n");
            while (br.ready()) {
                String strLine = br.readLine();
                String[] dataArr = strLine.split("\t");
                String chr = dataArr[0];
                String geneID = dataArr[1];
                String peakRegion = dataArr[2];
                String[] regionArr = peakRegion.split(";");
                double pm = Double.parseDouble(dataArr[4]);
                int peak_start = 0;
                int peak_end = 0;
                String blockStarts = "";
                String blockSizes = "";
                for (int i = 0; i < regionArr.length; i++) {
                    String[] range = regionArr[i].split("-");
                    int start = Integer.parseInt(range[0]);
                    int end = Integer.parseInt(range[1]);
                    if(peak_start ==0 || peak_start > start){
                        peak_start = start;
                    }
                    if(end > peak_end){
                        peak_end = end;
                    }
                }
                for(int i = 0; i < regionArr.length; i++){
                    String[] range = regionArr[i].split("-");
                    int start = Integer.parseInt(range[0]);
                    int end = Integer.parseInt(range[1]);
                    int blockStart = start-peak_start;
                    int blockSize = end - start + 1;
                    blockStarts = blockStarts + blockStart + ",";
                    blockSizes = blockSizes + blockSize + ",";
                }
                fw.write("chr" + chr + "\t" + peak_start + "\t" + peak_end + "\t" + geneID + "\t" + pm + "\t+\t" + peak_start + "\t" + peak_end + "\t0\t" + regionArr.length + "\t" + blockSizes + "\t" + blockStarts.substring(0,blockStarts.length()-1) + "\n");
            }
            br.close();
            fw.close();
        } catch (IOException io) {
            io.printStackTrace();
        }
    }

    /**
     * trans the peak records into bed format peak center format
     * @param peakfile output file of storeSimulatedPeaksRecord method
     * @param bedfile output BED file
     * @param centre peak center
     */
    public static void Trans2bed_peaksite(File peakfile, File bedfile, int centre){
        try {
            BufferedReader br = new BufferedReader(new FileReader(peakfile));
            FileWriter fw = new FileWriter(bedfile);
            br.readLine();//跳过表头
            fw.write("# chr\tStart\tEnd\tsite\tname\tscore\n");
            while (br.ready()) {
                String strLine = br.readLine();
                String[] dataArr = strLine.split("\t");
                String chr = dataArr[0];
                String geneID = dataArr[1];
                String peakRegion = dataArr[2];
                String[] regionArr = peakRegion.split(";");
                double pm = Double.parseDouble(dataArr[4]);
                String peak_string = null;
                int temp_length = 0;
                for (int i = 0; i < regionArr.length; i++) {
                    String[] range = regionArr[i].split("-");
                    int start = Integer.parseInt(range[0]);
                    int end = Integer.parseInt(range[1]);
                    int site = start + centre - temp_length -1;
                    temp_length =  temp_length + end - start +1;
                    if(temp_length >= centre) {
                        peak_string = "chr" + chr + "\t" + start + "\t" + end + "\t" + site + "\t" + geneID + "\t" + pm + "\n";
                        fw.write(peak_string);
                        break;
                    }
                }
            }
            br.close();
            fw.close();
        } catch (IOException io) {
            io.printStackTrace();
        }
    }
}
