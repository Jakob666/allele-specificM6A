package PeakSimulator;

import AseSeqSimulator.CommonMethod;
import AseSeqSimulator.TwoBitParser;

import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Peak2Seq {
    private HashMap<String, LinkedList<PeakRecord>> chrPeaks;
    private TwoBitParser twoBit;

    /**
     * get data information from peak calling bed file
     * @param bedFile bed File path
     */
    public Peak2Seq(String bedFile) {
        BedFileReader peakReader = new BedFileReader(bedFile);
        this.chrPeaks = peakReader.getChrPeaks();
    }

    /**
     * get genome sequence covered by m6A peak
     * @param twobitFile UCSC 2bit genome file
     * @throws Exception Exception
     */
    public void peakTrans2SeqString(String twobitFile)throws Exception{
        this.twoBit = new TwoBitParser(new File(twobitFile));
        int not_match = 0;
        int matchs = 0;
        for (Map.Entry<String, LinkedList<PeakRecord>> entry : chrPeaks.entrySet()) {
            String chr = entry.getKey();
            if (chr.equals("MT"))
                continue;
            LinkedList<PeakRecord> peakRecords = entry.getValue();
            this.twoBit.setCurrentSequence("chr"+chr);
            for(Iterator<PeakRecord> iterator = peakRecords.iterator(); iterator.hasNext();){
                PeakRecord peak = iterator.next();
                int start = peak.getPeakStart();
                int end  = peak.getPeakEnd();
                String strand = peak.getStrand();

                this.twoBit.reset();
                String fragment = this.twoBit.loadFragment(start-1,(end - start + 1));
                if (strand.equals("-")) {
                    fragment = CommonMethod.AntiChain(fragment);
                }

                // match m6A motif DRACH
                for(int i =0; i<fragment.length()-5; i++){
                    String motif = fragment.substring(i, i+5);
                    Pattern DRACH = Pattern.compile("[AGT][GA]AC[ACT]");
                    Matcher m = DRACH.matcher(motif);
                    if(m.matches()){
                        matchs++;
                    }else{
                        not_match++;
                    }
                }
            }
            twoBit.close();
        }
        System.out.println(matchs + "\t" + not_match);
    }
}
