package AseSeqSimulator;

import java.io.*;
import java.util.HashMap;
import java.util.LinkedList;

public class VcfReader {
    private String vcfFile;
    private HashMap<String, LinkedList<VcfRecord>> chrVcfs = new HashMap<>();

    public VcfReader(String vcfFile) {
        this.vcfFile = vcfFile;
    }

    /**
     * get vcf positions of each chromosome
     */
    public void getChrVcfSites() {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(
                    new InputStreamReader(new FileInputStream(new File(this.vcfFile)))
            );
            String line = "";
            String[] lineInfo;
            String chr, id, ref, alt;
            int vcfPosition;
            while (line != null) {
                line = bfr.readLine();
                if (line != null && !line.startsWith("#")) {
                    lineInfo = line.split("\t");
                    chr = lineInfo[0];
                    id = lineInfo[2];
                    ref = lineInfo[3];
                    alt = lineInfo[4];
                    vcfPosition = Integer.parseInt(lineInfo[1]);
                    VcfRecord vr = new VcfRecord(chr, vcfPosition, id, ref, alt);
                    LinkedList<VcfRecord> vcfs = this.chrVcfs.getOrDefault(chr, new LinkedList<VcfRecord>());
                    vcfs.add(vr);
                    this.chrVcfs.put(chr, vcfs);
                }
            }
            bfr.close();
        } catch (IOException ie) {
            ie.printStackTrace();
        } finally {
            if (bfr != null) {
                try {
                    bfr.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    public HashMap<String, LinkedList<VcfRecord>> getChrVcfs() {
        return this.chrVcfs;
    }
}

class VcfRecord {
    public int vcfSite;
    public String chrNum, vcfId, ref, alt;

    public VcfRecord(String chrNum, int vcfSite, String vcfId, String ref, String alt) {
        this.chrNum = chrNum;
        this.vcfSite = vcfSite;
        this.vcfId = vcfId;
        this.ref= ref;
        this.alt =alt;
    }
}
