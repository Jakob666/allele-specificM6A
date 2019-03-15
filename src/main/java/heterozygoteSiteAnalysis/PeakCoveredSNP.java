package heterozygoteSiteAnalysis;

import java.io.File;

public class PeakCoveredSNP {
    private File vcfFile;
    private File peakCallingRes;

    /**
     * Constructor
     * @param vcfRecordFile vcf record file
     * @param peakCallingRes m6a peak calling result
     */
    public PeakCoveredSNP(String vcfRecordFile, String peakCallingRes) {
        this.vcfFile = new File(vcfRecordFile);
        this.peakCallingRes = new File(peakCallingRes);
    }


}
