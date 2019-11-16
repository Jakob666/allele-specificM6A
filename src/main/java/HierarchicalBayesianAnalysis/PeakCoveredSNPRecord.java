package HierarchicalBayesianAnalysis;

import heterozygoteSiteAnalysis.PeakCoveredSNP;
import org.apache.log4j.Logger;

import java.io.File;

public class PeakCoveredSNPRecord {
    private String gtfFile, snpCallingFile, peakCallingFile, outputFile;
    private int readsInfimum;
    private Logger logger;

    public PeakCoveredSNPRecord(String gtfFile, String snpCallingFile, String peakCallingFile, String outputFile, int readsInfimum) {
        this.gtfFile = gtfFile;
        this.snpCallingFile = snpCallingFile;
        this.peakCallingFile = peakCallingFile;
        this.outputFile = outputFile;
        this.readsInfimum = readsInfimum;
        this.logger = initLog(new File(snpCallingFile).getParent());
    }

    /**
     * get m6A signal covered SNV sites
     */
    public void getPeakCoveredSNP() {
        PeakCoveredSNP pcs = new PeakCoveredSNP(this.gtfFile, this.snpCallingFile, this.peakCallingFile,
                                                this.outputFile, this.readsInfimum, this.logger);
        pcs.filterSNPAndPeak();
    }

    /**
     * initialize log4j Logger instance
     * @param logHome output directory of log file
     * @return Logger instance
     */
    private static Logger initLog(String logHome) {
        System.setProperty("log_home", logHome);
        return Logger.getLogger(PeakCoveredSNPRecord.class);
    }
}
