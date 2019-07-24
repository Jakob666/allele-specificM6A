package HierarchicalBayesianAnalysis;

import heterozygoteSiteAnalysis.PeakCoveredSNP;
import org.apache.log4j.Logger;

import java.io.File;

public class PeakCoveredSNPRecord {
    private String snpCallingFile, peakCallingFile, outputFile;
    private Logger logger;

    public PeakCoveredSNPRecord(String snpCallingFile, String peakCallingFile, String outputFile) {
        this.snpCallingFile = snpCallingFile;
        this.peakCallingFile = peakCallingFile;
        this.outputFile = outputFile;
        this.logger = initLog(new File(snpCallingFile).getParent());
    }

    /**
     * 将 peak calling与 SNP calling的结果整合，得到peak覆盖的SNP位点记录
     */
    public void getPeakCoveredSNP() {
        PeakCoveredSNP pcs = new PeakCoveredSNP(this.snpCallingFile, this.peakCallingFile, this.outputFile, this.logger);
        pcs.filterSNPAndPeak();
    }

    /**
     * 初始化log4j Logger 对象
     * @param logHome output directory of log file
     * @return Logger instance
     */
    private static Logger initLog(String logHome) {
        System.setProperty("log_home", logHome);
        return Logger.getLogger(PeakCoveredSNPRecord.class);
    }
}
