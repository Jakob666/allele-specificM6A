package CommonThreadClass;

import meripSeqPeakCalling.peakCaller;
import org.apache.log4j.Logger;

public class RunPeakCalling implements Runnable {
    private String ipBamFile, inputBamFile, gtfFile,  outputDir;
    private Logger logger;

    RunPeakCalling(String ipBamFile, String inputBamFile, String gtfFile, String outputDir, Logger logger) {
        this.inputBamFile = inputBamFile;
        this.ipBamFile = ipBamFile;
        this.gtfFile = gtfFile;
        this.outputDir = outputDir;
        this.logger = logger;
    }

    public void run() {
        getM6aPeak();
    }

    /**
     * get peak calling result
     */
    private void getM6aPeak() {
        String experimentName = "m6aPeak";
        peakCaller.peakCalling(gtfFile, outputDir, ipBamFile, inputBamFile, experimentName, logger);
    }
}
