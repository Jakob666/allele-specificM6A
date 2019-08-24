package CommonThreadClass;

import meripSeqPeakCalling.PeakCaller;
import org.apache.log4j.Logger;

/**
 * 子线程运行Peak calling
 */
public class RunPeakCalling implements Runnable {
    private String ipBamFile, inputBamFile, gtfFile,  outputDir;
    private Logger logger;

    /**
     * Constructor
     * @param ipBamFile IP比对结果
     * @param inputBamFile INPUT比对结果
     * @param gtfFile GTF文件
     * @param outputDir 输出目录
     * @param logger log4j对象
     */
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
        PeakCaller pc = new PeakCaller(this.gtfFile, this.ipBamFile, this.inputBamFile, this.outputDir, experimentName, this.logger);
        pc.peakCalling();
    }
}
