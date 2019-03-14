package meripSeqPeakCalling;

import java.io.IOException;

public class peakCaller {
    /**
     * peak calling with metpeak in R script
     * @param gtfFilePath GTF file path
     * @param outputDirPath name of output result directory
     */
    public void peakCalling(String gtfFilePath, String outputDirPath, String IPFileDir, String INPUTFileDir, String metPeakScript, String cellLine) {

        String[] params = new String[]{"Rscript ./metpeak_calling.r", gtfFilePath, IPFileDir, INPUTFileDir, outputDirPath, cellLine};
        String cmd = String.join(" ", params);
        System.out.println(cmd);

        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int res = p.waitFor();
            if (res != 0) {
                throw new RuntimeException("Run peak calling failed");
            }
        } catch (IOException | InterruptedException ie) {
            ie.printStackTrace();
        }

    }

}
