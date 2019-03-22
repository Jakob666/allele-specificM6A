package SamtoolsPileupSNPCalling;

public class SamtoolsPileupSNPCalling {
    public static void main(String[] args) {
        String refGenomeFile = args[0];
        String fastqFile = args[1];
        int execthread = Integer.parseInt(args[2]);



        ReadsMapping.alignment(refGenomeFile, fastqFile, execthread);
        // TODO: complete this procedure
//        SamtoolsProcessing.samFileProcess(aligmentResultFile, outputDir, sraNum, samtools);
//        AseInference.inferenceASE(refGenomeFile, sortedBamFile, samtools);
//        SnpFilter sf = new SnpFilter(refVcfFilePath, rawVcfFile);
//        sf.filterVcf();

    }
}
