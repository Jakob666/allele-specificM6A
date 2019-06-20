package ExomeSeqSnpCalling;

import org.apache.log4j.Logger;

import java.io.*;
import java.util.ArrayList;

public class ReadsAlignment {

    private String refGenomeFile, bwa, samtools, picard;
    private int execthread;
    private Logger logger;

    public ReadsAlignment(String bwa, String samtools, String picard, String refGenomeFile, int execthread, Logger logger) {
        this.refGenomeFile = refGenomeFile;
        this.bwa = bwa;
        this.samtools = samtools;
        this.picard = picard;
        this.execthread = execthread;
        this.logger = logger;
    }

    /**
     * 为参考基因组建立索引
     */
    public void genomeIndex() {
        String cmd = String.join(" ", new String[]{this.bwa, "index", "-a bwtsw",this.refGenomeFile});
//        File bwtFile = new File(this.refGenomeFile+".bwt");
//        File pacFile = new File(this.refGenomeFile+".pac");
//        File annFile = new File(this.refGenomeFile+".ann");
//        File ambFile = new File(this.refGenomeFile+".amb");
//        File saFile = new File(this.refGenomeFile+".sa");
//        if (bwtFile.exists() & pacFile.exists() & annFile.exists() & ambFile.exists() & saFile.exists()) {
//            logger.debug("index already established");
//            return;
//        }
        this.logger.debug("index reference genome " + this.refGenomeFile);
        this.logger.debug(cmd);
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                this.logger.error("reference genome index failed");
                System.exit(2);
            }
        } catch (IOException | InterruptedException ie) {
            this.logger.error(ie.getMessage());
            System.exit(2);
        }
    }

    /**
     * reads align to reference genome and preprocess
     * @return marked duplicated PCR file
     */
    public String alignmentToGenome(String mateFile1, String mateFile2, String alignmentSamFile) {
        File shellScript = new File(new File(mateFile1).getParent(), "alignment.sh");
        String alignmentBamFile, sortedBamFile, markedBamFile = null;
        try {
            BufferedWriter bfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(shellScript)));
            readsAlignment(mateFile1, mateFile2, alignmentSamFile, bfw);
            alignmentBamFile = samToBam(alignmentSamFile, bfw);
            sortedBamFile = sortBamFile(alignmentBamFile, bfw);
            markedBamFile = markDuplicates(sortedBamFile, bfw);
            bfw.close();

            Process p = Runtime.getRuntime().exec("chmod a+x " + shellScript.getAbsolutePath());
            int res = p.waitFor();
            if (res != 0) {
                logger.error("change mode failed");
                System.exit(2);
            }
            p = Runtime.getRuntime().exec(shellScript.getAbsolutePath());
            res = p.waitFor();
            if (res != 0) {
                logger.error("alignment failed");
                System.exit(2);
            }
        } catch (IOException | InterruptedException ie) {
            logger.error(ie.getMessage());
            System.exit(2);
        }
        boolean succ = shellScript.delete();
        if (!succ) {
            logger.error("fail to remove redundant file: " + shellScript);
        }

        return markedBamFile;
    }

    /**
     * reads align to reference genome
     */
    public void readsAlignment(String mateFile1, String mateFile2, String alignmentSamFile, BufferedWriter bfw) {
        String cmd;
        // 如果是单端测序，mateFile2为null
        if (mateFile2 == null) {
            cmd = String.join(" ", new String[] {this.bwa, "mem -t", Integer.toString(this.execthread),
                              "-M", this.refGenomeFile, mateFile1, ">", alignmentSamFile});
        } else {
            cmd = String.join(" ", new String[] {this.bwa, "mem -t", Integer.toString(this.execthread),
                              "-M", this.refGenomeFile, mateFile1, mateFile2, ">", alignmentSamFile});
        }
        logger.debug(cmd);
        try {
            bfw.write("#!/bin/bash\n");
            bfw.write(cmd);
            bfw.newLine();
        } catch (IOException ie) {
            logger.error(ie.getMessage());
            System.exit(2);
        }
    }

    /**
     * 将比对结果的SAM文件转换为BAM文件
     * @param alignmentSamFile sam文件名
     * @return BAM文件名
     */
    public String samToBam(String alignmentSamFile, BufferedWriter bfw) {
        String alignmentBamFile = alignmentSamFile.substring(0, alignmentSamFile.lastIndexOf(".")) + ".bam";
        String cmd = String.join(" ", new String[] {samtools, "view", "-@", Integer.toString(execthread), "-bS", alignmentSamFile, ">", alignmentBamFile});
        logger.debug(cmd);
        try {
            bfw.write(cmd);
            bfw.newLine();
            bfw.write("rm -f " + alignmentSamFile);
            bfw.newLine();
        } catch (IOException ie) {
            logger.error(ie.getMessage());
            System.exit(0);
        }

        return alignmentBamFile;
    }

    /**
     * 比对文件排序
     * @param alignmentBamFile 比对BAM文件
     * @return 排序后的BAM文件
     */
    public String sortBamFile(String alignmentBamFile, BufferedWriter bfw) {
        String sortedBamFile = alignmentBamFile.substring(0, alignmentBamFile.lastIndexOf(".")) + "_sort.bam";
        String cmd = String.join(" ", new String[] {samtools, "sort -@", Integer.toString(execthread) ,
                                 alignmentBamFile, "-o", sortedBamFile});
        logger.debug(cmd);
        try {
            bfw.write(cmd);
            bfw.newLine();
            bfw.write("rm -f " + alignmentBamFile);
            bfw.newLine();
        } catch (IOException ie) {
            logger.error(ie.getMessage());
            System.exit(2);
        }

        return sortedBamFile;
    }

    /**
     * Mark PCR重复
     * @param sortedBamFile 排序后的BAM文件
     * @return 标记重复的BAM文件
     */
    public String markDuplicates(String sortedBamFile, BufferedWriter bfw) {
        String deduplicatedBamFile = sortedBamFile.substring(0, sortedBamFile.lastIndexOf("_")) + "_markDup.bam";
        String metrics_file = new File(new File(this.refGenomeFile).getParent(), "output.metrics").getAbsolutePath();
        String cmd;
        if (picard.equals("picard"))
            cmd = "picard";
        else
            cmd = "java -jar" + picard;

        cmd += " MarkDuplicates I=" + sortedBamFile + " O=" + deduplicatedBamFile + " CREATE_INDEX=true " +
                "VALIDATION_STRINGENCY=SILENT M="+ metrics_file +" TMP_DIR=./";

        logger.debug(cmd);
        try {
            bfw.write(cmd);
            bfw.newLine();
            bfw.write("rm -f " + sortedBamFile);
            bfw.newLine();
        } catch (IOException ie) {
            logger.error(ie.getMessage());
            System.exit(2);
        }

        return deduplicatedBamFile;
    }

    /**
     * 若同一个样本多组重复，将BAM文件进行合并
     * @param bamFileList BAM文件列表
     * @param mergedBamFile 合并的BAM文件
     */
    public void mergeBamFiles(ArrayList<String> bamFileList, String mergedBamFile) {
        File shellScript = new File(new File(mergedBamFile).getParent(), "merge.sh");
        String cmd = String.join(" ", new String[] {samtools, "merge", mergedBamFile, String.join(" ", bamFileList)});
        logger.debug(cmd);
        try {
            BufferedWriter bfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(shellScript)));
            bfw.write("#!/bin/bash\n");
            bfw.write(cmd);
            bfw.newLine();
            for (String bamFile: bamFileList) {
                bfw.write("rm -f " + bamFile.substring(0, bamFile.lastIndexOf("bam")) + "*");
                bfw.newLine();
            }
            bfw.close();
            Process process = Runtime.getRuntime().exec("chmod a+x " + shellScript.getAbsolutePath());
            int res = process.waitFor();
            if (res != 0) {
                logger.error("change mode failed");
                System.exit(2);
            }
            process = Runtime.getRuntime().exec(shellScript.getAbsolutePath());
            res = process.waitFor();
            if (res != 0) {
                logger.error("merge failed");
                System.exit(2);
            }
        } catch (IOException | InterruptedException ie) {
            logger.error(ie.getMessage());
            System.exit(2);
        }
        boolean succ = shellScript.delete();
        if (!succ) {
            logger.error("fail to remove redundant file: " + shellScript);
        }
    }
}
