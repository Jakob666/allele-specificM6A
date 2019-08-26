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
        File bwtFile = new File(this.refGenomeFile+".bwt");
        File pacFile = new File(this.refGenomeFile+".pac");
        File annFile = new File(this.refGenomeFile+".ann");
        File ambFile = new File(this.refGenomeFile+".amb");
        File saFile = new File(this.refGenomeFile+".sa");
        if (bwtFile.exists() & pacFile.exists() & annFile.exists() & ambFile.exists() & saFile.exists()) {
            logger.debug("index already established");
            return;
        }
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
        String alignmentBamFile, sortedBamFile, markedBamFile = null;

        alignmentBamFile = readsAlignment(mateFile1, mateFile2, alignmentSamFile);
        sortedBamFile = sortBamFile(alignmentBamFile);
        markedBamFile = markDuplicates(sortedBamFile);

        return markedBamFile;
    }

    /**
     * reads align to reference genome
     */
    public String readsAlignment(String mateFile1, String mateFile2, String alignmentSamFile) {
        String cmd;
        String alignmentBamFile = alignmentSamFile.substring(0, alignmentSamFile.lastIndexOf(".")) + ".bam";
        File logFile = new File(new File(alignmentBamFile).getParentFile().getParent(), "logout.log");
        // 如果是单端测序，mateFile2为null
        if (mateFile2 == null) {
            cmd = String.join(" ", new String[] {this.bwa, "mem -t", Integer.toString(this.execthread),
                              "-M", this.refGenomeFile, mateFile1});
        } else {
            cmd = String.join(" ", new String[] {this.bwa, "mem -t", Integer.toString(this.execthread),
                              "-M", this.refGenomeFile, mateFile1, mateFile2});
        }
        cmd = cmd + " 1>"+ alignmentSamFile +" 2>" + logFile.getAbsolutePath();
        this.logger.debug(cmd);
        BufferedWriter bfw = null;
        String outputDir = new File(alignmentBamFile).getParent();
        File script = new File(outputDir, "alignment.sh");
        try {
            bfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(script)));
            bfw.write("#!/bin/bash\n");
            bfw.write(cmd);
            bfw.newLine();
            bfw.write(this.samtools + " view -bS -h "+alignmentSamFile +" -o " + alignmentBamFile);
            bfw.newLine();
            bfw.write("rm -f " + alignmentSamFile);
            bfw.newLine();
            bfw.close();
            Process p = Runtime.getRuntime().exec("chmod a+x " + script.getAbsolutePath());
            int res = p.waitFor();
            if (res != 0) {
                this.logger.error("chmod failed");
                System.exit(2);
            }
            p = Runtime.getRuntime().exec(script.getAbsolutePath());
            res = p.waitFor();
            if (res != 0) {
                this.logger.error("alignment failed");
                System.exit(2);
            }
            script.delete();
        } catch (IOException | InterruptedException ie) {
            logger.error(ie.getMessage());
            System.exit(2);
        } finally {
            if (bfw != null) {
                try {
                    bfw.close();
                } catch (IOException e) {
                    this.logger.error(e.getMessage());
                    System.exit(2);
                }
            }
        }

        return alignmentBamFile;
    }

    /**
     * 比对文件排序
     * @param alignmentBamFile 比对BAM文件
     * @return 排序后的BAM文件
     */
    public String sortBamFile(String alignmentBamFile) {
        String sortedBamFile = alignmentBamFile.substring(0, alignmentBamFile.lastIndexOf(".")) + "_sort.bam";
        String cmd;
        if (this.picard.equals("picard"))
            cmd = "picard";
        else
            cmd = "java -jar " + this.picard;
        cmd += " AddOrReplaceReadGroups I=" + alignmentBamFile + " O=" + sortedBamFile + " SO=coordinate RGID=id RGLB=library " +
                "RGPL=platform RGPU=machine RGSM=sample TMP_DIR=./tmp";
        this.logger.debug("grouping alignment reads");
        this.logger.debug(cmd);

        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int res = p.waitFor();
            if (res != 0) {
                this.logger.error("sorting BAM files failed");
                System.exit(2);
            }
            p = Runtime.getRuntime().exec("rm -f " + alignmentBamFile);
            res = p.waitFor();
            if (res != 0) {
                this.logger.error("delete bam file failed");
                System.exit(2);
            }
        } catch (IOException | InterruptedException ie) {
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
    public String markDuplicates(String sortedBamFile) {
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
            Process p = Runtime.getRuntime().exec(cmd);
            int res = p.waitFor();
            if (res != 0) {
                this.logger.error("marking duplications failed");
                System.exit(2);
            }
            p = Runtime.getRuntime().exec("rm -f " + sortedBamFile);
            res = p.waitFor();
            if (res != 0) {
                this.logger.error("delete bam file failed");
                System.exit(2);
            }
        } catch (IOException | InterruptedException ie) {
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
                bfw.write("rm -f " + bamFile.substring(0, bamFile.lastIndexOf("ba")) + "*");
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
