package ReadsMapping;

import SamtoolsPileupSNPCalling.AseInference;
import SamtoolsPileupSNPCalling.SnpFilter;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ReadsMapping {

    /**
     * complete snp calling process
     * @param refGenomeFilePath reference sequence file path
     * @param gtfFile GTF file for generate genome index
     * @param refVcfFile reference VCF file, recommend dbSNP
     * @param prefix output result file prefix
     * @param fastqDir directory of the INPUT fastq file
     * @param outputDir output result directory
     * @param picardJar executable picard tool jar
     * @param gatkJar executable gatk tool jar
     * @param samtools executable samtools tool file
     * @param execThread number of working thread
     * @param logger Logger instance
     */
    public static void snpCalling(String refGenomeFilePath, String gtfFile, String refVcfFile, String prefix, String fastqDir,
                                  String outputDir, String picardJar, String gatkJar, String samtools, int execThread, Logger logger) {

        File refGenomePath = new File(refGenomeFilePath);
        String refGenomeDir = refGenomePath.getParent();
        int readLength = 0;

        File fastqDirectory = new File(fastqDir);
        File outputTargetDir = new File(outputDir);
        if (!outputTargetDir.exists()) {
            boolean res = outputTargetDir.mkdir();
            if (!res) {
                logger.error("fail to establish output directory " + outputDir);
                System.exit(2);
            }
        }

        File[] fastqFiles = fastqDirectory.listFiles();
        if (fastqFiles != null) {
            readLength = getFastqReadLength(fastqFiles[0].getAbsolutePath(), logger) - 1;

            logger.debug("start SNP calling");
            File outputSamFile = new File(refGenomeDir, "Aligned.out.sam");
            if (!outputSamFile.exists())
                readsMapping(refGenomeDir, refGenomeFilePath, gtfFile, fastqFiles, readLength, execThread, logger);
            String bamFile = filterMappedReads(refGenomeFilePath, picardJar, gatkJar, samtools, outputDir, prefix, logger);
            // get SNP reads abundance with samtools
            String readsCountFile = AseInference.inferenceASE(refGenomeFilePath, bamFile, samtools, logger);
            SnpFilter sf = new SnpFilter(refVcfFile, readsCountFile, logger);
            sf.filterVcf();
            logger.debug("SNP calling complete");
        }
    }

    /**
     * mapping RNA-Seq reads to reference genome using STAR
     * @param refGenomeDir directory of reference genome file
     * @param execThread number of thread which execute reads mapping process
     * @param fastqFiles File list of fastq file
     * @param genomeFileName name of reference genome file without path
     */
    public static void readsMapping(String refGenomeDir, String genomeFileName, String gtfFile, File[] fastqFiles,
                                    int readLength, int execThread, Logger log) {
        File genomeIdxFile = new File(refGenomeDir, "SAindex");
        log.debug("reads mapping start.");
        if (!genomeIdxFile.exists()) {
            genomeIndex(refGenomeDir, genomeFileName, gtfFile, execThread, log);
        }
        log.debug("1-pass alignment");
        alignment(refGenomeDir, fastqFiles, execThread, log);
        twoPassMapping(refGenomeDir, genomeFileName, readLength, execThread, log);
        log.debug("2-pass alignment");
        alignment(refGenomeDir, fastqFiles, execThread, log);
        log.debug("reads mapping complete.");
    }

    /**
     * adding read group information, sorting, marking duplicates and indexing. The reads are from SAM file output by readsMapping method
     * @param refGenomeFilePath reference genome file path
     * @param picardJar executable picard jar package
     * @param gatkJar executable gatk jar package
     * @param outputDir the output result directory which contains STAR alignment output file
     * @param prefix prefix of output result
     * @param log Logger instance
     */
    public static String filterMappedReads(String refGenomeFilePath, String picardJar, String gatkJar, String samtools,
                                         String outputDir, String prefix, Logger log) {
        String refGenomeDir = new File(refGenomeFilePath).getParent();
        File starSamFile = new File(refGenomeDir, "Aligned.out.sam");
        File sortedBamFile = new File(outputDir, prefix + "_sorted.bam");
        File deduplicatedBamFile = new File(outputDir, prefix + "_deduplicated.bam");
        File splitBamFile = new File(outputDir, prefix + "_split.bam");
        log.debug("filtering alignment reads");
        if (splitBamFile.exists()) {
            log.debug("already existed final bam file " + splitBamFile.getAbsolutePath());
            return splitBamFile.getAbsolutePath();
        }
        if (deduplicatedBamFile.exists()) {
            log.debug("already existed deduplicated bam file " + deduplicatedBamFile.getAbsolutePath());
            // 如果大于2Gb就不SplitNCigar操作了，容易挂
            if (deduplicatedBamFile.length() > new Long("2147483648")) {
                log.debug("Split Reads is not performed by default if the resulting file is too large.");
                return deduplicatedBamFile.getAbsolutePath();
            }

            refGenomeDict(picardJar, refGenomeFilePath, log);
            createFastaiFile(samtools, refGenomeFilePath, log);
            readsTrimReassign(gatkJar, refGenomeFilePath, deduplicatedBamFile.getAbsolutePath(), splitBamFile.getAbsolutePath(), log);
            log.debug("filtration complete");
            return splitBamFile.getAbsolutePath();
        }
        if (sortedBamFile.exists()) {
            log.debug("already existed sorted bam file " + sortedBamFile.getAbsolutePath());
            dropDuplicateReads(picardJar, sortedBamFile.getAbsolutePath(), deduplicatedBamFile.getAbsolutePath(),log);
            if (deduplicatedBamFile.length() > new Long("2147483648")) {
                log.debug("Split Reads is not performed by default if the resulting file is too large.");
                return deduplicatedBamFile.getAbsolutePath();
            }
            refGenomeDict(picardJar, refGenomeFilePath, log);
            createFastaiFile(samtools, refGenomeFilePath, log);
            readsTrimReassign(gatkJar, refGenomeFilePath, deduplicatedBamFile.getAbsolutePath(), splitBamFile.getAbsolutePath(), log);
            log.debug("filtration complete");
            return splitBamFile.getAbsolutePath();
        }
        readsGroup(picardJar, starSamFile.getAbsolutePath(), sortedBamFile.getAbsolutePath(), log);
        if (starSamFile.exists()) {
            String cmd = "rm -f " + starSamFile;
            try {
                Process p = Runtime.getRuntime().exec(cmd);
                int exitval = p.waitFor();
                if (exitval != 0) {
                    log.error("can not clean up redundant file " + starSamFile);
                }
            }catch (IOException | InterruptedException ie) {
                log.error(ie.getMessage());
            }
        }
        dropDuplicateReads(picardJar, sortedBamFile.getAbsolutePath(), deduplicatedBamFile.getAbsolutePath(),log);
        if (deduplicatedBamFile.length() > new Long("2147483648"))
            return deduplicatedBamFile.getAbsolutePath();

        refGenomeDict(picardJar, refGenomeFilePath, log);
        createFastaiFile(samtools, refGenomeFilePath, log);
        readsTrimReassign(gatkJar, refGenomeFilePath, deduplicatedBamFile.getAbsolutePath(), splitBamFile.getAbsolutePath(), log);
        log.debug("filtration complete");

        return splitBamFile.getAbsolutePath();
    }

    /**
     * Deprecated!
     * variant calling using gatk and output raw VCF file.
     * @param gatkLocalJar install directory which contains gatk jar package
     * @param genomeFileName reference genome file path
     * @param outputDir output result directory
     * @param prefix sra ID number of this data
     */
    public static void variantCalling(String gatkLocalJar, String genomeFileName, String outputDir, String prefix, Logger log) {
        File splitBamFile = new File(outputDir, prefix + "_split.bam");
        File outputVCF = new File(outputDir, prefix + "RawSNP.vcf");

        String cmd = gatkLocalJar + " HaplotypeCaller -R " + genomeFileName + " --dont-use-soft-clipped-bases true -I " +
                     splitBamFile.getAbsolutePath() + " -stand-call-conf 20.0 -O " + outputVCF.getAbsolutePath();
        log.debug("generate haplotype vcf file");
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                log.error("raw SNP calling failed");
                System.exit(2);
            }
        } catch (IOException | InterruptedException ie) {
            log.error(ie.getMessage());
            System.exit(2);
        }
    }

    /**
     * Deprecated!
     * variant calling using gatk and output raw VCF file
     * @param gatkLocalJar install directory which contains gatk jar package
     * @param genomeFileName reference genome file path
     * @param outputDir output result directory
     * @param prefix sra ID number of this data
     */
    public static String variantFilter(String gatkLocalJar, String genomeFileName, String outputDir, String prefix, Logger log) {
        File rawVcfFile = new File(outputDir, prefix + "RawSNP.vcf");
        File refinedVCF = new File(outputDir, prefix + "RefinedSNP.vcf");
        String cmd = gatkLocalJar + " VariantFiltration -R " + genomeFileName + " -V " + rawVcfFile.getAbsolutePath() +
                     " -window 35 -cluster 3 --filter-name FS -filter \"FS > 30.0\" --filter-name QD -filter \"QD < 2.0\" -O "
                     + refinedVCF.getAbsolutePath();
        log.debug("filter raw SNP sites");

        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                log.error("filter raw vcf failed");
                System.exit(2);
            }

            p = Runtime.getRuntime().exec("rm -f " + rawVcfFile.getAbsolutePath());
            exitVal = p.waitFor();
            if (exitVal != 0) {
                log.error("remove previous file failed.");
                System.exit(2);
            }

        } catch (IOException | InterruptedException ie) {
            log.error(ie.getMessage());
            System.exit(2);
        }

        return refinedVCF.getAbsolutePath();
    }

    /**
     * index the reference genome with STAR command "STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles hg19.fa --runThreadN"
     * this STAR command will generate two file named "Genome" and "SAindex" for further 1-pass mapping
     * @param refGenomeDir  directory of the reference genome file
     * @param genomeFilePath    reference genome file Path
     * @param execThread    number of working thread
     */
    public static void genomeIndex(String refGenomeDir, String genomeFilePath, String gtfFile, int execThread, Logger log) {
        String idxCmd = "STAR --runMode genomeGenerate --genomeDir " + refGenomeDir + " --genomeFastaFiles " + genomeFilePath +
                        " --runThreadN " + execThread + " --limitGenomeGenerateRAM=124544990592 --sjdbGTFfile " + gtfFile;

        log.debug("index reference genome file " + genomeFilePath);
        log.debug(idxCmd);
        List<String> processList = new ArrayList<String>();
        try{
            Process p = Runtime.getRuntime().exec(idxCmd);

            // waitFor() method return 0 if Process works well
            int exitVal = p.waitFor();
            HashMap<String, String> cmdOutput = getCommandOutput(p);
            String debugInfo = cmdOutput.get("debug");
            String errorInfo = cmdOutput.get("error");
            if (!debugInfo.equals(""))
                log.debug(debugInfo);
            if (!errorInfo.equals(""))
                log.error(errorInfo);
            p.destroy();

            if (exitVal != 0) {
                System.exit(2);
            }
        }catch (IOException e){
            e.printStackTrace();
        }catch (InterruptedException e) {
            log.error("Genome index failed.");
        }
    }

    /**
     * 1-pass mapping with STAR command "STAR --genomeDir $genomeDir --readFilesIn mate1.fq mate2.fq --runThreadN <n>"
     * this command generate Aligned.out.sam files which record the first-time's alignment result and SJ.out.tab file
     * @param refGenomeDir  directory of the index reference genome file
     * @param fastqFiles  File instance list of fastq files
     * @param execThread number of working thread
     * @param log Logger instance
     */
    public static void alignment(String refGenomeDir, File[] fastqFiles, int execThread, Logger log) {

        String cmd = "STAR --genomeDir " + refGenomeDir + " --readFilesIn ";
        String[] sb = new String[fastqFiles.length];
        for (int i = 0; i < fastqFiles.length; i++) {
            sb[i] = fastqFiles[i].getAbsolutePath();
        }
        cmd += String.join(",", sb) + " --runThreadN " + execThread;
        log.debug(cmd);
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            HashMap<String, String> cmdOutput = getCommandOutput(p);
            String debugInfo = cmdOutput.get("debug");
            String errorInfo = cmdOutput.get("error");
            if (!debugInfo.equals(""))
                log.debug(debugInfo);
            if (!errorInfo.equals(""))
                log.error(errorInfo);
            p.destroy();

            if (exitVal != 0) {
                log.error("reads mapping failed");
                System.exit(2);
            }
        }catch (IOException | InterruptedException ie){
            log.error(ie.getMessage());
        }
    }

    /**
     * 2-pass mapping with STAR command "STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles hg19.fa --sjdbFileChrStartEnd /path/to/1pass/SJ.out.tab --sjdbOverhang 75 --runThreadN <n>"
     * the command result in a new Genome, SA and SAindex file and sjdbList.out.tab file
     * @param refGenomeDir directory of the index reference genome file
     * @param genomeFileName reference genome file name
     * @param execThread number of working thread
     */
    public static void twoPassMapping(String refGenomeDir, String genomeFileName, int readLength, int execThread, Logger log) {

        File onePassSJOutput = new File(refGenomeDir, "SJ.out.tab");
        String cmd = "STAR --runMode genomeGenerate --genomeDir " + refGenomeDir + " --genomeFastaFiles " + genomeFileName +
                     " --sjdbFileChrStartEnd " + onePassSJOutput + " --sjdbOverhang " + readLength + " --runThreadN "
                     + execThread + " --limitGenomeGenerateRAM 124544990592";
        log.debug("2-pass mapping procedure...");
        log.debug(cmd);
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            HashMap<String, String> cmdOutput = getCommandOutput(p);
            String debugInfo = cmdOutput.get("debug");
            String errorInfo = cmdOutput.get("error");
            if (!debugInfo.equals(""))
                log.debug(debugInfo);
            if (!errorInfo.equals(""))
                log.error(errorInfo);
            p.destroy();

            if (exitVal != 0) {
                System.exit(2);
            }
        }catch (IOException e){
            e.printStackTrace();
        }catch (InterruptedException e) {
            log.error("2-pass mapping failed.");
        }
    }

    /**
     * adding read group information with picard tool
     * @param picardJar install directory which contains picard jar package
     * @param starSamFile the sam file output by STAR
     * @param sortedBamFile picard tool output sorted bam file
     */
    public static void readsGroup(String picardJar, String starSamFile, String sortedBamFile, Logger log) {
        String outputDir = new File(sortedBamFile).getParent();
        File tmp_dir = new File(outputDir, "tmp");
        if (!tmp_dir.exists())
            tmp_dir.mkdir();
        String cmd;
        if (picardJar.equals("picard"))
            cmd = "picard";
        else
            cmd = "java -jar " + picardJar;
        cmd += " AddOrReplaceReadGroups I=" + starSamFile + " O=" + sortedBamFile + " SO=coordinate RGID=id RGLB=library " +
                "RGPL=platform RGPU=machine RGSM=sample TMP_DIR=" + tmp_dir.getAbsolutePath();
        log.debug("grouping alignment reads");
        log.debug(cmd);
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            HashMap<String, String> cmdOutput = getCommandOutput(p);
            String debugInfo = cmdOutput.get("debug");
            if (!debugInfo.equals(""))
                log.debug(debugInfo);
            p.destroy();

            if (exitVal != 0) {
                log.error("picard reads grouping process failed.");
                System.exit(2);
            }
        }catch (IOException | InterruptedException ie){
            log.error(ie.getMessage());
        }
    }

    /**
     * drop duplicates and indexing reads
     * @param picardJar install directory which contains picard jar package
     * @param sortedBamFile picard tool output sorted bam file
     * @param deduplicatedBamFile picard tool output deduplicated bam file
     */
    public static void dropDuplicateReads(String picardJar, String sortedBamFile, String deduplicatedBamFile, Logger log) {
        String outputDir = new File(sortedBamFile).getParent();
        File tmp_dir = new File(outputDir, "tmp");
        if (!tmp_dir.exists())
            tmp_dir.mkdir();
        String cmd;
        if (picardJar.equals("picard"))
            cmd = "picard";
        else
            cmd = "java -jar" + picardJar;

        cmd += " MarkDuplicates I=" + sortedBamFile + " O=" + deduplicatedBamFile + " CREATE_INDEX=true " +
                "VALIDATION_STRINGENCY=SILENT M=output.metrics TMP_DIR=" + tmp_dir.getAbsolutePath();
        log.debug("marking duplicated alignment reads");
        log.debug(cmd);
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            HashMap<String, String> cmdOutput = getCommandOutput(p);
            String debugInfo = cmdOutput.get("debug");
            if (!debugInfo.equals(""))
                log.debug(debugInfo);
            if (exitVal != 0) {
                log.error("picard reads deduplicate process failed.");
                System.exit(2);
            }

            p = Runtime.getRuntime().exec("rm -f " + sortedBamFile);
            exitVal = p.waitFor();
            if (exitVal != 0) {
                log.error("remove previous file failed.");
                System.exit(2);
            }
            p.destroy();
        }catch (IOException | InterruptedException ie){
            log.error(ie.getMessage());
        }
    }

    /**
     *
     * @param gatkJar install directory which contains gatk jar package
     * @param refGenomeFile reference genome file path
     * @param dedupBamFile deduplicated bam file
     * @param splitBamFile split N cigar reads bam file
     * @param log Logger instance
     */
    public static void readsTrimReassign(String gatkJar, String refGenomeFile, String dedupBamFile,
                                         String splitBamFile, Logger log) {
        String cmd;
        if (gatkJar.equals("gatk"))
            cmd = "gatk SplitNCigarReads";
        else
            cmd = "java -jar" + gatkJar + " -T  SplitNCigarReads ";

        cmd = cmd + " -R " + refGenomeFile + " -I " + dedupBamFile + " -O " + splitBamFile;
        log.debug("split N cigar sequence, may take a long time");
        log.debug(cmd);
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            HashMap<String, String> cmdOutput = getCommandOutput(p);
            String debugInfo = cmdOutput.get("debug");
            if (!debugInfo.equals(""))
                log.debug(debugInfo);

            if (exitVal != 0) {
                log.error("gatk reads split process failed.");
                System.exit(2);
            }

            p = Runtime.getRuntime().exec("rm -f " + dedupBamFile);
            exitVal = p.waitFor();
            if (exitVal != 0) {
                log.error("remove previous file failed.");
                System.exit(2);
            }
            p.destroy();
        } catch (IOException | InterruptedException ie) {
            log.error(ie.getMessage());
            System.exit(2);
        }
    }

    /**
     * create sequence genome with picard tool for further gatk reads trimming
     * @param picardJarDir install directory which contains picard jar package
     * @param refGenomeFastaPath reference sequence file path
     */
    public static void refGenomeDict(String picardJarDir, String refGenomeFastaPath, Logger log) {
        File referenceFile = new File(refGenomeFastaPath);
        File refSequenceDir = referenceFile.getParentFile();
        int dot = referenceFile.getName().lastIndexOf('.');
        File refSequenceDict = new File(refSequenceDir.getAbsolutePath(), referenceFile.getName().substring(0, dot) + ".dict");
        if (refSequenceDict.exists()) {
            return;
        }

        String cmd = "java -jar " + picardJarDir + " CreateSequenceDictionary R=" + refGenomeFastaPath + " O=" +
                     refSequenceDict.getAbsolutePath();
        log.debug("create reference sequence dictionary");
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                log.error("create sequence dictionary failed");
                System.exit(2);
            }
            p.destroy();
        } catch (IOException | InterruptedException ie){
            log.error(ie.getMessage());
        }
    }

    /**
     * create fai file for further gatk trimming
     * @param refGenomeFastaFile reference sequence file path
     */
    public static void createFastaiFile(String samtool, String refGenomeFastaFile, Logger log) {
        File idxFile = new File(refGenomeFastaFile + ".fai");
        if (idxFile.exists()) {
            return;
        }
        String cmd = samtool + " faidx " + refGenomeFastaFile + " -o " + idxFile.getAbsolutePath();
        log.debug("generate fai file");
        log.debug(cmd);
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                log.error("create sequence idx failed");
                System.exit(2);
            }
            p.destroy();
        } catch (IOException | InterruptedException ie){
            log.error(ie.getMessage());
        }
    }

    /**
     * get sequence length from fastq data
     * @param fastqFile fastq file path
     * @return length of reads in fastq file
     */
    public static int getFastqReadLength(String fastqFile, Logger log) {
        int length = 0;
        try {
            BufferedReader reader = new BufferedReader(
                    new InputStreamReader(new FileInputStream(new File(fastqFile)))
            );

            String line = "";
            while (line.isEmpty()) {
                line = reader.readLine();
            }
            String[] splitInfo = line.split(" ");
            String lengthInfo = splitInfo[splitInfo.length -1];

            Pattern pattern = Pattern.compile("\\d+");
            Matcher lengthMatch = pattern.matcher(lengthInfo);

            while (lengthMatch.find()) {
                length = Integer.parseInt(lengthMatch.group());
            }

        } catch (FileNotFoundException e) {
            log.error("fastqFile : [" + fastqFile + "] not found");
            System.exit(2);
        } catch (IOException e) {
            log.error("IO error");
        }
        log.debug("get reads length in file " + fastqFile + ", length: " + length);
        return length;
    }

    /**
     * get command output infomation
     * @param process Process instance
     * @return HashMap
     */
    private static HashMap<String, String> getCommandOutput(Process process) {
        StringBuilder cmdOutputDebug = new StringBuilder();
        StringBuilder cmdOutputError = new StringBuilder();
        BufferedReader bfIn = new BufferedReader(new InputStreamReader(process.getInputStream()));
        BufferedReader bferr = new BufferedReader(new InputStreamReader(process.getErrorStream()));
        HashMap<String, String> output = new HashMap<>();

        String line = null;
        try {
            while ((line = bfIn.readLine()) != null) {
                cmdOutputDebug.append(line).append("\n");
            }


            while ((line = bferr.readLine()) != null) {
                cmdOutputError.append(line).append("\n");
            }
        } catch (IOException ie) {
            cmdOutputError.append(ie.getMessage()).append("\n");
        } finally {
            closeStream(bfIn);
            closeStream(bferr);
        }
        output.put("debug", cmdOutputDebug.toString());
        output.put("error", cmdOutputError.toString());
        return output;
    }

    /**
     * close input stream
     * @param stream input stream instance
     */
    private static void closeStream(Closeable stream) {
        if (stream != null) {
            try {
                stream.close();
            } catch (IOException ie) {
                ie.printStackTrace();
            }
        }
    }
}
