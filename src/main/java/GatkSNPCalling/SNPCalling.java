package GatkSNPCalling;

import org.apache.log4j.Logger;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class SNPCalling {

    private static Logger log;

    /**
     * complete snp calling process
     * @param refGenomeFilePath reference sequence file path
     * @param fastqDir directory of the INPUT fastq file
     * @param outputDir output result directory
     * @param picardJarDir path of picard tool executive jar
     * @param gatkJarDir path of gatk tool executive file
     */
    public static void snpCalling(String refGenomeFilePath, String fastqDir, String outputDir, String picardJarDir, String gatkJarDir, int execThread, Logger logger) {

        log = logger;
        File refGenomePath = new File(refGenomeFilePath);
        String refGenomeDir = refGenomePath.getParent();
        int readLength = 0;

        File fastqDirectory = new File(fastqDir);
        File outputTargetDir = new File(outputDir);
        if (!outputTargetDir.exists()) {
            boolean res = outputTargetDir.mkdir();
        }

        File[] fastqFiles = fastqDirectory.listFiles();
        if (fastqFiles != null) {
            for (File f : fastqFiles) {
                if (f.getName().endsWith(".fastq")) {
                    readLength = getFastqReadLength(f.getAbsolutePath()) - 1;
                    int dot = f.getName().lastIndexOf('.');
                    String prefix = f.getName().substring(0, dot);

                    if (new File(outputDir, prefix + "_RefinedSNP.vcf").exists())
                        continue;
                    logger.debug(f.getAbsolutePath() + " start SNP calling");
                    readsMapping(refGenomeDir, refGenomePath.getAbsolutePath(), f.getAbsolutePath(), readLength, execThread, logger);
                    filterMappedReads(refGenomeFilePath, picardJarDir, outputDir, prefix);
                    refGenomeDict(picardJarDir, refGenomeFilePath);
                    createFastaiFile(refGenomeFilePath);
                    String bamFile = readsTrimReassign(gatkJarDir, refGenomeFilePath, outputDir, prefix);
                    variantCalling(gatkJarDir, refGenomeFilePath, outputDir, prefix);
                    String gatkVcfFile = variantFilter(gatkJarDir, refGenomeFilePath, outputDir, prefix);
                    readsCount(refGenomeFilePath, bamFile, gatkJarDir, gatkVcfFile);
                    logger.debug("SNP calling complete");
                }
            }
            if (fastqFiles.length > 1)
                mergeVcfFiles(outputDir, gatkJarDir);
        }

    }

    /**
     * mapping RNA-Seq reads to reference genome using STAR
     * @param refGenomeDir directory of reference genome file
     * @param execThread number of thread which execute reads mapping process
     * @param fastqFile absolute path of fastq file
     * @param genomeFileName name of reference genome file without path
     */
    public static void readsMapping(String refGenomeDir, String genomeFileName, String fastqFile, int readLength, int execThread, Logger log) {
        File genomeIdxFile = new File(refGenomeDir, "SAindex");
        log.debug("reads mapping start.");
        if (!genomeIdxFile.exists()) {
            genomeIndex(refGenomeDir, genomeFileName, execThread, log);
        }
        onePassMapping(refGenomeDir, fastqFile, execThread, log);
        twoPassMapping(refGenomeDir, genomeFileName, readLength, execThread, log);
        alignment(refGenomeDir, fastqFile, execThread, log);
        log.debug("reads mapping complete.");
    }

    /**
     * adding read group information, sorting, marking duplicates and indexing. The reads are from SAM file output by readsMapping method
     * @param picardJarDir install directory which contains picard jar package
     * @param outputDir the output result directory which contains STAR alignment output file
     */
    public static void filterMappedReads(String refGenomeFilePath, String picardJarDir, String outputDir, String prefix) {
        String refGenomeDir = new File(refGenomeFilePath).getParent();
        File starSamFile = new File(refGenomeDir, "Aligned.out.sam");
        File sortedBamFile = new File(outputDir, prefix + "_sorted.bam");
        File deduplicatedBamFile = new File(outputDir, prefix + "_deduplicated.bam");
        log.debug("filtering alignment reads");
        readsGroup(picardJarDir, starSamFile.getAbsolutePath(), sortedBamFile.getAbsolutePath());
        dropDuplicateReads(picardJarDir, sortedBamFile.getAbsolutePath(), deduplicatedBamFile.getAbsolutePath());
        log.debug("filtration complete");
    }

    /**
     *
     * @param gatkLocalJar install directory which contains gatk jar package
     * @param refGenomeFile reference genome file path
     * @param outputDir output result directory
     */
    public static String readsTrimReassign(String gatkLocalJar, String refGenomeFile, String outputDir, String prefix) {
        File deduplicatedBamFile = new File(outputDir, prefix + "_deduplicated.bam");
        File trimmedBamFile = new File(outputDir, prefix + "_trimmed.bam");
        String cmd = gatkLocalJar + " SplitNCigarReads -R " + refGenomeFile + " -I " + deduplicatedBamFile.getAbsolutePath() +
                     " -O " + trimmedBamFile.getAbsolutePath();
        log.debug("split N cigar sequence");

        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                log.error("gatk reads trim process failed.");
                System.exit(2);
            }

            p = Runtime.getRuntime().exec("rm -f " + deduplicatedBamFile.getAbsolutePath());
            exitVal = p.waitFor();
            if (exitVal != 0) {
                log.error("remove previous file failed.");
                System.exit(2);
            }
        } catch (IOException | InterruptedException ie) {
            ie.printStackTrace();
        }

        return trimmedBamFile.getAbsolutePath();
    }

    /**
     * variant calling using gatk and output raw VCF file
     * @param gatkLocalJar install directory which contains gatk jar package
     * @param genomeFileName reference genome file path
     * @param outputDir output result directory
     * @param prefix sra ID number of this data
     */
    public static void variantCalling(String gatkLocalJar, String genomeFileName, String outputDir, String prefix) {
        File trimmedBamFile = new File(outputDir, prefix + "_trimmed.bam");
        File outputVCF = new File(outputDir, prefix + "RawSNP.vcf");
        String cmd = gatkLocalJar + " HaplotypeCaller -R " + genomeFileName + " -I " +
                     trimmedBamFile.getAbsolutePath() + " -O " + outputVCF.getAbsolutePath();
        log.debug("generate haplotype vcf file");

        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                log.error("raw SNP calling failed");
                System.exit(2);
            }
        } catch (IOException | InterruptedException ie) {
            ie.printStackTrace();
        }
    }

    /**
     * variant calling using gatk and output raw VCF file
     * @param gatkLocalJar install directory which contains gatk jar package
     * @param genomeFileName reference genome file path
     * @param outputDir output result directory
     * @param prefix sra ID number of this data
     */
    public static String variantFilter(String gatkLocalJar, String genomeFileName, String outputDir, String prefix) {
        File rawVcfFile = new File(outputDir, prefix + "RawSNP.vcf");
        File refinedVCF = new File(outputDir, prefix + "RefinedSNP.vcf");
        String cmd = gatkLocalJar + " VariantFiltration -R " + genomeFileName + " -V " + rawVcfFile.getAbsolutePath() +
                     " -window 35 -cluster 3 -O " + refinedVCF.getAbsolutePath();
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
            ie.printStackTrace();
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
    public static void genomeIndex(String refGenomeDir, String genomeFilePath, int execThread, Logger log) {
        String idxCmd = "STAR --runMode genomeGenerate --genomeDir " + refGenomeDir + " --genomeFastaFiles " + genomeFilePath +
                        " --runThreadN " + execThread + " --limitGenomeGenerateRAM=124544990592";

        log.debug("index reference genome file " + genomeFilePath);
        List<String> processList = new ArrayList<String>();
        try{
            Process p = Runtime.getRuntime().exec(idxCmd);
            // return result of command can be obtained by using getInputStream() method in class Process
            BufferedReader input = new BufferedReader(new InputStreamReader(p.getInputStream()));

            // waitFor() method return 0 if Process works well
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                System.exit(2);
            }
            String line;
            while ((line = input.readLine()) != null) {
                processList.add(line);
            }
            input.close();
        }catch (IOException e){
            e.printStackTrace();
        }catch (InterruptedException e) {
            log.error("Genome index failed.");
        }

        try {
            FileWriter fw = new FileWriter(refGenomeDir + "gen_idx.log");
            for (String line : processList) {
                fw.write(line);
            }
            fw.close();
        }catch (IOException e){
            e.printStackTrace();
        }
    }

    /**
     * 1-pass mapping with STAR command "STAR --genomeDir $genomeDir --readFilesIn mate1.fq mate2.fq --runThreadN <n>"
     * this command generate Aligned.out.sam files which record the first-time's alignment result and SJ.out.tab file
     * @param refGenomeDir  directory of the index reference genome file
     * @param fastqFile  absolute file path of fastq files
     * @param execThread number of working thread
     */
    public static void onePassMapping(String refGenomeDir, String fastqFile, int execThread, Logger log) {

        String cmd = "STAR --genomeDir " + refGenomeDir + " --readFilesIn " + fastqFile + " --runThreadN " + execThread;
        log.debug("1-pass mapping procedure...");

        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                System.exit(2);
            }
        }catch (IOException e){
            e.printStackTrace();
        }catch (InterruptedException e) {
            log.error("1-pass mapping failed.");
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
                     " --sjdbFileChrStartEnd " + onePassSJOutput + " --sjdbOverhang " + readLength + " --runThreadN " + execThread + " --limitGenomeGenerateRAM 124544990592";
        log.debug("2-pass mapping procedure...");
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
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
     * align reads to reference genome with command "STAR --genomeDir $genomeDir --readFilesIn mate1.fq mate2.fq --runThreadN <n>"
     * the function of this command is same as 1-pass mapping
     * @param refGenomeDir directory of the index reference genome file
     * @param execThread number of working thread
     */
    public static void alignment(String refGenomeDir, String fastqFile, int execThread, Logger log) {

        String cmd = "STAR --genomeDir " + refGenomeDir + " --readFilesIn " + fastqFile + " --runThreadN " +
                execThread;
        log.debug("output the final alignment result...");
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                log.error("Reads alignment process failed");
                System.exit(2);
            }
        }catch (IOException e){
            e.printStackTrace();
        }catch (InterruptedException e) {
            System.out.println("alignment result output failed.");
        }
    }

    /**
     * adding read group information with picard tool
     * @param picardJarDir install directory which contains picard jar package
     * @param starSamFile the sam file output by STAR
     * @param sortedBamFile picard tool output sorted bam file
     */
    public static void readsGroup(String picardJarDir, String starSamFile, String sortedBamFile) {
        String cmd = "java -jar " + picardJarDir + " AddOrReplaceReadGroups I=" + starSamFile + " O=" + sortedBamFile +
                     " SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample TMP_DIR=./tmp";
        log.debug("grouping alignment reads");

        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                log.error("picard reads grouping process failed.");
                System.exit(2);
            }
        }catch (IOException e){
            e.printStackTrace();
        }catch (InterruptedException e) {
            System.out.println("picard reads grouping process failed.");
        }
    }

    /**
     * drop duplicates and indexing reads
     * @param picardJarDir install directory which contains picard jar package
     * @param sortedBamFile picard tool output sorted bam file
     * @param deduplicatedBamFile picard tool output deduplicated bam file
     */
    public static void dropDuplicateReads(String picardJarDir, String sortedBamFile, String deduplicatedBamFile) {
        String cmd = "java -jar " + picardJarDir + " MarkDuplicates I=" + sortedBamFile + " O=" + deduplicatedBamFile +
                     " CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics";
        log.debug("marking duplicated alignment reads");

        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
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
        }catch (IOException e){
            e.printStackTrace();
        }catch (InterruptedException e) {
            System.out.println("picard reads deduplicate process failed.");
        }
    }

    /**
     * create sequence genome with picard tool for further gatk reads trimming
     * @param picardJarDir install directory which contains picard jar package
     * @param refGenomeFastaPath reference sequence file path
     */
    public static void refGenomeDict(String picardJarDir, String refGenomeFastaPath) {
        File referenceFile = new File(refGenomeFastaPath);
        File refSequenceDir = referenceFile.getParentFile();
        int dot = referenceFile.getName().lastIndexOf('.');
        File refSequenceDict = new File(refSequenceDir.getAbsolutePath(), referenceFile.getName().substring(0, dot) + ".dict");
        if (refSequenceDict.exists()) {
            return;
        }

        String cmd = "java -jar " + picardJarDir + " CreateSequenceDictionary R=" + refGenomeFastaPath + " O=" + refSequenceDict.getAbsolutePath();
        log.debug("create reference sequence dictionary");
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                log.error("create sequence dictionary failed");
                System.exit(2);
            }
        }catch (IOException e){
            e.printStackTrace();
        }catch (InterruptedException e) {
            System.out.println("create sequence dictionary failed");
        }
    }

    /**
     * create fai file for further gatk trimming
     * @param refGenomeFastaFile reference sequence file path
     */
    public static void createFastaiFile(String refGenomeFastaFile) {
        File idxFile = new File(refGenomeFastaFile + ".fai");
        if (idxFile.exists()) {
            return;
        }
        String cmd = "samtools faidx " + refGenomeFastaFile;
        log.debug("generate fai file");
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                log.error("create sequence idx failed");
                System.exit(2);
            }
        }catch (IOException e){
            e.printStackTrace();
        }catch (InterruptedException e) {
            System.out.println("create sequence idx failed");
        }
    }

    /**
     * get sequence length from fastq data
     * @param fastqFile fastq file path
     * @return length of reads in fastq file
     */
    public static int getFastqReadLength(String fastqFile) {
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
            System.out.println("fastqFile : [" + fastqFile + "] not found");
            System.exit(2);
        } catch (IOException e) {
            System.out.println("IO error");
        }

        return length;
    }

    /**
     * SNP site reads count by gatk ASEReadCounter tool
     * @param refGenomeFilePath reference genome file
     * @param bamFile gatk alignment bam file
     * @param gatk gatk executive file
     * @param gatkVcfFile gatk output VCF file
     */
    public static void readsCount(String refGenomeFilePath, String bamFile, String gatk, String gatkVcfFile) {
        SnpReadCount src = new SnpReadCount(refGenomeFilePath, bamFile, gatk, gatkVcfFile, log);
        src.snpReads();
    }

    /**
     * merge the output vcf files in one for further heterozygote site analysis
     * @param outputDir directory stores vcf calling result
     * @param gatkLocaljar install directory which contains gatk jar package
     */
    public static void mergeVcfFiles(String outputDir, String gatkLocaljar) {
        File outputDirectory = new File(outputDir);
        File[] snpCallingRes = outputDirectory.listFiles();
        File mergeVcfFile = new File(outputDir, "result.vcf");
        String cmd = gatkLocaljar + " MergeVcfs ";
        StringBuilder sb = new StringBuilder();

        if (snpCallingRes != null) {
            for (File f : snpCallingRes) {
                if (f.getName().endsWith("Select.vcf")) {
                    sb.append(" -I=");
                    sb.append(f.getAbsolutePath());
                }
            }
            cmd = cmd + sb.toString();
            cmd += " -O=" + mergeVcfFile.getAbsolutePath();
            System.out.println(cmd);

            try {
                Process p = Runtime.getRuntime().exec(cmd);
                int exitVal = p.waitFor();
                if (exitVal == 0) {
                    for (File f : snpCallingRes) {
                        if (f.getName().endsWith(".vcf") | f.getName().endsWith(".idx") | f.getName().endsWith(".bai")){
                            p = Runtime.getRuntime().exec("rm -f " + f.getAbsolutePath());
                            exitVal = p.waitFor();
                        }
                    }
                }
            } catch (IOException | InterruptedException ie) {
                ie.printStackTrace();
            }
        }
    }
}
