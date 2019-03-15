package SNPCalling;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class SNPCaller {

    public static void main(String[] args) {

        String genomeFilePath = args[0];
        String FastqDir = args[1];
        String outputDir = args[2];
        String picardDir = args[3];
        String gatkLocalJar = args[4];
        int execThread = Integer.parseInt(args[5]);

        snpCalling(genomeFilePath, FastqDir, outputDir, picardDir, gatkLocalJar, execThread);
    }

    /**
     * complete snp calling process
     * @param refGenomeFilePath reference sequence file path
     * @param fastqDir directory of the INPUT fastq file
     * @param outputDir output result directory
     * @param picardJarDir path of picard tool executive jar
     * @param gatkJarDir path of gatk tool executive file
     */
    public static void snpCalling(String refGenomeFilePath, String fastqDir, String outputDir, String picardJarDir, String gatkJarDir, int execThread) {
        File refGenomePath = new File(refGenomeFilePath);
        String refGenomeDir = refGenomePath.getParent();
        int readLength = 0;

        File fastqDirectory = new File(fastqDir);
        File outputCellDir = new File(outputDir);
        if (!outputCellDir.exists()) {
            boolean res = outputCellDir.mkdir();
        }

        String outputCellDirPath = outputCellDir.getAbsolutePath();

        File[] fastqFiles = fastqDirectory.listFiles();
        for (File f : fastqFiles) {
            if (f.getName().endsWith(".fastq")) {
                readLength = getFastqReadLength(f.getAbsolutePath()) - 1;
                int dot = f.getName().lastIndexOf('.');
                String sraNum = f.getName().substring(0, dot);

                if (new File(outputCellDir, sraNum + "_RefinedSNP.vcf").exists())
                    continue;

                readsMapping(refGenomeDir, refGenomePath.getAbsolutePath(), f.getAbsolutePath(), readLength, execThread);
                filterMappedReads(refGenomeDir, picardJarDir, outputCellDirPath, sraNum);
                refGenomeDict(picardJarDir, refGenomeFilePath);
                createFastaiFile(refGenomeFilePath);
                readsTrimReassign(gatkJarDir, refGenomeFilePath, outputCellDirPath, sraNum);
                variantCalling(gatkJarDir, refGenomeFilePath, outputCellDirPath, sraNum);
                variantFilter(gatkJarDir, refGenomeFilePath, outputCellDirPath, sraNum);
            }
        }
    }

    /**
     * mapping RNA-Seq reads to reference genome using STAR
     * @param refGenomeDir directory of reference genome file
     * @param execThread number of thread which execute reads mapping process
     * @param fastqFile absolute path of fastq file
     * @param genomeFileName name of reference genome file without path
     */
    public static void readsMapping(String refGenomeDir, String genomeFileName, String fastqFile, int readLength, int execThread) {
        File genomeIdxFile = new File(refGenomeDir, "SAindex");

        if (!genomeIdxFile.exists()) {
            genomeIndex(refGenomeDir, genomeFileName, execThread);
        }
        onePassMapping(refGenomeDir, fastqFile, execThread);
        twoPassMapping(refGenomeDir, genomeFileName, readLength, execThread);
        alignment(refGenomeDir, fastqFile, execThread);
    }

    /**
     * adding read group information, sorting, marking duplicates and indexing. The reads are from SAM file output by readsMapping method
     * @param picardJarDir install directory which contains picard jar package
     * @param outputDir the output result directory which contains STAR alignment output file
     */
    public static void filterMappedReads(String refGenomeDir, String picardJarDir, String outputDir, String sraNum) {
        File starSamFile = new File(refGenomeDir, "Aligned.out.sam");
        File sortedBamFile = new File(outputDir, sraNum + "_sorted.bam");
        File deduplicatedBamFile = new File(outputDir, sraNum + "_deduplicated.bam");
        readsGroup(picardJarDir, starSamFile.getAbsolutePath(), sortedBamFile.getAbsolutePath());
        dropDuplicateReads(picardJarDir, sortedBamFile.getAbsolutePath(), deduplicatedBamFile.getAbsolutePath());
    }

    /**
     *
     * @param gatkLocalJar install directory which contains gatk jar package
     * @param refGenomeFile reference genome file path
     * @param outputDir output result directory
     */
    public static void readsTrimReassign(String gatkLocalJar, String refGenomeFile, String outputDir, String sraNum) {
        File deduplicatedBamFile = new File(outputDir, sraNum + "_deduplicated.bam");
        File trimmedBamFile = new File(outputDir, sraNum + "_trimmed.bam");
        String cmd = gatkLocalJar + " SplitNCigarReads -R " + refGenomeFile + " -I " + deduplicatedBamFile.getAbsolutePath() +
                     " -O " + trimmedBamFile.getAbsolutePath();
        System.out.println(cmd);

        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                System.out.println("gatk reads trim process failed.");
                System.exit(2);
            }

            p = Runtime.getRuntime().exec("rm -f " + deduplicatedBamFile.getAbsolutePath());
            exitVal = p.waitFor();
            if (exitVal != 0) {
                System.out.println("remove previous file failed.");
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
     * @param sraNum sra ID number of this data
     */
    public static void variantCalling(String gatkLocalJar, String genomeFileName, String outputDir, String sraNum) {
        File trimmedBamFile = new File(outputDir, sraNum + "_trimmed.bam");
        File outputVCF = new File(outputDir, sraNum + "RawSNP.vcf");
        String cmd = gatkLocalJar + " HaplotypeCaller -R " + genomeFileName + " -I " +
                     trimmedBamFile.getAbsolutePath() + " -O " + outputVCF.getAbsolutePath();
        System.out.println(cmd);

        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                System.out.println("raw SNP calling failed");
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
     * @param sraNum sra ID number of this data
     */
    public static void variantFilter(String gatkLocalJar, String genomeFileName, String outputDir, String sraNum) {
        File rawVcfFile = new File(outputDir, sraNum + "RawSNP.vcf");
        File refinedVCF = new File(outputDir, sraNum + "RefinedSNP.vcf");
        String cmd = gatkLocalJar + " VariantFiltration -R " + genomeFileName + " -V " + rawVcfFile.getAbsolutePath() +
                     " -window 35 -cluster 3 -O " + refinedVCF.getAbsolutePath();
        System.out.println(cmd);

        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                System.out.println("filter raw vcf failed");
                System.exit(2);
            }

            p = Runtime.getRuntime().exec("rm -f " + rawVcfFile.getAbsolutePath());
            exitVal = p.waitFor();
            if (exitVal != 0) {
                System.out.println("remove previous file failed.");
                System.exit(2);
            }

        } catch (IOException | InterruptedException ie) {
            ie.printStackTrace();
        }
    }

    /**
     * index the reference genome with STAR command "STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles hg19.fa --runThreadN"
     * this STAR command will generate two file named "Genome" and "SAindex" for further 1-pass mapping
     * @param refGenomeDir  directory of the reference genome file
     * @param genomeFilePath    reference genome file Path
     * @param execThread    number of working thread
     */
    public static void genomeIndex(String refGenomeDir, String genomeFilePath, int execThread) {
        String idxCmd = "STAR --runMode genomeGenerate --genomeDir " + refGenomeDir + " --genomeFastaFiles " + genomeFilePath +
                        " --runThreadN " + execThread + " --limitGenomeGenerateRAM=124544990592";

        System.out.println(idxCmd);
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
            System.out.println("Genome index failed.");
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
    public static void onePassMapping(String refGenomeDir, String fastqFile, int execThread) {

        String cmd = "STAR --genomeDir " + refGenomeDir + " --readFilesIn " + fastqFile + " --runThreadN " + execThread;
        System.out.println(cmd);

        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                System.exit(2);
            }
        }catch (IOException e){
            e.printStackTrace();
        }catch (InterruptedException e) {
            System.out.println("Genome index failed.");
        }
    }

    /**
     * 2-pass mapping with STAR command "STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles hg19.fa --sjdbFileChrStartEnd /path/to/1pass/SJ.out.tab --sjdbOverhang 75 --runThreadN <n>"
     * the command result in a new Genome, SA and SAindex file and sjdbList.out.tab file
     * @param refGenomeDir directory of the index reference genome file
     * @param genomeFileName reference genome file name
     * @param execThread number of working thread
     */
    public static void twoPassMapping(String refGenomeDir, String genomeFileName, int readLength, int execThread) {

        File onePassSJOutput = new File(refGenomeDir, "SJ.out.tab");
        String cmd = "STAR --runMode genomeGenerate --genomeDir " + refGenomeDir + " --genomeFastaFiles " + genomeFileName +
                     " --sjdbFileChrStartEnd " + onePassSJOutput + " --sjdbOverhang " + readLength + " --runThreadN " + execThread + " --limitGenomeGenerateRAM 124544990592";
        System.out.println(cmd);
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                System.exit(2);
            }
        }catch (IOException e){
            e.printStackTrace();
        }catch (InterruptedException e) {
            System.out.println("Genome index failed.");
        }
    }

    /**
     * align reads to reference genome with command "STAR --genomeDir $genomeDir --readFilesIn mate1.fq mate2.fq --runThreadN <n>"
     * the function of this command is same as 1-pass mapping
     * @param refGenomeDir directory of the index reference genome file
     * @param execThread number of working thread
     */
    public static void alignment(String refGenomeDir, String fastqFile, int execThread) {

        String cmd = "STAR --genomeDir " + refGenomeDir + " --readFilesIn " + fastqFile + " --runThreadN " +
                execThread;
        System.out.println(cmd);
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                System.out.println("Reads alignment process failed");
                System.exit(2);
            }
        }catch (IOException e){
            e.printStackTrace();
        }catch (InterruptedException e) {
            System.out.println("Genome index failed.");
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
        System.out.println(cmd);

        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                System.out.println("picard reads grouping process failed.");
                System.exit(2);
            }
        }catch (IOException e){
            e.printStackTrace();
        }catch (InterruptedException e) {
            System.out.println("Genome index failed.");
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
        System.out.println(cmd);

        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                System.out.println("picard reads deduplicate process failed.");
                System.exit(2);
            }

            p = Runtime.getRuntime().exec("rm -f " + sortedBamFile);
            exitVal = p.waitFor();
            if (exitVal != 0) {
                System.out.println("remove previous file failed.");
                System.exit(2);
            }
        }catch (IOException e){
            e.printStackTrace();
        }catch (InterruptedException e) {
            System.out.println("Genome index failed.");
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
        System.out.println(cmd);
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                System.out.println("create sequence dictionary failed");
                System.exit(2);
            }
        }catch (IOException e){
            e.printStackTrace();
        }catch (InterruptedException e) {
            System.out.println("Genome index failed.");
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
        System.out.println(cmd);
        try {
            Process p = Runtime.getRuntime().exec(cmd);
            int exitVal = p.waitFor();
            if (exitVal != 0) {
                System.out.println("create sequence idx failed");
                System.exit(2);
            }
        }catch (IOException e){
            e.printStackTrace();
        }catch (InterruptedException e) {
            System.out.println("Genome index failed.");
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
}
