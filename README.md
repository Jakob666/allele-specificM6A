# Renlab.m6A_allele 1.0

## HARDWARE/SOFTWARE REQUIREMENTS
* Java 1.8+
* Windows / Linux / Mac OS

## INSTALLATION
* clone the repo,
```
git clone https://github.com/Jakob666/allele-specificM6A.git
```
* target JAR package
```
cd ./allele-specificM6A
```
make sure the directory contains `renlabm6a_allele.jar`.

## USAGE
### Overview
#### Tools Introduction
renlabm6a_allele.jar provides the following tools:

Tool | Function
---|---
AseGeneDetection|detect allele-specific expression (ASE) genes (one sample test)
AsmPeakDetection|detect allele-specific modification (ASM) m6A signals  (one sample test)
SampleSpecificASE|detect sample-specific ASE genes  (two sample test)
SampleSpecificASM|detect sample-specific ASM m6A signals  (two sample test)
AseSeqSimulator|generate simulation data for test

#### parameters description
1. **AseGeneDetection**
   - `-vcf/--vcf_file` : required, VCF format file generate by RNA-seq or MeRIP-seq data SNP calling process
   - `-wes/--wes_vcf_file` : optional, VCF format file generate by WES data SNP calling process
   - `-g/--gtf` : required, GTF annotation file
   - `-db/--dbsnp` : optional, big scale SNV annotation data set, like dbsnp, 1000Genome etc. (the file format is described below)
   - `-o/--output` : optional, ASE gene test output file, default `./aseGene.txt`
   - `-df/--degree_of_freedom` : optional, degree of freedom of inverse-Chi-square distribution, default 10
   - `-rc/--reads_coverage` : optional, reads coverage threshold using for filter RNA-seq or MeRIP-seq data SNVs in VCF file (aim for reducing FDR), default 10
   - `-bc/--bkg_coverage` :optional, reads coverage threshold using for filter WES data SNVs in VCF file (aim for reducing FDR), default 30
   - `-s/--sampling`: optional, MCMC sampling times, larger than 500, default 10000
   - `-b/--burn` : optional, MCMC burn-in times, more than 100 and less than sampling times. Default 2000
   - `-t/--thread` : optional, thread number for running test. Default 2
   - `-h/--help` : help message of AseGeneDetection

2. **AsmPeakDetection**
   - `-bed/--peak_bed_file` : required, peak calling output result in BED format
   - `-vcf/--vcf_file` : required, VCF format file generate by RNA-seq or MeRIP-seq data SNP calling process
   - `-wes/--wes_vcf_file` : optional, VCF format file generate by WES data SNP calling process
   - `-g/--gtf` : required, GTF annotation file
   - `-db/--dbsnp` : optional, big scale SNV annotation data set, like dbsnp, 1000Genome etc. (the file format is described below)
   - `-o/--output` : optional, ASM m6A signal test output file, default `./asmPeak.txt`
   - `-df/--degree_of_freedom` : optional, degree of freedom of inverse-Chi-square distribution, default 10
   - `-rc/--reads_coverage` : optional, reads coverage threshold using for filter RNA-seq or MeRIP-seq data SNVs in VCF file (aim for reducing FDR), default 10
   - `-bc/--bkg_coverage` :optional, reads coverage threshold using for filter WES data SNVs in VCF file (aim for reducing FDR), default 30
   - `-s/--sampling`: optional, MCMC sampling times, larger than 500, default 10000
   - `-b/--burn` : optional, MCMC burn-in times, more than 100 and less than sampling times. Default 2000
   - `-t/--thread` : optional, thread number for running test. Default 2
   - `-h/--help` : help message of AsmPeakDetection

3. **SampleSpecificASE**
   - `-s1_vcf/--sample1_vcf_file` : required, VCF format file generate by sample1 RNA-seq or MeRIP-seq data SNP calling process
   - `-s2_vcf/--sample2_vcf_file` : required, VCF format file generate by sample2 RNA-seq or MeRIP-seq data SNP calling process
   - `-s1_wes/--sample1_wes_vcf_file` : optional, VCF format file generate by sample1 WES data SNP calling process
   - `-s2_wes/--sample2_wes_vcf_file` : optional, VCF format file generate by sample2 WES data SNP calling process
   - `-g/--gtf` : required, GTF annotation file
   - `-db/--dbsnp` : optional, big scale SNV annotation data set, like dbsnp, 1000Genome etc. (the file format is described below)
   - `-o/--output` : optional, ASE gene test output file, default `./sampleSpecificASE.txt`
   - `-df/--degree_of_freedom` : optional, degree of freedom of inverse-Chi-square distribution, default 10
   - `-rc/--reads_coverage` : optional, reads coverage threshold using for filter RNA-seq or MeRIP-seq data SNVs in VCF file (aim for reducing FDR), default 10
   - `-bc/--bkg_coverage` :optional, reads coverage threshold using for filter WES data SNVs in VCF file (aim for reducing FDR), default 30
   - `-s/--sampling`: optional, MCMC sampling times, larger than 500, default 10000
   - `-b/--burn` : optional, MCMC burn-in times, more than 100 and less than sampling times. Default 2000
   - `-t/--thread` : optional, thread number for running test. Default 2
   - `-h/--help` : help message of SampleSpecificASE
4. **SampleSpecificASM**
   - `-s1_bed/--sample1_bed_file` : required, sample1 peak calling output result in BED format
   - `-s2_bed/--sample2_bed_file` : required, sample2 peak calling output result in BED format
   - `-s1_vcf/--sample1_vcf_file` : required, VCF format file generate by sample1 RNA-seq or MeRIP-seq data SNP calling process
   - `-s2_vcf/--sample2_vcf_file` : required, VCF format file generate by sample2 RNA-seq or MeRIP-seq data SNP calling process
   - `-s1_wes/--sample1_wes_vcf_file` : optional, VCF format file generate by sample1 WES data SNP calling process
   - `-s2_wes/--sample2_wes_vcf_file` : optional, VCF format file generate by sample2 WES data SNP calling process
   - `-g/--gtf` : required, GTF annotation file
   - `-db/--dbsnp` : optional, big scale SNV annotation data set, like dbsnp, 1000Genome etc. (the file format is described below)
   - `-o/--output` : optional, ASE gene test output file, default `./sampleSpecificASM.txt`
   - `-df/--degree_of_freedom` : optional, degree of freedom of inverse-Chi-square distribution, default 10
   - `-rc/--reads_coverage` : optional, reads coverage threshold using for filter RNA-seq or MeRIP-seq data SNVs in VCF file (aim for reducing FDR), default 10
   - `-bc/--bkg_coverage` :optional, reads coverage threshold using for filter WES data SNVs in VCF file (aim for reducing FDR), default 30
   - `-s/--sampling`: optional, MCMC sampling times, larger than 500, default 10000
   - `-b/--burn` : optional, MCMC burn-in times, more than 100 and less than sampling times. Default 2000
   - `-t/--thread` : optional, thread number for running test. Default 2
   - `-h/--help` : help message of SampleSpecificASM

5. **AseSeqSimulator**
   - `-g/--gtf` : required, GTF annotation file
   - `-t/--twobit` : required, UCSC 2bit format file
   - `-v/--vcf_file` : optional, vcf file path used for generate SNP, default null
   - `-db/--dbsnp` : optional, big scale SNV annotation data set, like dbsnp, 1000Genome etc. (the file format is described below)
   - `-o/--output` : optional, simulate data output directory, default `./AseSeqRead`
   - `-ls/--library_size` : optional, cDNA library size, default 1000000
   - `-rl/--read_length` : optional, sequencing reads length, default 75
   - `-pl/--peak_length` : optional, m6A signal peak length, default 250
   - `-pn/--peakNum` : optional, maximum m6A signal number on a gene, default -1 (no restrict two this parameter)
   - `-min_mut/--minimum_mutation` : optional, minimum mutation sites number on fragment under a m6A signal, default 0
   - `-max_mut/--maximum_mutation` : optional, maximum mutation sites number on fragment under a m6A signal, default 3
   - `-min_cover/--minimum_coverage` : optional, minimum reads coverage when generate RNA-seq data, default 10
   - `-max_cover/--maximum_coverage` : optional, maximum reads coverage when generate RNA-seq data, default 70
   - `-al/--ase_infimum` : optional, ASE ratio infimum, default 0.55
   - `-ah/--ase_supremum` : optional, ASE ratio supremum, default 0.85
   - `-fm/--fragment_mean` : optional, mean length of fragments, default 250
   - `-ft/--fragment_theta` : optional, standard deviation of fragment length, default 25
   - `-gp/--gene_proportion` : optional, the proportion of genes selected in the total number of genes on a particular chromosome, default 0.2
   - `-mp/--mutate_proportion` : optional, the proportion of mutated genes in the total number of selected genes on a particular chromosome, default 0.4
   - `-ol/--overlap` : optional, if the random select gene overlapped with each other, default false
   - `-se/--single_end` : optional, if true, single-end reads, otherwise pair-end. Default true
   - `-pe/--pcr_error` : optional, probability of PCR sequencing error, default 0.005
   - `-h/--help` : help message of AseSeqSimulator

### 1. Allele-specific expression (ASE) gene detection (one sample test)
**data dependency**:
1. VCF format file generate by SNV calling process of `RNA-seq data` or `MeRIP-seq INPUT data` (required, the format of the file is described below)
2. VCF format file generate by SNV calling process of `WES data`(optional)
3. GTF file (required)
4. Big scale SNV annotation data set, like dbsnp, 1000Genome etc. (optional, the format of the file is described below). If `WES data` is supported, this parameter will be **ignore**.

**examples**:\
suppose here exists files below:
1. human genome GTF file `/path/to/Homo_sapiens.GRCh38.93.chr.gtf`
2. VCF format file generate by RNA data `/path/to/rna_filtered.vcf`
3. VCF format file generate by DNA data `/path/to/wes_filtered.vcf`
4. large-scale SNV data set, dbSNP, in VCF format `/path/to/dbsnp.vcf`

* detect ASE gene only by using VCF data generate by `RNA-seq` or `MeRIP-seq INPUT`
```
# command
java -jar ./renlabm6a_allele.jar AseGeneDetection 
     -g /path/to/Homo_sapiens.GRCh38.93.chr.gtf 
     -vcf /path/to/rna_filtered.vcf 
     -o /path/to/output_file 
     -t 6
```
* detect ASE gene by using VCF data generate by `RNA-seq` or `MeRIP-seq INPUT` and `WES`
```
# command
java -jar ./renlabm6a_allele.jar AseGeneDetection 
     -g /path/to/Homo_sapiens.GRCh38.93.chr.gtf 
     -vcf /path/to/rna_filtered.vcf
     -wes /path/to/wes_filtered.vcf
     -o /path/to/output_file 
     -t 6
```
* detect ASE gene by using VCF data generate by `RNA-seq` or `MeRIP-seq INPUT` and large-scale SNV data set
```
# command
java -jar ./renlabm6a_allele.jar AseGeneDetection 
     -g /path/to/Homo_sapiens.GRCh38.93.chr.gtf 
     -vcf /path/to/rna_filtered.vcf
     -db /path/to/dbsnp.vcf
     -o /path/to/output_file 
     -t 6
```


### 2. Allele-specific modification (ASM) m6A signal detection (one sample test)
**data dependency**:
1. GTF format file
2. VCF format file generate by SNV calling process of `RNA-seq data` or `MeRIP-seq INPUT data` (required, the format of the file is described below)
3. BED format peak calling result generate by `MeRIP-seq data` (required, the format of the file is described below)
4. VCF format file generate by SNV calling process of `WES data`(optional)
5. large-scale SNV annotation data set, like dbsnp, 1000Genome etc. (optional, the format of the file is described below). If `WES data` is supported, this parameter will be **ignore**.

**examples**:\
suppose here exists files below:
1. human genome GTF file `/path/to/Homo_sapiens.GRCh38.93.chr.gtf`
2. VCF format file generate by RNA data `/path/to/rna_filtered.vcf`
3. BED format file generate by peak calling process `/path/to/peak.bed`
4. VCF format file generate by DNA data `/path/to/wes_filtered.vcf`
5. large-scale SNV data set, dbSNP, in VCF format `/path/to/dbsnp.vcf`

* detect ASM m6A signal only by using VCF data generate by `RNA-seq` or `MeRIP-seq INPUT`
```
# command
java -jar ./renlabm6a_allele.jar AsmPeakDetection 
     -g /path/to/Homo_sapiens.GRCh38.93.chr.gtf 
     -bed /path/to/peak.bed 
     -vcf /path/to/rna_filtered.vcf 
     -o /path/to/output_file 
     -t 6
```
* detect ASM m6A signal by using VCF data generate by `RNA-seq` or `MeRIP-seq INPUT` and `WES`
```
# command
java -jar ./renlabm6a_allele.jar AsmPeakDetection 
     -g /path/to/Homo_sapiens.GRCh38.93.chr.gtf 
     -bed /path/to/peak.bed 
     -vcf /path/to/rna_filtered.vcf
     -wes /path/to/wes_filtered.vcf
     -o /path/to/output_file
     -t 6
```
* detect ASM m6A signal by using VCF data generate by `RNA-seq` or `MeRIP-seq INPUT` and large-scale SNV data set
```
# command
java -jar ./renlabm6a_allele.jar AsmPeakDetection 
     -g /path/to/Homo_sapiens.GRCh38.93.chr.gtf 
     -bed /path/to/peak.bed 
     -vcf /path/to/rna_filtered.vcf
     -db /path/to/dbsnp.vcf
     -o /path/to/output_file
     -t 6
```

### 3. Sample-specific ASE gene detection (two sample test)
**data dependency**:
1. paired sample VCF format files generate by SNV calling process of `RNA-seq data` or `MeRIP-seq INPUT data` (required, the format of the file is described below)
2. paired sample VCF format files generate by SNV calling process of `WES data`(optional)
3. GTF file (required)
4. Big scale SNV annotation data set, like dbsnp, 1000Genome etc. (optional, the format of the file is described below). If `WES data` is supported, this parameter will be **ignore**.

**examples**:\
suppose here exists files below:
1. human genome GTF file `/path/to/Homo_sapiens.GRCh38.93.chr.gtf`
2. paired sample VCF format files generate by RNA data `/path/to/sample1_rna_filtered.vcf` & `/path/to/sample2_rna_filtered.vcf`
3. paired sample VCF format files generate by WES data `/path/to/sample1_wes_filtered.vcf` & `/path/to/sample2_wes_filtered.vcf`
4. large-scale SNV data set, dbSNP, in VCF format `/path/to/dbsnp.vcf`

* detect sample-specific ASE gene only by using VCF data generate by `RNA-seq` or `MeRIP-seq INPUT`
```
# command
java -jar ./renlabm6a_allele.jar SampleSpecificASE 
     -g /path/to/Homo_sapiens.GRCh38.93.chr.gtf 
     -s1_vcf /path/to/sample1_rna_filtered.vcf 
     -s2_vcf /path/to/sample2_rna_filtered.vcf
     -o /path/to/output_file 
     -t 6
```
* detect sample-specific ASE gene by using VCF data generate by `RNA-seq` or `MeRIP-seq INPUT` and `WES`
```
# command
java -jar ./renlabm6a_allele.jar SampleSpecificASE 
     -g /path/to/Homo_sapiens.GRCh38.93.chr.gtf 
     -s1_vcf /path/to/sample1_rna_filtered.vcf 
     -s2_vcf /path/to/sample2_rna_filtered.vcf
     -s1_wes /path/to/sample1_wes_filtered.vcf
     -s2_wes /path/to/sample2_wes_filtered.vcf
     -o /path/to/output_file 
     -t 6
```
* detect sample-specific ASE gene by using VCF data generate by `RNA-seq` or `MeRIP-seq INPUT` and large-scale SNV data set
```
# command
java -jar ./renlabm6a_allele.jar SampleSpecificASE 
     -g /path/to/Homo_sapiens.GRCh38.93.chr.gtf 
     -s1_vcf /path/to/sample1_rna_filtered.vcf 
     -s2_vcf /path/to/sample2_rna_filtered.vcf
     -db /path/to/dbsnp.vcf
     -o /path/to/output_file 
     -t 6
```

### 4. Sample-specific ASM m6A signal detection (two sample test)
**data dependency**:
1. GTF format file
2. paired sample VCF format files generate by SNV calling process of `RNA-seq data` or `MeRIP-seq INPUT data` (required, the format of the file is described below)
3. paired sample BED format peak calling results generate by `MeRIP-seq data` (required, the format of the file is described below)
4. paired sample VCF format files generate by SNV calling process of `WES data`(optional)
5. large-scale SNV annotation data set, like dbsnp, 1000Genome etc. (optional, the format of the file is described below). If `WES data` is supported, this parameter will be **ignore**.

**examples**:\
suppose here exists files below:
1. human genome GTF file `/path/to/Homo_sapiens.GRCh38.93.chr.gtf`
2. paired sample VCF format files generate by RNA data `/path/to/sample1_rna_filtered.vcf` & `/path/to/sample2_rna_filtered.vcf`
3. paired sample BED format files generate by peak calling process `/path/to/sample1_peak.bed` & `/path/to/sample2_peak.bed`
4. paired sample VCF format files generate by WES data `/path/to/sample1_wes_filtered.vcf` & `/path/to/sample2_wes_filtered.vcf`
5. large-scale SNV data set, dbSNP, in VCF format `/path/to/dbsnp.vcf`

* detect sample-specific ASM m6A signal only by using VCF data generate by `RNA-seq` or `MeRIP-seq INPUT`
```
# command
java -jar ./renlabm6a_allele.jar SampleSpecificASM 
     -g /path/to/Homo_sapiens.GRCh38.93.chr.gtf 
     -s1_bed /path/to/sample1_peak.bed
     -s2_bed /path/to/sample2_peak.bed 
     -s1_vcf /path/to/sample1_rna_filtered.vcf 
     -s2_vcf /path/to/sample2_rna_filtered.vcf 
     -o /path/to/output_file 
     -t 6
```
* detect sample-specific ASM m6A signal by using VCF data generate by `RNA-seq` or `MeRIP-seq INPUT` and `WES`
```
# command
java -jar ./renlabm6a_allele.jar SampleSpecificASM 
     -g /path/to/Homo_sapiens.GRCh38.93.chr.gtf 
     -s1_bed /path/to/sample1_peak.bed
     -s2_bed /path/to/sample2_peak.bed 
     -s1_vcf /path/to/sample1_rna_filtered.vcf 
     -s2_vcf /path/to/sample2_rna_filtered.vcf
     -s1_wes /path/to/sample1_wes_filtered.vcf
     -s2_wes /path/to/sample2_wes_filtered.vcf
     -o /path/to/output_file
     -t 6
```
* detect sample-specific ASM m6A signal by using VCF data generate by `RNA-seq` or `MeRIP-seq INPUT` and large-scale SNV data set
```
# command
java -jar ./renlabm6a_allele.jar SampleSpecificASM 
     -g /path/to/Homo_sapiens.GRCh38.93.chr.gtf 
     -s1_bed /path/to/sample1_peak.bed
     -s2_bed /path/to/sample2_peak.bed 
     -s1_vcf /path/to/sample1_rna_filtered.vcf 
     -s2_vcf /path/to/sample2_rna_filtered.vcf
     -db /path/to/dbsnp.vcf
     -o /path/to/output_file
     -t 6
```

### 5. Simulated data generation
generate simulated data for test,
**data dependency**:
1. GTF format file (required)
2. two bit file (required)

**examples**:\
suppose here exists files below:
1. human genome GTF file `/path/to/Homo_sapiens.GRCh38.93.chr.gtf`
2. human genome 2bit file `/path/to/hg38.2bit`

simulated data parameters:
1. sequencing read length: `read_length`, integer
2. simulated m6A signal length: `peak length`, integer
3. library size: `library_size`, integer
4. min/max_cover: minimum and maximum gene reads coverage of simulated MeRIP-seq data, integer, denote as `coverage_infimum` and `coverage_supremum` respectively
5. dep: sequencing depth of simulated WES data, `depth`, integer
6. al/ah: infimum and supremum of major allele frequency, float, `al ≤ ah < 1`
7. gp: the proportion of mutated genes in the total number of selected genes, float, `gp < 1`
* generate MeRIP-seq simulated data
```
# command
java -jar ./renlabm6a_allele.jar AseSeqSimulator
     -g /path/to/Homo_sapiens.GRCh38.93.chr.gtf
     -t /path/to/hg38.2bit
     -o /path/to/output_dir
     -rl read_length
     -pl peak_length
     -ls library_size
     -min_cover coverage_infimum
     -max_cover coverage_supremum
     -al major_allele_frequency_infimum
     -ah major_allele_frequency_supremum
     -gp mutated_gene_proportion
```
* generate MeRIP-seq and WES simulated data
```
# command
java -jar ./renlabm6a_allele.jar AseSeqSimulator
     -g /path/to/Homo_sapiens.GRCh38.93.chr.gtf
     -t /path/to/hg38.2bit
     -o /path/to/output_dir
     -rl read_length
     -pl peak_length
     -ls library_size
     -dep WES_sequencing_depth
     -min_cover coverage_infimum
     -max_cover coverage_supremum
     -al major_allele_frequency_infimum
     -ah major_allele_frequency_supremum
     -gp mutated_gene_proportion
```

## FORMAT DECLARATION
### 1. VCF generate by SNV calling of RNA and DNA sequencing data
At least 7 columns,
* \#CHROM: chromosome number, `1,2,3,...X,Y,MT`
* POS: mutation position
* ID: mutation ID, default `.`
* REF: reference nucleotide
* ALT: alternative nucleotide
* QUAL: quality score
* FILTER: `PASS` if SNV sites were filtered, default `.`
* INFO: **if** VCF file contains no `FORMAT` and `Sample` columns or there's no `AD` item in `FORMAT` columns, there **must** exists `DP4` item in `INFO` column, like `DP4=9,8,3,6` denotes as the number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases, respectively.  
* FORMAT (optional): make sure if there's no `DP4` item in `INFO` column, this data field **must** contains `AD` item, which represents allelic depths for the ref and alt alleles in the order listed.\
If the records in VCF file exist `DP4` and `AD` items at once, data in `AD` item will be given priority.
* Sample (optional): the corresponding data in `FORMAT` field.

The two examples bellow are satisfied to the above requirements,
* example 1
> \#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample\
> 1	3025531	.	T	A	64.28	PASS	AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=60.00;QD=32.14;SOR=0.693	GT:**AD**:DP:GQ:PL	1/1:0,2:2:6:76,6,0\
> 1	3037125	.	A	C	68.28	PASS	AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=60.00;QD=34.14;SOR=0.693	GT:**AD**:DP:GQ:PL	1/1:0,2:2:6:80,6,0\
> 1	5170624	.	A	G	434.6	SnpCluster	AC=1;AF=0.500;AN=2;BaseQRankSum=2.479;DP=17;ExcessHet=3.0103;FS=5.315;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=25.56;ReadPosRankSum=-1.640;SOR=0.662	GT:**AD**:DP:GQ:PL	0/1:5,12:17:99:442,0,142\
> 1	85864585	.	T	A,C	771.02	PASS	AC=1,1;AF=0.500,0.500;AN=2;DP=20;ExcessHet=3.0103;FS=0.000;MLEAC=1,1;MLEAF=0.500,0.500;MQ=60.00;QD=31.86;SOR=1.022	GT:**AD**:DP:GQ:PL	1/2:0,5,14:19:99:788,579,564,209,0,167

* example 2
> \#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample\
> 1	959078	.	G	T	12.4875	PASS	DP=35;VDB=0.198714;SGB=-0.662043;RPB=0.471274;MQB=1;MQSB=1;BQB=1;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;**DP4**=15,11,3,6;MQ=20	GT:PL:DP:SP	0/1:46,0,139:35:6\
> 1	1426307	.	G	C	106	PASS	DP=83;VDB=0.1007;SGB=-0.693145;RPB=0.908855;MQB=1;MQSB=1;BQB=1;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;**DP4**=15,27,21,20;MQ=20	GT:PL:DP:SP	0/1:139,0,145:83:7\
> 1	1426367	.	A	G	99	PASS	DP=57;VDB=0.189985;SGB=-0.693021;RPB=0.998851;MQB=1;MQSB=1;BQB=1;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;**DP4**=18,12,14,13;MQ=20	GT:PL:DP:SP	0/1:132,0,134:57:2\
> 1	2430430	.	T	A	36.5982	PASS	DP=51;VDB=0.923739;SGB=-0.688148;RPB=0.994672;MQB=1;MQSB=1;BQB=1;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;**DP4**=15,21,5,10;MQ=20	GT:PL:DP:SP	0/1:70,0,149:51:1\
> 1	13892877	.	C	A	20.4138	PASS	DP=37;VDB=0.558487;SGB=-0.670168;RPB=0.149569;MQB=1;MQSB=1;BQB=1;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;**DP4**=15,12,5,5;MQ=20	GT:PL:DP:SP	0/1:54,0,142:37:0


### 2. Large-scale mutation data set
Must contains chromosome number and mutation position, e.g.
> \#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\
> MT	16463	rs878856034	A	G	.	.	RS=878856034;RSPOS=16463;dbSNPBuildID=147;SSR=0;SAO=0;VP=0x050000020005000002000100;GENEINFO=ND6:4541;WGT=1;VC=SNV;R5;ASP\
> MT	16468	rs879224782	T	C	.	.	RS=879224782;RSPOS=16468;dbSNPBuildID=147;SSR=0;SAO=0;VP=0x050000020005000002000100;GENEINFO=ND6:4541;WGT=1;VC=SNV;R5;ASP\
> MT	16474	rs878872875	G	C	.	.	RS=878872875;RSPOS=16474;dbSNPBuildID=147;SSR=0;SAO=0;VP=0x050000020005000002000100;GENEINFO=ND6:4541;WGT=1;VC=SNV;R5;ASP\
> MT	16482	rs878935154	A	G	.	.	RS=878935154;RSPOS=16482;dbSNPBuildID=147;SSR=0;SAO=0;VP=0x050000020005000002000100;GENEINFO=ND6:4541;WGT=1;VC=SNV;R5;ASP

### 3. BED format file
Contains fields below, more details see [BED format demonstration UCSC](http://genome.ucsc.edu/FAQ/FAQformat#format1)
* \# chr: chromosome number, `1,2,3,...X,Y,MT`
* chromStart: m6A signal start position on chromosome
* chromEnd: m6A signal end position on chromosome
* name: ESMBEL gene ID
* score: significant score(adjusted p.value), generate by peak calling tools, less than `0.05`
* strand: `+` or `-`
* thickStart: The starting position at which the feature is drawn thickly
* thickEnd: The ending position at which the feature is drawn thickly
* itemRgb: An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. 
* blockCount: sub-block number of the m6A signal peak, integer, `≥1`
* blockSizes: block size of each sub-block, separate by `,`
* blockStarts: block start position on chromosome, separate by `,`

> \# chr	chromStart	chromEnd	name	score	strand	thickStart	thickEnd	itemRgb	blockCount	blockSizes	blockStarts\
> 1	9747647	9747845	ENSMUSG00000097893	7.1e-05	+	9747647	9747845	0	1	198,	0\
> 1	16105773	16105923	ENSMUSG00000025921	4.9e-05	+	16105773	16105923	0	1	150,	0\
> 1	33739519	33739819	ENSMUSG00000004768	0.0032	+	33739519	33739819	0	1	300,	0\
> 1	34180162	34180463	ENSMUSG00000026131	0.00022	+	34180162	34180463	0	1	301,	0\
> 1	34306583	34307612	ENSMUSG00000026131	0.00038	+	34306583	34307612	0	2	68,283,	0,746
