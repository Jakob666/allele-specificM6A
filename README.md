# Renlab.m6A_allele 1.0

## HARDWARE/SOFTWARE REQUIREMENTS
* Java 1.8
* 64 bit Linux

## INSTALLATION
* clone the repo,
```
git clone https://github.com/Jakob666/allele-specificM6A.git
```
* target JAR package
```
cd ./allele-specificM6A
```
make sure the directory contains `renlab.m6a_allele-1.0.jar`.

## USAGE
### 1. Allele-specific expression (ASE) gene detection
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
java -cp ./renlab.m6a_allele-1.0.jar HierarchicalBayesianAnalysis.AseGeneDetection 
     -g /path/to/Homo_sapiens.GRCh38.93.chr.gtf 
     -vcf /path/to/rna_filtered.vcf 
     -o output_file_name 
```
* detect ASE gene by using VCF data generate by `RNA-seq` or `MeRIP-seq INPUT` and `WES`
```
# command
java -cp ./renlab.m6a_allele-1.0.jar HierarchicalBayesianAnalysis.AseGeneDetection 
     -g /path/to/Homo_sapiens.GRCh38.93.chr.gtf 
     -vcf /path/to/rna_filtered.vcf
     -wes /path/to/wes_filtered.vcf
     -o output_file_name 
```
* detect ASE gene by using VCF data generate by `RNA-seq` or `MeRIP-seq INPUT` and large-scale SNV data set
```
# command
java -cp ./renlab.m6a_allele-1.0.jar HierarchicalBayesianAnalysis.AseGeneDetection 
     -g /path/to/Homo_sapiens.GRCh38.93.chr.gtf 
     -vcf /path/to/rna_filtered.vcf
     -db /path/to/dbsnp.vcf
     -o output_file_name 
```


### 2. Allele-specific modification (ASM) m6A signal detection 
**data dependency**:
1. VCF format file generate by SNV calling process of `RNA-seq data` or `MeRIP-seq INPUT data` (required, the format of the file is described below)
2. BED format peak calling result generate by `MeRIP-seq data` (required, the format of the file is described below)
2. VCF format file generate by SNV calling process of `WES data`(optional)
4. large-scale SNV annotation data set, like dbsnp, 1000Genome etc. (optional, the format of the file is described below). If `WES data` is supported, this parameter will be **ignore**.

**examples**:\
suppose here exists files below:
2. VCF format file generate by RNA data `/path/to/rna_filtered.vcf`
2. BED format file generate by peak calling process `/path/to/peak.bed`
3. VCF format file generate by DNA data `/path/to/wes_filtered.vcf`
4. large-scale SNV data set, dbSNP, in VCF format `/path/to/dbsnp.vcf`

* detect ASE gene only by using VCF data generate by `RNA-seq` or `MeRIP-seq INPUT`
```
# command
java -cp ./renlab.m6a_allele-1.0.jar HierarchicalBayesianAnalysis.AsmPeakDetection 
     -bed /path/to/peak.bed 
     -vcf /path/to/rna_filtered.vcf 
     -o output_file_name 
```
* detect ASE gene by using VCF data generate by `RNA-seq` or `MeRIP-seq INPUT` and `WES`
```
# command
java -cp ./renlab.m6a_allele-1.0.jar HierarchicalBayesianAnalysis.AsmPeakDetection 
     -bed /path/to/peak.bed 
     -vcf /path/to/rna_filtered.vcf
     -wes /path/to/wes_filtered.vcf
     -o output_file_name 
```
* detect ASE gene by using VCF data generate by `RNA-seq` or `MeRIP-seq INPUT` and large-scale SNV data set
```
# command
java -cp ./renlab.m6a_allele-1.0.jar HierarchicalBayesianAnalysis.AsmPeakDetection 
     -bed /path/to/peak.bed 
     -vcf /path/to/rna_filtered.vcf
     -db /path/to/dbsnp.vcf
     -o output_file_name 
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
* **INFO**: must contains `DP4` information, like `DP4=9,8,3,6` denotes as the number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases, respectively. \
If there is no `DP4` information, use `AD` information (Allelic depths for the ref and alt alleles in the order listed) of each record to generate `DP4` information and add into `INFO` column, e.g. `GT:AD:DP:GQ:PL	0/1:4,4:8:99:155,0,149` equals `DP4=4,0,4,0` 

### 2. Large-scale mutation data set
Must contains chromosome number and mutation position, e.g.
> MT	16463	rs878856034	A	G	.	.	RS=878856034;RSPOS=16463;dbSNPBuildID=147;SSR=0;SAO=0;VP=0x050000020005000002000100;GENEINFO=ND6:4541;WGT=1;VC=SNV;R5;ASP
  MT	16468	rs879224782	T	C	.	.	RS=879224782;RSPOS=16468;dbSNPBuildID=147;SSR=0;SAO=0;VP=0x050000020005000002000100;GENEINFO=ND6:4541;WGT=1;VC=SNV;R5;ASP
  MT	16474	rs878872875	G	C	.	.	RS=878872875;RSPOS=16474;dbSNPBuildID=147;SSR=0;SAO=0;VP=0x050000020005000002000100;GENEINFO=ND6:4541;WGT=1;VC=SNV;R5;ASP
  MT	16482	rs878935154	A	G	.	.	RS=878935154;RSPOS=16482;dbSNPBuildID=147;SSR=0;SAO=0;VP=0x050000020005000002000100;GENEINFO=ND6:4541;WGT=1;VC=SNV;R5;ASP
