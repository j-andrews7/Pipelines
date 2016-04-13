# Coding Variants Pipeline
**Up to date as of 04/13/2016.**  
This is an imitation of Liv's coding variant calling pipeline, albeit with some improvements to increase sensitivity (hypothetically). This also began as a comparison between the samtools and VarScan variant callers as well, but after the analysis, it seemed the best bet was to simply merge the results from the two callers, as they have fairly high overlap.

This was done on the CHPC cluster, so all of the `export`, `source`, and `module load/remove` statements are to load the various software necessary to run the command(s) that follow. If you're running this locally and the various tools needed are located on your `PATH`, you can ignore these.

> Bash scripts are submitted on the cluster with the `qsub` command. Check the [CHPC wiki](http://mgt2.chpc.wustl.edu/wiki119/index.php/Main_Page) for more info on cluster commands and what software is available. All scripts listed here should be accessible to anyone in the Payton Lab, i.e., **you should be able to access everything in my scratch folder and run the scripts from there if you so choose.**

All necessary scripts should be here: **N:\Bioinformatics\Jareds_Code**  
They are also in `/scratch/jandrews/bin/` or `/scratch/jandrews/Bash_Scripts/` on the cluster as well as stored on my local PC and external hard drive.  

An _actual_ workflow (Luigi, Snakemake, etc) could easily be made for this with a bit of time, maybe I'll get around to it at some point.

**IMPORTANT**: Be sure to sort and index BAMs before beginning this. There are obviously named bash scripts in the above code folders to do so. All Python scripts here use **Python 3**, so be sure you have the appropriate version installed/set.

**Software Requirements:**
- [BEDOPS](http://bedops.readthedocs.org/en/latest/index.html)
- [Samtools, BCFtools, and HTSlib](http://www.htslib.org/)  
  - These should be available on the CHPC cluster.
- [Python3](https://www.python.org/downloads/)
  - Use an [anaconda environment](http://mgt2.chpc.wustl.edu/wiki119/index.php/Python#Anaconda_Python) if on the CHPC cluster. (also useful for running various versions of python locally).  
- [bedtools](http://bedtools.readthedocs.org/en/latest/)
  - Also available on the CHPC cluster.
- [VarScan2](http://dkoboldt.github.io/varscan/)
- [GATK](https://www.broadinstitute.org/gatk/download/)
- [Picard Tools](http://broadinstitute.github.io/picard/)


---

##### 1A.) Call variants with VarScan.
samtools: mpileup piped to varscan to call variants with filters for read depth (5) and quality (15).

**Bash script (var_call_varscan.sh):**
```Bash
#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N Var_Call_VarScan
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=24gb

module load samtools-1.2
module load bcftools-1.2

for file in /scratch/jandrews/Data/RNA_Seq/ALIGNED_BAMs/Batch28/*.bam; do

	samtools mpileup -t DP -f /scratch/jandrews/Ref/hg19.fa $file | java -Xmx15g -jar /scratch/jandrews/bin/VarScan.v2.3.9.jar mpileup2cns --min-coverage 5 --min-avg-qual 15 --variants 1 --output-vcf 1 >"$file".vs.vcf  &
	
done
wait
module remove samtools-1.2
module remove bcftools-1.2
```

##### 1B i.) Call variants with bcftools (samtools).
**Bash script (var_call_bcf.sh):**
```Bash
#!/bin/sh
# give the job a name to help keep track of running jobs (optional)
#PBS -N Var_Call_RNA
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=8:00:00,vmem=16gb

module load samtools-1.2
module load bcftools-1.2

for file in /scratch/jandrews/Data/RNA_Seq/ALIGNED_BAMs/Batch9/*.bam; do

	samtools mpileup -u -t DP -f /scratch/jandrews/Ref/hg19.fa $file | bcftools call -cv -O v - | /scratch/jandrews/bin/vcfutils.pl varFilter -D100 > $file.vcf &
	
done
wait
module remove samtools-1.2
module remove bcftools-1.2
```


##### 1B ii.) 
Filter the BCFtools VCFs for the same quality and read-depth that VarScan used.

```Bash
module load bcftools-1.2
for f in *.vcf; do
    base=${f##*/}
    bcftools filter -i 'DP>=5 & QUAL>=15' "$f" > ${base%.*}.filtered.vcf
done
module remove bcftools-1.2
```

At this point, I also did a comparison between the SNP arrays with these callers, which is explained [below](#comparison) after the rest of the pipeline. 


##### 2.) Clean VCFs.
```Bash
for file in *.vcf; do
	base=${file%%_*} ;
	(sed '/_g/d' "$file" | sed '/chrM/d' | sed '/chrY/d') > "$base"_RNAseq_VS_clean.vcf ;
done
```

##### 3.) Sort all VCFs.
```Bash
for f in *.vcf; do
   	base=${f%%_*} ;
    vcf-sort "$f" > ${base}_RNAseq_VS_sorted.vcf ;
done
```

##### 4.) Zip and index each VCF.
```Bash
for f in *.vcf; do
	bgzip -c "$f" > "$f".gz;
done

for f in *.gz; do
	tabix -p vcf "$f" ;
done
```

##### 5.) Merge VCFs.
First, those from VarScan with each other. Then those from BCFTools with each other. `-m none` means multiallelic records will be split to separate lines. Doing so is rather important, as if two samples have different variant alleles at the same position, only one is reported as having the variant if multiallelic records are allowed. Alternatively, setting `-m both` should create a multiallelic record, which may be wanted at times. No idea what the default is, BCFtools docs don't mention.

**For VarScan:**
```Bash
module load bcftools-1.2
bcftools merge -O v -m none --force-samples -i ADP:sum *.gz > merge.vcf
module remove bcftools-1.2
```

**For BCFTools:**
```Bash
module load bcftools-1.2
bcftools merge -O v -m none -i DP:sum *.gz > merge.vcf
module remove bcftools-1.2
```

**TO DO - ROLL IN CODING VARIANTS FROM CHIP-SEQ DATA AFTER FUNSEQ ANNOTATION AROUND HERE.**
General idea:
- Parse FunSeq output to get coding variants for each sample.
- Backwards intersect these with the ChIP-seq VCF for the given sample to get the genotype and other info that FunSeq utterly destroyed for each variant.
- Take all the resulting VCFs and merge them, somewhat similar to the above commands.
- Include this file in the combining steps below.

##### 6.) Fix headers.
Text editor style because I was too lazy to write something.

##### 7.) Combine the BCFTools and VarScan files.
This step was an **enormous** hassle to figure out.

**i. Create a sequence dict for reference genome**  
This has to be done for GATK/Picard tools to work properly. Downloading the [bundle](ftp://ftp.broadinstitute.org/bundle/2.8/) for your genome of interest might be easier than doing this really.  
**Bash script (create_ref_dict.sh):**
```Bash
#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N CREATE_REF_DICT
#PBS -m e
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=36gb
#PBS -q old

module load java
java -jar /scratch/jandrews/bin/picard-tools-2.2.1/picard.jar CreateSequenceDictionary R= /scratch/jandrews/Ref/hg19.fa O= /scratch/jandrews/Ref/hg19.dict
module remove java
```

**ii. Sort VCFs to same order as reference genome.**  
I still don't get why the tools can't figure this out on their own, but again, a necessity.  

```Bash
perl /scratch/jandrews/bin/vcfsorter.pl /scratch/jandrews/Ref/hg19.dict merge_vs.vcf > merge_vs.sorted.vcf 2>STDERR
perl /scratch/jandrews/bin/vcfsorter.pl /scratch/jandrews/Ref/hg19.dict merge_samtools.vcf > merge_samtools.sorted.vcf 2>STDERR
```

**iii. Combine the samtools and VarScan VCFs**  
**Bash script (combine_merged_vcfs.sh):**
```Bash
#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N COMBINE_VCFs
#PBS -m e
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=36gb
#PBS -q old

java -Xmx15g -jar /scratch/jandrews/bin/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar \
        -T CombineVariants \
        -R /scratch/jandrews/Ref/hg19.fa \
        --variant:bcftools /scratch/jandrews/Data/Variant_Calling/Coding/FINAL/merge_samtools.sorted.vcf  \
        --variant:varscan /scratch/jandrews/Data/Variant_Calling/Coding/FINAL/merge_vs.sorted.vcf \
        -o /scratch/jandrews/Data/Variant_Calling/Coding/FINAL/Coding_Variants_Combined.vcf \
        -genotypeMergeOptions UNIQUIFY
```


##### 8.) Annotate with VEP. 
This is essentially impossible to get working on the cluster due to how perl is set up on it, so install and run locally. Be sure to use the GrCH37 cache `--port 3337` for hg19, not GrCH38. Motif info is pulled from JASPAR mainly, it seems.  
```Bash
perl ~/bin/ensembl-tools-release-82/scripts/variant_effect_predictor/variant_effect_predictor.pl --everything --vcf --format vcf --fork 2 --symbol --cache --port 3337 -i Coding_Variants_Combined.vcf -o Coding_Variants_Combined.annot.vcf
```

**Fin.** Can do whatever with it at this point. I did an intersect with the GM TF ChIP-Seq data for kicks.

##### 9.) Intersect with GM TF ChIP-Seq data.
bedtools intersect -wa -wb -a /scratch/jandrews/Data/Variant_Calling/Coding/FINAL/Coding_Variants_Combined.annot.vcf -b GM12878_TF151_names_final.bed > Coding_Variants_Combined_Annot_GM_ChIP_TFs_isec.bed

---

### COMPARISON OF SNP ARRAYS, SAMTOOLS, & VARSCAN RESULTS <a name="comparison"></a>
I made the **mistake of assuming** the calls on the arrays (1,2,3,0 aka AA,AB,BB,no call) had allele A as the reference allele and B as the variant. This is apparently not the case. Don't make my mistake, only use the het calls to determine how well your variant calling pipeline is doing vs the arrays.

This picks up after step 1 of the initial variant calling as noted above.

##### 2.) Parse SNP array.
This will create VCFs for each sample column in the summary SNP array file. 
**Bash script (parse_snp_array.sh):**
```Bash
#!/bin/sh
# give the job a name to help keep track of running jobs (optional)
#PBS -N PARSE_SNP_ARRAYS
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=16gb

export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

python3 /scratch/jandrews/bin/parse_snp_array.py /scratch/jandrews/Data/Variant_Calling/SNP_Arrays/GenomeWideSNP_6.na35.annot.csv /scratch/jandrews/Data/Variant_Calling/SNP_Arrays/LYMPHOMASNP77_GTYPE_2014.txt
```

##### 3.) Scrub files.
Remove any potential garbage chromosomes from both array VCFs and those from samtools/varscan.
`(sed '/_g/d' file.vcf | sed '/chrM/d' | sed '/chrY/d') > output.vcf`

or

```Bash
for file in *.vcf; do
	sed -i '/_g/d' "$file";
	sed -i '/chrM/d' "$file";
	sed -i '/chrY/d' "$file";
done
```


##### 4A.) Filter variants.
Figure out which variants on the SNP array contain sufficient read depth (>=5) for each sample to actually be called by samtools or varscan. 
	
First, get read depth at each SNP array variant for each sample.
**Bash script (array_variant_cov.sh):**
```Bash
#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N SNP_ARRAY_COV
#PBS -m e
#PBS -l nodes=1:ppn=8,walltime=8:00:00,vmem=16gb
#PBS -q old

module load bedtools2

for fold in /scratch/jandrews/Data/RNA_Seq/ALIGNED_BAMs/Batch*/; do

	cd "$fold"
	for f in *.bam; do
		echo "$f"
		base=${f%%_*}
		bedtools multicov -bams "$f" -bed /scratch/jandrews/Data/Variant_Calling/Coding/Array_Comparisons/SNP_Array_Ind_Samples/Het_Only/"${base}"_array.vcf > /scratch/jandrews/Data/Variant_Calling/Coding/Array_Comparisons/Cov_At_Each_Array_Variant_In_RNASeq/Het_Only/"${base}"_array_coverage_min5DP.vcf &
	done

done
wait
module remove bedtools2
```

##### 4B.) Actually filter for read-depth. 
Change `$11` as needed (ie. if read count column is column 9, change to `$9`)
```Bash
for file in *.vcf; do 
	base=${file%%_*} ; 
	cat $file | awk '$11 > 4' > "$base"_array_coverage_min5DP.vcf ; 
done
```

##### 5.) Sort each filtered SNP array VCF.
```Bash
for f in *.vcf; do
	base=${f%%.*} ;
	vcf-sort "$f" > ${base}_sorted.vcf ;
done
```

##### 6.) Remove last column (counts) from filtered SNP array VCFs.
```Bash
for f in *.vcf; do
	base=${f%%.*} ;
	cut -f1-10 "$f" > "$base"_cut.vcf ;
done
```

##### 7.) Re-add header.
Have to do this for each filtered SNP array, as bedtools multicov removes it.
```Bash
for f in *.vcf; do
	base=${f%%_*} ;
	sed -i '1s/^/#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	'$base'_ARRAY\n/' "$f" ;
	sed -i '1s/^/##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n/' "$f" ;
	sed -i '1s/^/##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n/' "$f" ;
	sed -i '1s/^/##reference=hg19\n/' "$f" ;
	sed -i '1s/^/##source=parse_snp_array.py\n/' "$f" ;
	sed -i '1s/^/##fileformat=VCFv4.2\n/' "$f" ;
done
```

##### 8.) Zip each VCF.  
```Bash
for f in *.vcf; do
	bgzip -c "$f" > "$f".gz ;
done
```

##### 9. Index each VCF.  
```Bash
for f in *.gz; do
	tabix -p vcf "$f" ;
done
```

##### 10.) Organize.  
Create two directories - one for array/BCFtools intersections, one for array/VarScan intersections. Copy files into each as appropriate.


##### 11.) Intersect files for each sample.
(checks by position only - `-c all` option does this.)
```Bash
module load bcftools-1.2
for f in *_RNAseq*.gz; do
	base=${f%%_*} ;
	bcftools isec -c all -p "$base" "$f" "$base"_array_coverage_min5DP_het.vcf.gz ;
done
module remove bcftools-1.2
```

In this example, 0000.vcf will be records unique to first provided file. 0001.vcf will be records unique to second provided file. 0002.vcf and 0003.vcf are those shared between them.


##### 12.) Count lines for each file, ignoring the headers. Do whatever with that information. Venn diagrams or something.
`grep -v '#' -c file.vcf `


---

### To figure out insert sizes for PINDEL
PINDEL is for figuring out more about indels, inversions, etc, but it's really not meant to run on RNA-Seq data. Regardless, you have to know the insert sizes for your samples to run it, so the below scripts will determine them for you. Our samples have an average/median insert size of 200 bp, though the histograms look more like majority are 150 bp.

**Bash script(get_insert_sizes.sh):**
```Bash
#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N CHECK_INS_SIZES
#PBS -m e
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=32gb

module load java
module load R

for fold in /scratch/jandrews/Data/RNA_Seq/ALIGNED_BAMs/Batch*/; do

		cd "$fold"
		for f in *.bam; do
			echo "$f"
			base=${f%%_*}
			java -jar /export/picard-tools-2.0.1/picard.jar	CollectInsertSizeMetrics INPUT="$f" OUTPUT=/scratch/jandrews/Data/RNA_Seq/Insert_Metrics/"$base"_metrics.txt HISTOGRAM_FILE=/scratch/jandrews/Data/RNA_Seq/Insert_Metrics/"$base"_hist.pdf ;
		wait
		done	
wait
done

module remove R
module remove java
```
