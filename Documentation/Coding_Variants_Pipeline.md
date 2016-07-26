# Coding Variants Pipeline
**Up to date as of 07/26/2016.**  
**Author: jared.andrews07@gmail.com**  

---

This is an imitation of Liv's variant calling pipelines, albeit with some improvements to increase sensitivity and stringency  (hypothetically). Liv tended to treat the variants from RNA-seq and ChIP-seq separately all the way through, and while I do call and filter them slightly differently, I think it's easier to merge them together in the end. This also began as a comparison between the samtools and VarScan variant callers as well, but after the analysis, it seemed the best bet was to simply merge the results from the two callers for the RNA-seq data, as they have fairly high overlap, keeping only the variants called in both.

This was done on the CHPC cluster, so all of the `export`, `source`, and `module load/remove` statements are to load the various software necessary to run the command(s) that follow. If you're running everything locally and the various tools needed are located on your `PATH`, you can ignore these.

> Bash scripts are submitted on the cluster with the `qsub` command. Check the [CHPC wiki](http://mgt2.chpc.wustl.edu/wiki119/index.php/Main_Page) for more info on cluster commands and what software is available. All scripts listed here should be accessible to anyone in the Payton Lab.

All necessary scripts should be here: **N:\Bioinformatics\Jareds_Code**  
They are also in `/scratch/jandrews/bin/` or `/scratch/jandrews/Bash_Scripts/` on the cluster as well as stored on my local PC and external hard drive.  

An _actual_ workflow (Luigi, Snakemake, etc) could easily be made for this with a bit of time, maybe I'll get around to it at some point.

**IMPORTANT**: Be sure to sort and index BAMs before beginning this. There are obviously named bash scripts in the above code folders to do so. All Python scripts here use **Python 3**, so be sure you have the appropriate version installed/set.

---

**Software Requirements:**
- [BEDOPS](http://bedops.readthedocs.org/en/latest/index.html)
- [Samtools, BCFtools, and HTSlib](http://www.htslib.org/)  
  - These should be available on the CHPC cluster.
- [Python3](https://www.python.org/downloads/)
  - Use an [anaconda environment](http://mgt2.chpc.wustl.edu/wiki119/index.php/Python#Anaconda_Python) if on the CHPC cluster. (Also useful for running various versions of python locally).  
- [bedtools](http://bedtools.readthedocs.org/en/latest/)
  - Also available on the CHPC cluster.
- [VarScan2](http://dkoboldt.github.io/varscan/)
- [GATK](https://www.broadinstitute.org/gatk/download/)
- [Picard Tools](http://broadinstitute.github.io/picard/)
- [Perl](https://www.perl.org/)
- [Variant Effect Predictor](http://useast.ensembl.org/info/docs/tools/vep/index.html)

---

**Sections:**  
I recommend going through the first three of these in order, though the RNA-seq & ChIP-seq data can be processed in parallel up to a certain point.
- [Variant Calling from RNA-seq Data](#variant-calling-from-rna-seq-data)
- [Variant Calling from ChIP-seq Data](#variant-calling-from-chip-seq-data)
- [Merging Variants from ChIP-seq/RNA-seq Data](#merging-all-variants)
- [Creating Mutational Signatures](#create-mutational-signatures)
- [Intersections with Genomic Features](#intersections-with-features-of-interest)

---

## Variant Calling From RNA-seq
This section will call variants from RNA-seq data using both VarScan and bcftools (samtools). 

#### 0.) Mark duplicates if you haven't already.
PCR duplicates and such can affect variant calling in a big way, apparently. We mark them here so that `samtools mpileup` will ignore them (which it does by default once they're marked).

**Bash script (mark_duplicates.sh):**  
```Bash
#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N MARK_DUPS
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=32gb

module load java
module load R

for fold in /scratch/jandrews/Data/RNA_Seq/ALIGNED_BAMs/Batch*/; do

	cd "$fold"
	for f in *.sorted.bam; do
		echo "$f"
		base=${f%%.*}
		java -jar /export/picard-tools-2.0.1/picard.jar	MarkDuplicates INPUT="$f" OUTPUT="$base".dups_rmvd.bam METRICS_FILE="$base".dup_metrics.txt ;
	wait
	done	
wait
done

module remove R
module remove java
```

#### 1A.) Call variants with VarScan.
samtools: mpileup piped to varscan to call variants with filters for read depth (10) and quality (15). Get to submit this guy for every file, fun fun. Could be done in a better way, but I'm lazy and don't want to wait.

**Bash script (var_call_varscan.sh):**
```Bash
#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N Var_Call_VarScan
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=48gb

module load samtools

for file in /scratch/jandrews/Data/RNA_Seq/ALIGNED_BAMs/Batch28/*dups_marked.bam; do

	samtools mpileup -t DP -f /scratch/jandrews/Ref/hg19.fa $file | java -Xmx15g -jar /scratch/jandrews/bin/VarScan.v2.3.9.jar mpileup2cns --min-coverage 10 --min-avg-qual 15 --variants 1 --output-vcf 1 >"$file".vs.vcf  &
	
done
wait
module remove samtools
```

#### 1B i.) Call variants with bcftools (samtools).
**Bash script (var_call_bcf.sh):**
```Bash
#!/bin/sh
# give the job a name to help keep track of running jobs (optional)
#PBS -N Var_Call_RNA
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=8:00:00,vmem=16gb

module load samtools
module load bcftools

for file in /scratch/jandrews/Data/RNA_Seq/ALIGNED_BAMs/Batch9/*.bam; do

	samtools mpileup -u -t DP -f /scratch/jandrews/Ref/hg19.fa $file | bcftools call -cv -O v - | /scratch/jandrews/bin/vcfutils.pl varFilter -D100 > $file.vcf &
	
done
wait
module remove samtools
module remove bcftools
```


#### 1B ii.) 
Filter the BCFtools VCFs for the same quality and read-depth that VarScan used.

```Bash
module load bcftools
for f in *.vcf; do
    base=${f##*/}
    bcftools filter -i 'DP>=10 & QUAL>=15' "$f" > ${base%.*}.filtered.vcf
done
module remove bcftools
```

At this point, I also did a comparison between the SNP arrays with these callers, which is explained [below](#comparison) after the rest of the pipeline. 


#### 2.) Clean VCFs.
```Bash
for file in *.vcf; do
	base=${file%%_*} ;
	(sed '/_g\|chrM\|chrY\|chrX\|chr23/d' "$file") > "$base"_RNAseq_VS_clean.vcf ;
done
```

#### 3.) Sort all VCFs.
```Bash
for f in *.vcf; do
	base=${f%%_*} ;
	vcf-sort "$f" > ${base}_RNAseq_VS_sorted.vcf ;
done
```

#### 4.) Zip and index each VCF.
```Bash
for f in *.vcf; do
	bgzip -c "$f" > "$f".gz;
done

for f in *.gz; do
	tabix -p vcf "$f" ;
done
```



#### 5.) Merge VCFs.
First, those from VarScan with each other. Then those from BCFTools with each other. `-m none` means multiallelic records will be split to separate lines. Doing so is rather important, as if two samples have different variant alleles at the same position, only one is reported as having the variant if multiallelic records are allowed. Alternatively, setting `-m both` should create a multiallelic record, which may be wanted at times. No idea what the default is, BCFtools docs don't mention.

**For VarScan:**
```Bash
module load bcftools
bcftools merge -O v -m none --force-samples -i ADP:sum *.gz > merge.vcf
module remove bcftools
```

**For BCFTools:**
```Bash
module load bcftools
bcftools merge -O v -m none -i DP:sum *.gz > merge.vcf
module remove bcftools
```

#### 6.) Fix headers.
Text editor style because I was too lazy to write something.

Now set these files aside and let's work on the ChIP-seq data.

---

## Variant Calling from ChIP-Seq Data
Here's how I identify variants from the ChIP-seq data. This is fairly stringent due to the low coverage, but hopefully it reduces false positives and ensures we don't waste time trying to validate mutations that aren't around.

**First**, call and scrub variant files the same as we did for the RNA-seq data above, but do **not** filter for read-depth or quality yet.

#### 1.) Organize.
Copy all VCF.gz and VCF.gz.tbi files for a given sample into an empty directory.  

#### 2.) Merge the VCFs for each sample.
This *sums* read depth (DP) and *averages* quality values (QV) for each file (histone mark) used for the sample.  
**Bash script (merge_samp_vcfs.sh)**:
```Bash
#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N MERGE_VCFs
#PBS -m e
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=32gb
#PBS -q old

module load bcftools-1.2
for fold in /scratch/jandrews/Data/Variant_Calling/Non_Coding/NEW_VCFS_wQVs/IND_SAMP_MERGED/*/; do
	cd "$fold"
	base=${PWD##*/}
	mkdir "$base"_IND_VCFS
	mv *vcf* "$base"_IND_VCFS
	cd "$base"_IND_VCFS/
	bcftools merge -O v -m none -i DP:sum,QV:avg,DP4:sum *.gz > ../"$base"_variants.vcf
done

module remove bcftools-1.2
```

#### 3.) Fix header of the resulting merged VCF for each sample.
Shortens the column headers to (sample_mark). This also fixes other header inconcistencies like case, extra underscores, etc.  
**Python script (shorten_vcf_header.py)**:
```Bash
export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda
for fold in /scratch/jandrews/Data/Variant_Calling/Non_Coding/VCFs_With_Quals/Ind_Samp_Merged/*/; do
	cd "$fold"
	python3 /scratch/jandrews/bin/shorten_vcf_header.py ${PWD##*/}_variants.vcf ${PWD##*/}_variants_fixed.vcf
done
```

#### 4.) Zip and index the VCF for each sample, then filter for >=10 summed DP.
Those that don't meet the read depth cutoff are removed.  

**Bash script (process_samp_vcfs.sh):**  
```Bash
#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N PROCESS_VCFs
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=4gb

module load bcftools-1.2

for fold in /scratch/jandrews/Data/Variant_Calling/Non_Coding/VCFs_With_Quals/Ind_Samp_Merged/*/; do
	cd "$fold"
	bgzip -c ${PWD##*/}_variants.vcf > ${PWD##*/}_variants.vcf.gz
	tabix -p vcf ${PWD##*/}_variants.vcf.gz
    bcftools filter -i 'DP>=10' ${PWD##*/}_variants.vcf.gz > ${PWD##*/}_variants_filtered.vcf
done

module remove bcftools-1.2
```

#### 5.) Remove single-mark variants from the VCF for each sample. 
We're only interested in those that occur in multiple data types (histone marks) for a given sample so that we can be sure they aren't artifacts of a specific histone mark, etc. As such, this removes those that only occur in a single histone mark (or data type).  

**NOTE:** _Format of samples in the header must be **'sample_mark'**._

**Python script (filter_nc_vcf.py):**
```Bash
export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda
for fold in /scratch/jandrews/Data/Variant_Calling/Non_Coding/NEW_VCFs_wQVs/IND_SAMP_MERGED/*/; do
	cd "$fold"
	python3 /scratch/jandrews/bin/filter_nc_vcf.py ${PWD##*/}_variants_filtered.vcf ${PWD##*/}_variants_filtered_multimark.vcf
done
```

---

## Merging All Variants  
Now we can put the variants from the RNA-seq and ChIP-seq data together. Though Liv previously treated these sets separately, it becomes a *real* hassle to maintain several files for each sample. I have instead merged these sets for each sample, yielding a final, single VCF for each. 

#### 1.) Pool files.
Stick the VS, BCF, and filtered, multimark variant files from ChIP-seq data into the same folder.

#### 2.) Fix headers.
Be sure the file names are `sample.restofname.vcf` for each sample. 

#### 3.) Combine the BCFTools, VarScan, and noncoding files.
This step was an **enormous** hassle to figure out. This will yield a single file for each sample with **all** variants in the sample.

**i. Create a sequence dict for reference genome**  
This has to be done for GATK/Picard tools to work properly.   
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
#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N SORT_VCF_FOR_GATK
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=16gb
#PBS -q old

module load java

for f in /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/MERGED/*.vcf; do
	java -Xmx8g -jar /scratch/jandrews/bin/picard-tools-2.2.1/picard.jar SortVcf \
	I="$f" \
	O="$f".sorted \
	SEQUENCE_DICTIONARY=/scratch/jandrews/Ref/hg19.dict 
	rename .vcf.sorted .sorted.vcf "$f".sorted
done

module remove java
```

**iii. Combine the samtools, VarScan, and noncoding VCFs for each sample.**  
You have to change the `samp` line below to match the sample for each file. Be sure to check all samples were combined properly. If a sample doesn't have all three sets, adjust as necessary to combine.

**Bash script (combine_samp_vcfs.sh):**
```Bash
#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N COMBINE_VCFs
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=16gb

samp=FL120

module load java
java -Xmx8g -jar /scratch/jandrews/bin/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar \
	-T CombineVariants \
	-R /scratch/jandrews/Ref/hg19.fa \
	--variant:bcftools /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/scratch/"$samp".RNAseq_BCF.sorted.vcf  \
	--variant:varscan /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/scratch/"$samp".RNAseq_VS.sorted.vcf \
	--variant:chip_seq /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/scratch/"$samp".variants_filtered_multimark.sorted.vcf \
	-o /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/scratch/"$samp".Combined.vcf \
	-genotypeMergeOptions UNIQUIFY 

module remove java
```

#### 4.) Filter those only found in one data set.  
This will remove SNPs found only by one of the callers in the RNA-seq data and not corroborated by the ChIP-seq data. Those found in only the ChIP-seq set are retained, as they had to be called in two ChIP histone markers to even be included in that set.

**Bash script (filter_singleset_vcfs.sh):**  
```Bash
#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N FILTER_SINGLESET_VCFs
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=16gb

# This script filters variants found by only one of the two callers used. A variant has to be found in two datasets (already done for chip_seq),
# so the set would have to be "bcftools-varscan", "chip_seq", "varscan-chip_seq", or "bcftools-chip_seq" to be included in the output.

module load bcftools

for f in /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/*Combined.vcf; do
    bcftools filter -e 'INFO/set="bcftools" || INFO/set="varscan"' "$f" > "$f".filtered ;
    rename .vcf.filtered .multiset.vcf "$f".filtered
done

module remove bcftools
```

#### 5.) Sort the files for GATK.
I usually move these files into a new directory just so things are a bit easier to handle, as we probably have a ton of files at this point. Just need the plain `.VCF` files, not compressed or anything. GATK is dumb in that it requires the VCFs sorted similarly to the reference dict used for it (which we created previously). So we have to sort the VCFs again to ensure they are in the correct order.

> This can be annoying as hell for the ChIP_seq only files. For whatever reason, the order of the contigs in the header are sometimes not changed, so you may have to reorder them manually for the first file or two, then trying to sort again seems to work.

```Bash
module load java
for f *.vcf; do
	perl /scratch/jandrews/bin/vcfsorter.pl /scratch/jandrews/Ref/hg19.dict "$f">"$f".sorted 2>STDERR
done
```

**Note:** If you have trouble with this step, you can try using Picard SortVcf, which can sort to a reference dict, though I had trouble getting it to play nice.

#### 6.) Merge the combined files.  
Now we can merge the VCFs for each sample into a single big file. I usually move these files into a new directory just so things are a bit easier to handle, as we probably have a ton of files at this point. We can use the samples that only used ChIP-seq data as well, though we may have to edit the header for consistency after merging. Have to stick full paths to each `VCF` here, which is annoying, but you can use a loop and `echo` along with copy and pasting to make it painless.

*To get arguments for `CombineVariants` GATK tool:*

```Bash
for f in *.vcf; 
	do base=${f%%.*}; 
	echo --variant:"$base" /scratch/jandrews/Data/Variant_Calling/Coding_Noncoding_Merged/COMBINED_VARS_FINAL/MERGED/"$f" ; 
done
```

**Bash script (combine_all_vcfs.sh)** 
Not even bothering to paste this one, since it's just a million file paths and the CombineVariants GATK tool again. You'll probably want ot manually clean up the header after this - it's likely to be a bit messy.


#### 7.) Filter out common variants.
This will remove the common variants (those with a MAF >0.01 in dbSNP build 146). I originally used the internal VEP filters that use the 1000 genomes project allele frequencies, but found that they don't have all the ones that dbSNP does. Since I use hg19 from UCSC as my reference genome, I had to download the common snps (146) track from UCSC as a `bed` file through the table browser to ensure correct positions. Bedtools is *very* memory in-efficient when using a large file for `-b` as we are here, hence why I go ahead and use an interactive session with a ton of memory.

```Bash
qsub -I -l nodes=1:ppn=1,walltime=4:00:00,vmem=128gb 

module load bedtools2

for f in *.vcf.gz; do 
	bedtools intersect -v -header -a "$f" -b /scratch/jandrews/Ref/dbSNP146_common_variants.bed.gz > "${f%%.*}".dbSNP146_common_rmvd.vcf; 
done
```

#### 8.) Annotate with VEP. 
This is essentially impossible to get working on the cluster due to how perl is set up on it, so install and run locally. Be sure to use the GrCH37 cache `--port 3337` for hg19, not GrCH38. Motif info is pulled from JASPAR mainly, it seems.  

```Bash
for f in *.gz; 
	do perl ~/bin/ensembl-tools-release-82/scripts/variant_effect_predictor/variant_effect_predictor.pl --fork 2 --check_existing --biotype --gencode_basic --hgvs --canonical --uniprot --variant_class --gmaf --maf_1kg --maf_esp --polyphen b --regulatory --sift b --species homo_sapiens --symbol --cache --port 3337 --vcf --stats_file "$f".stats.html --input_file "$f" -o "$f".VEP_Anno; 
	rename .vcf.gz.VEP_Anno .VEP_Anno.vcf "$f".VEP_Anno; 
done
```

After this point, you're on your own. You have individual sample file and one enormous master table.

---

## Create mutational signatures  
This section tries to identify mutational signatures from our variants, allowing us to make some inferences as to their generating processes. See the [Alexandrov paper](http://www.nature.com/nature/journal/v500/n7463/full/nature12477.html) for some more details about this idea.

#### 1.) Create frequency matrix for SNVs.
We'll use this matrix to generate the mutational signatures for our samples. Throw all the `<sample>_Combined.vcf` files you want to include into a directory by themselves. 

**Python script (make_trinuc_matrix.py):**  
```Bash

python /scratch/jandrews/bin/make_trinuc_matrices.py -i /vcf_directory -o FLDLCLL_CCCB.txt -r /scratch/jandrews/Ref/hg19.fa
```

#### 2.) Create mutational signatures.
I use R studio and the [SomaticSignatures](http://www.bioconductor.org/packages/devel/bioc/vignettes/SomaticSignatures/inst/doc/SomaticSignatures-vignette.html) package for this. It is relatively straight-forward, so just read the link and you'll be able to figure it out. Just import the matrix file we created in the last step, convert it to a matrix, and plug it into the commands above. A script could easily be written to do this if you were so inclined.

Compare the output figures to the [COSMIC mutational signatures](http://cancer.sanger.ac.uk/cosmic/signatures) or do whatever you want with them. 

#### 3.) (Optional) Run it again for the samples grouped by cell type.
If you want to try to show clear differences between the cell types, you can merge the VCFs for each cell type using `bcftools` and just run it on those.


---

## Variant Caller Comparisons
I made the **mistake of assuming** the calls on the arrays (1,2,3,0 aka AA,AB,BB,no call) had allele A as the reference allele and B as the variant. This is apparently not the case. Don't make my mistake, only use the het calls to determine how well your variant calling pipeline is doing vs the arrays.

This picks up after step 1 of the initial variant calling as noted above.

#### 1.) Parse SNP array.
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

python /scratch/jandrews/bin/parse_snp_array.py /scratch/jandrews/Data/Variant_Calling/SNP_Arrays/GenomeWideSNP_6.na35.annot.csv /scratch/jandrews/Data/Variant_Calling/SNP_Arrays/LYMPHOMASNP77_GTYPE_2014.txt
```

#### 2.) Scrub files.
Remove any chromosomes we don't want/care about from both array VCFs and those from samtools/varscan.
```Bash
sed -i '/_g\|chrM\|chrY/d' file.vcf > output.vcf
```


#### 3A.) Filter variants.
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

#### 3B.) Actually filter for read-depth. 
Change `$11` as needed (ie. if read count column is column 9, change to `$9`)
```Bash
for file in *.vcf; do 
	base=${file%%_*} ; 
	cat $file | awk '$11 > 4' > "$base"_array_coverage_min5DP.vcf ; 
done
```

#### 4.) Sort each filtered SNP array VCF.
```Bash
for f in *.vcf; do
	base=${f%%.*} ;
	vcf-sort "$f" > ${base}_sorted.vcf ;
done
```

#### 5.) Remove last column (counts) from filtered SNP array VCFs.
```Bash
for f in *.vcf; do
	base=${f%%.*} ;
	cut -f1-10 "$f" > "$base"_cut.vcf ;
done
```

#### 6.) Re-add header.
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

#### 7.) Zip & index each VCF.  
```Bash
for f in *.vcf; do
	bgzip -c "$f" > "$f".gz ;
	tabix -p vcf "$f".gz ;
done
```

#### 8.) Organize.  
Create two directories - one for array/BCFtools intersections, one for array/VarScan intersections. Copy files into each as appropriate.


#### 9.) Intersect files for each sample.
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


##### 10.) Count lines for each file, ignoring the headers. Do whatever with that information. Venn diagrams or something.
`grep -v '#' -c file.vcf `


---

## To figure out insert sizes for PINDEL
PINDEL is for figuring out more about indels, inversions, etc, but it's really not meant to run on RNA-Seq data. Regardless, you have to know the insert sizes for your samples to run it, so the below script will determine them for you. Our samples have an average/median insert size of 200 bp, though the histograms look more like majority are 150 bp.

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

## Intersections with features of interest

Intersect the resulting multimark, filtered variants file for each sample with SEs, regular enhancers, and TSSs to determine the percentage of SNVs for file found in each. Make sure to always sort and remove dups first, as it'll skew your results otherwise. Some of the MMPIDs and such may have duplicates due to the way they were acquired, hence why this is necessary.

This same principle can be applied to intersections with whatever you want (TF motifs, regulatory elements, etc). Here I show intersections with **super enhancers, regular enhancers (FAIRE-positive, non-TSS MMPIDs that aren't located within SEs), TSSs, and TSSs located within SEs.**

#### 1.) Intersect TSSs with SEs to get those that lie within/outside the SEs. 
The TSS file here can be found on the Payton Shared Drive in the master files folder (and probably about ten other places as well).

```Bash
sort -k1,1 -k2,2n 2kbTSStranscripts_gencode19_protein_coding_positions.bed \
| uniq - > 2kbTSStranscripts_gencode19_protein_coding_positions_uniq.bed

module load bedtools2

bedtools intersect -wa -a 2kbTSStranscripts_gencode19_protein_coding_positions.bed -b ../SES/ALL_SES_POSITIONS_SORTED.bed > 2kbTSStranscripts_gencode19_protein_coding_inSEs.bed
bedtools intersect -v -a 2kbTSStranscripts_gencode19_protein_coding_positions.bed -b ../SES/ALL_SES_POSITIONS_SORTED.bed > 2kbTSStranscripts_gencode19_protein_coding_outsideSEs.bed
```  

  
#### 2.) Intersect FAIRE-positive MMPID positions with SEs to remove those that overlap  
Sort and remove dups (shouldn't be any in MMPID positions, though intersects may result in some - an MMPID in two SEs).
```
sort -k1,1 -k2,2n MMPID_NonTSS_FAIRE_POSITIVE_POSITIONS.bed | uniq - > MMPID_NonTSS_FAIRE_POSITIVE_POSITIONS_UNIQ.bed

bedtools intersect -wa -a MMPID_NonTSS_FAIRE_POSITIVE_POSITIONS_UNIQ.bed -b ../SES/ALL_SES_POSITIONS_SORTED.bed \
| uniq - > MMPID_NonTSS_FAIRE_POSITIVE_inSEs.bed

bedtools intersect -v -a MMPID_NonTSS_FAIRE_POSITIVE_POSITIONS_UNIQ.bed -b ../SES/ALL_SES_POSITIONS_SORTED.bed \
| uniq - > MMPID_NonTSS_FAIRE_POSITIVE_outsideSEs.bed
```

#### 3.) Intersect the variants with features of interest
This intersects the variants with the SEs, MMPIDs, and TSS positions created in previous steps. This will create a summary file for each sample with counts for each intersection. This filters out variants found near the Ig loci as well, though it does so in a lazy way and just removes all the variants in the region rather than only those that overlap with the Ig genes. Can edit  `/scrach/jandrews/Ref/Ig_Loci.bed` to change these positions. 

The script below is set up to intersect only the ChIP-seq variants, but you can edit it to use your combined variants file for each sample instead. This is also using bed files, but we can use `VCFs` since `VEP` doesn't break the format like `FunSeq` does, so `bedtools` should be able to handle them fine.

**Bash script (isec_ind_samp_variants_wMMPIDs_SEs_TSS.sh):**
```Bash
#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N ISEC_EVERYTHING
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=8gb

module load bedtools2

rm /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/all_samples_isecs_summary.txt
touch /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/all_samples_isecs_summary.txt

for f in /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/*positions.bed; do

	base=${f##*/}

	sort -k1,1 -k2,2n "$f" \
	| uniq - > /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_sorted.bed
	echo -en ${base%.*}' funseq annotated variants: ' >> /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/all_samples_isecs_summary.txt
	var_count=$(wc -l < $f)
	echo -e $var_count >> /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/all_samples_isecs_summary.txt

	bedtools intersect -v -a /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_sorted.bed -b /scratch/jandrews/Ref/Ig_Loci.bed \
	| uniq - > /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_sorted_funseq_noIgLoci.bed

	echo -en ${base%.*}' funseq annotated variants (Ig filtered): ' >> /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/all_samples_isecs_summary.txt
	var_count=$(wc -l < /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_sorted_funseq_noIgLoci.bed)
	echo -e $var_count >> /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/all_samples_isecs_summary.txt

	bedtools intersect -wa -a /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_sorted_funseq_noIgLoci.bed -b /scratch/jandrews/Data/Variant_Calling/Non_Coding/SE_TSS_MMPID_INTERSECTS/MMPIDS/MMPID_NonTSS_FAIRE_POSITIVE_outsideSEs.bed \
	| uniq - > /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_funseq_inMMPIDS_outsideSEs_noIgLoci.bed

	echo -en ${base%.*}' funseq annotated variants in regular enhancers (Ig filtered): ' >> /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/all_samples_isecs_summary.txt
	var_count=$(wc -l < /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_funseq_inMMPIDS_outsideSEs_noIgLoci.bed)
	echo -e $var_count >> /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/all_samples_isecs_summary.txt

	bedtools intersect -wa -a /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_sorted_funseq_noIgLoci.bed -b /scratch/jandrews/Data/Variant_Calling/Non_Coding/SE_TSS_MMPID_INTERSECTS/SES/ALL_SES_POSITIONS_SORTED.bed \
	| uniq - > /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_funseq_inSEs_noIgLoci.bed

	echo -en ${base%.*}' funseq annotated variants in super enhancers (Ig filtered): ' >> /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/all_samples_isecs_summary.txt
	var_count=$(wc -l < /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_funseq_inSEs_noIgLoci.bed)
	echo -e $var_count >> /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/all_samples_isecs_summary.txt

	bedtools intersect -wa -a /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_sorted_funseq_noIgLoci.bed -b /scratch/jandrews/Data/Variant_Calling/Non_Coding/SE_TSS_MMPID_INTERSECTS/TSS/2kbTSStranscripts_gencode19_protein_coding_positions_uniq.bed \
	| uniq - > /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_funseq_in2kbTSS_noIgLoci.bed

	echo -en ${base%.*}' funseq annotated variants in 2kb TSSs (Ig filtered): ' >> /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/all_samples_isecs_summary.txt
	var_count=$(wc -l < /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_funseq_in2kbTSS_noIgLoci.bed)
	echo -e $var_count >> /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/all_samples_isecs_summary.txt

	bedtools intersect -wa -a /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_sorted_funseq_noIgLoci.bed -b /scratch/jandrews/Data/Variant_Calling/Non_Coding/SE_TSS_MMPID_INTERSECTS/TSS/2kbTSStranscripts_gencode19_protein_coding_inSEs.bed \
	| uniq - > /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_funseq_in2kbTSS_inSEs_noIgLoci.bed

	echo -en ${base%.*}' funseq annotated variants in 2kb TSSs in SEs (Ig filtered): ' >> /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/all_samples_isecs_summary.txt
	var_count=$(wc -l < /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/${base%.*}_funseq_in2kbTSS_inSEs_noIgLoci.bed)
	echo -e $var_count >> /scratch/jandrews/Data/Variant_Calling/Non_Coding/Annotated_Results/FLDL_Only_Ind_Filtering/all_samples_isecs_summary.txt

done
module remove bedtools2
```
