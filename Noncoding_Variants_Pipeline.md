# Noncoding Variant Pipeline
##### Up to date as of 02/24/2015.
This is an imitation of Liv's Noncoding Variant Calling Pipeline. Newer versions of the software, but ideally similar results. Much of this could likely be done more elegantly/efficiently, but I'm no wizard.   

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
  - Use an [anaconda environment](http://mgt2.chpc.wustl.edu/wiki119/index.php/Python#Anaconda_Python) if on the CHPC cluster (also useful for running various versions of python locally).  
- [bedtools](http://bedtools.readthedocs.org/en/latest/)
  - Also available on the CHPC cluster.


---
##### 1.) Call variants and generate VCFs for all BAM files for all samples.  
The BAMs should be broken up into batches as shown here if you want them to run **asynchronously**. If time isn't important to you, just delete the `&` and `wait` from the script below and throw all the BAMs in the same folder. Just be sure to allocate enough time to the script, as it'll take a while (days) when done this way.
>The `&` at the end of a command in a bash script means it **will not wait for the command to complete** before moving to the next.   

This is great when you can allocate enough processors and RAM to handle all the commands running at once, but samtools/BCFtools eat a fair amount of RAM. Rather than sit in the queue for hours waiting for the 48 processors and 128 GB of RAM you requested to handle all your files, you're usually better off breaking them into batches and submitting multiple times as shown here.

**Bash script (var_call_bcf.sh)**:  
```Bash
#!/bin/sh
# give the job a name to help keep track of running jobs (optional)
#PBS -N Var_Call
#PBS -m e
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=12gb

module load samtools-1.2
module load bcftools-1.2

for file in /scratch/jandrews/Data/ChIP_Seq/BAMs/K27AC/Batch1/*.bam; do

	samtools mpileup -u -t DP -f /scratch/jandrews/Ref/hg19.fa $file \
	| bcftools call -cv -O v \
	| /scratch/jandrews/bin/vcfutils.pl varFilter -D100 > $file.vcf &
	
done
wait
module remove samtools-1.2
module remove bcftools-1.2
```  

##### 2.) Add qual value for each sample to the INFO field of its VCF.
INFO fields can have operations done on them during the merge, whereas for the QUAL score, the max value across samples are taken (while we want the average). So we add a custom info field.

```Bash
export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

for file in *.vcf; do 
	python3 /scratch/jandrews/bin/add_qual_to_info_vcf.py "$file" ; 
done
```

##### 3.) Remove garbage chromosomes from VCF.
```Bash
for file in *.vcf; do
	sed -i '/_g\|chrM\|chrY/d' "$file";
done
```

##### 4.) Sort VCFs.
```Bash
for f in *.vcf; do
	vcf-sort "$f" > "$f".sorted ;
done
```

##### 5.) Zip and index each VCF.
```Bash
for f in *.vcf; do
	bgzip -c "$f" > "$f".gz;
done

for f in *.gz; do
	tabix -p vcf "$f" ;
done
```

##### 6.) Copy all VCF.gz and VCF.gz.tbi files for a given sample into an empty directory.  

##### 7.) Merge the VCFs for each sample.
This sums read depth (DP) and averages quality values (QV) for each file (histone mark) used for the sample.  
**Bash script (merge_samp_vcfs.sh)**:
```Bash
#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N MERGE_VCFs
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=8gb

module load bcftools-1.2
for fold in /scratch/jandrews/Data/Variant_Calling/Non_Coding/VCFs_With_Quals/Ind_Samp_Merged/*/; do
	cd "$fold"
	base=${PWD##*/}
	cd ${PWD##*/}_Ind_VCFs/
	bcftools merge -O v -m none -i DP:sum,QV:avg,DP4:sum *.gz > ../"$base"_variants.vcf
done

module remove bcftools-1.2
```

##### 8.) Fix header of the resulting merged VCF for each sample.
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

##### 9.) Zip and index the VCF for each sample, then filter for >=10 summed DP.
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

##### 10.) Remove single-mark variants from the VCF for each sample. 
We're only interested in those that occur in multiple data types (histone marks) for a given sample so that we can be sure they aren't artifacts of a specific histone mark, etc. As such, this removes those that only occur in a single histone mark (or data type).  

**NOTE:** _Format of samples in the header must be **'sample_mark'**._

**Python script (filter_nc_vcf.py):**
```Bash
export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda
for fold in /scratch/jandrews/Data/Variant_Calling/Non_Coding/VCFs_With_Quals/Ind_Samp_Merged/*/; do
	cd "$fold"
	python3 /scratch/jandrews/bin/filter_nc_vcf.py ${PWD##*/}_variants_filtered.vcf ${PWD##*/}_variants_filtered_multimark.vcf
done
```

##### 11.) Convert VCFs to BED files for easy intersecting later. 
Will use these to determine those variants removed by FunSeq after annotation and generate a list of all variants, annotated or not.

```Bash
export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda
for fold in /scratch/jandrews/Data/Variant_Calling/Non_Coding/VCFs_With_Quals/Ind_Samp_Merged/*/; do
	cd "$fold"
	bedops vcf2bed --do-not-split-alt-alleles < ${PWD##*/}_variants_filtered_multimark.vcf > ${PWD##*/}_variants_filtered_multimark.bed
done
```

##### 12.) Annotate with [FunSeq2](http://funseq2.gersteinlab.org/analysis) using a MAF of 0.01 to remove known SNPs found in the 1000 genomes project. 
Use BED output, as the VCF output isn't formatted properly. Upload all VCFs for the samples you want to compare at once. Click the green button to add additional files. As far as I know, there is no limit. Download the output files.


##### 13.) Parse the FunSeq output (Output.BED). 
This will remove any variants that are said to affect coding sequences by FunSeq. It prints two files for each sample - one with the positions of the variants, the other with the full lines for the sample variants that meet the above requirements. To exclude those found in any of the other samples provided to FunSeq as well, use the `-uniq` option.

**Python script (parse_funseq.py):**  
`python3 /scratch/jandrews/bin/parse_funseq.py FLDL_filtered_multitype_funseq_variants.bed`

  
  
# Intersections with features of interest

Intersect the resulting multimark, filtered variants for each sample with SEs, regular enhancers, and TSSs to determine percentage of SNVs for file found in each. Make sure to always sort and remove dups first, as it'll skew your results otherwise. Some of the MMPIDs and such may have duplicates due to the way they were acquired, hence why this is necessary.

This same principle can be applied to intersections with whatever you want (TF motifs, regulatory elements, etc). Here I show intersections with **super enhancers, regular enhancers (FAIRE-positive MMPIDs that aren't located within SEs), TSSs, and TSSs located within SEs.**

##### 1.) Intersect TSSs with SEs to get those that lie within/outside the SEs. 
The TSS file here can be found on the Payton Shared Drive in the master files folder (and probably about ten other places as well).

```Bash
sort -k1,1 -k2,2n 2kbTSStranscripts_gencode19_protein_coding_positions.bed \
| uniq - > 2kbTSStranscripts_gencode19_protein_coding_positions_uniq.bed

module load bedtools2

bedtools intersect -wa -a 2kbTSStranscripts_gencode19_protein_coding_positions.bed -b ../SES/ALL_SES_POSITIONS_SORTED.bed > 2kbTSStranscripts_gencode19_protein_coding_inSEs.bed
bedtools intersect -v -a 2kbTSStranscripts_gencode19_protein_coding_positions.bed -b ../SES/ALL_SES_POSITIONS_SORTED.bed > 2kbTSStranscripts_gencode19_protein_coding_outsideSEs.bed
```  

  
##### 2.) Intersect FAIRE-positive MMPID positions with SEs to remove those that overlap  
Sort and remove dups (shouldn't be any in MMPID positions, though intersects may result in some - an MMPID in two SEs).
```
sort -k1,1 -k2,2n MMPID_NonTSS_FAIRE_POSITIVE_POSITIONS.bed | uniq - > MMPID_NonTSS_FAIRE_POSITIVE_POSITIONS_UNIQ.bed

bedtools intersect -wa -a MMPID_NonTSS_FAIRE_POSITIVE_POSITIONS_UNIQ.bed -b ../SES/ALL_SES_POSITIONS_SORTED.bed \
| uniq - > MMPID_NonTSS_FAIRE_POSITIVE_inSEs.bed

bedtools intersect -v -a MMPID_NonTSS_FAIRE_POSITIVE_POSITIONS_UNIQ.bed -b ../SES/ALL_SES_POSITIONS_SORTED.bed \
| uniq - > MMPID_NonTSS_FAIRE_POSITIVE_outsideSEs.bed
```

##### 3.) Intersect the variants with features of interest
This intersects the funseq annotated variants with the SEs, MMPIDs, and TSS positions created in previous steps. This will create a summary file for each sample with counts for each intersection. This filters out variants found near the Ig loci as well, though it does so in a lazy way and just removes all the variants in the region rather than only those that overlap with the Ig genes. Can edit  `/scrach/jandrews/Ref/Ig_Loci.bed` to change these positions.  

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

#### Done.
Check the summary file for a sample by sample breakdown of the intersections.  
>"The only difference between science and screwing around is writing things down."
