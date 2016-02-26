Up to date as of 01/15/2015.
Jared's imitation of Liv's Coding variant calling pipeline. This also began as a comparison between the samtools and varscan variant callers as well,
but after the analysis, it seemed the best bet was to simply merge the results from the two callers, as they have fairly high overlap.

All necessary scripts should be in Jared's code folder: N:\Bioinformatics\Jareds_Code

IMPORTANT: Be sure to sort and index BAMs before beginning this. Bash scripts in code folder to do so.

1A.) samtools: mpileup piped to varscan to call variants with filter for read depth (5) and quality (15).
Bash script (var_call_varscan.sh):
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


1B i.) Samtools piped into BCFtools (old calling method that Liv used - '-c' option).
Bash script (var_call_bcf_orig.sh):
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


or with new method (supposedly better for multiallelic/rare variant finding).


Bash script (var_call_bcf_new.sh):
#!/bin/sh

# give the job a name to help keep track of running jobs (optional)
#PBS -N Var_Call_BCF_NEW
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=8:00:00,vmem=16gb

module load samtools-1.2
module load bcftools-1.2

for file in /scratch/jandrews/Data/RNA_Seq/ALIGNED_BAMs/Batch9/*.bam; do

	base=${f##*/}
	samtools mpileup -u -t DP -f /scratch/jandrews/Ref/hg19.fa $file | bcftools call -mv -O v - | /scratch/jandrews/bin/vcfutils.pl varFilter -D100 > /scratch/jandrews/Data/Variant_Calling/Coding/BCFTools/New/VCFs/${base%.*}.vcf &
	
done
wait
module remove samtools-1.2
module remove bcftools-1.2


1B ii.) Filter the BCFtools VCFs for the same quality and read-depth that VarScan used.

module load bcftools-1.2
for f in *.vcf; do
    base=${f##*/}
    bcftools filter -i 'DP>=5 & QUAL>=15' "$f" > ${base%.*}.filtered.vcf
done
module remove bcftools-1.2


###-COMPARISON OF SNP ARRAYS, SAMTOOLS, & VARSCAN RESULTS-###
I made the mistake of assuming the calls on the arrays (1,2,3,0 aka AA,AB,BB,no call) had allele A as the reference allele and B as the variant. This is apparently
not the case. Don't make my mistake, only use the het calls to determine how well your variant calling pipeline is doing vs the arrays.

	2.) Parse SNP array to create VCFs for each sample column in the summary SNP array file. Can sometimes call python script from command line on cluster, but it'll run out
	of memory if too many people are running stuff on the login nodes.
	Bash script (parse_snp_array.sh):
	#!/bin/sh
	# give the job a name to help keep track of running jobs (optional)
	#PBS -N PARSE_SNP_ARRAYS
	#PBS -m e
	#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=16gb

	export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
	source activate anaconda

	python3 /scratch/jandrews/bin/parse_snp_array.py /scratch/jandrews/Data/Variant_Calling/SNP_Arrays/GenomeWideSNP_6.na35.annot.csv /scratch/jandrews/Data/Variant_Calling/SNP_Arrays/LYMPHOMASNP77_GTYPE_2014.txt


	3.) Remove any potential garbage chromosomes from both array VCFs and those from samtools/varscan.
	(sed '/_g/d' file.vcf | sed '/chrM/d' | sed '/chrY/d') > output.vcf

	or

	for file in *.vcf; do
		sed -i '/_g/d' "$file";
		sed -i '/chrM/d' "$file";
		sed -i '/chrY/d' "$file";
	done


	4A.) Figure out which variants on the SNP array contain sufficient read depth (>=5) for each sample to actually be called by samtools or varscan. 
	First, get read depth at each SNP array variant for each sample.
	Bash script (array_variant_cov.sh):
	#!/bin/sh

	# give the job a name to help keep track of running jobs (optional)
	#PBS -N SNP_ARRAY_COV
	#PBS -m e
	#PBS -l nodes=1:ppn=8,walltime=8:00:00,vmem=16gb

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


	4B.) Actually filter for read-depth. Change $11 as needed (ie. if read count column is column 9, change to $9)
	for file in *.vcf; do 
		base=${file%%_*} ; 
		cat $file | awk '$11 > 4' > "$base"_array_coverage_min5DP.vcf ; 
	done


	5.) Sort each filtered SNP array VCF.
	for f in *.vcf; do
    	base=${f%%.*} ;
    	vcf-sort "$f" > ${base}_sorted.vcf ;
	done


	6.) Cut last column (counts) from filtered SNP array VCFs.
	for f in *.vcf; do
		base=${f%%.*} ;
		cut -f1-10 "$f" > "$base"_cut.vcf ;
	done


	7.) Re-add header to each filtered SNP array, as bedtools multicov removes it.
	for f in *.vcf; do
		base=${f%%_*} ;
		sed -i '1s/^/#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	'$base'_ARRAY\n/' "$f" ;
		sed -i '1s/^/##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n/' "$f" ;
		sed -i '1s/^/##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n/' "$f" ;
		sed -i '1s/^/##reference=hg19\n/' "$f" ;
		sed -i '1s/^/##source=parse_snp_array.py\n/' "$f" ;
		sed -i '1s/^/##fileformat=VCFv4.2\n/' "$f" ;
	done


	8.) Zip each VCF.
	for f in *.vcf; do
		bgzip -c "$f" > "$f".gz ;
	done


	9. Index each VCF.
	for f in *.gz; do
		tabix -p vcf "$f" ;
	done


	10.) Create two directories - one for array/BCFtools intersections, one for array/VarScan intersections.
	Copy files into each as appropriate.


	11.) Intersect files for each sample (checks by position only - '-c all' option does this.).
	module load bcftools-1.2
	for f in *_RNAseq*.gz; do
		base=${f%%_*} ;
		bcftools isec -c all -p "$base" "$f" "$base"_array_coverage_min5DP_het.vcf.gz ;
	done
	module remove bcftools-1.2


	In this example, 0000.vcf will be records unique to first provided file. 0001.vcf will be records unique to second provided file. 
	0002.vcf and 0003.vcf are those shared between them.


	12.) Count lines for each file, ignoring the headers. Do whatever with that information. Venn diagrams or something.
	grep -v '#' -c file.vcf 



2.) Remove any potential garbage chromosomes from all VCFs.
for file in *.vcf; do
	base=${file%%_*} ;
	(sed '/_g/d' "$file" | sed '/chrM/d' | sed '/chrY/d') > "$base"_RNAseq_VS_clean.vcf ;
done


3.) Sort all VCFs.
for f in *.vcf; do
   	base=${f%%_*} ;
    vcf-sort "$f" > ${base}_RNAseq_VS_sorted.vcf ;
done


4.) Zip and index each VCF.
for f in *.vcf; do
	bgzip -c "$f" > "$f".gz;
done

for f in *.gz; do
	tabix -p vcf "$f" ;
done


5.) Merge VCFs from VarScan with each other. Then those from BCFTools with each other. '-m none' means multiallelic records will be split to separate lines. 
Doing so is rather important, as if two samples have a different variant alleles at the same position, only one is reported as having the variant if multiallelic
records are allowed. Alternatively, setting "-m both" should create a multiallelic record, which may be wanted at times. No idea what the default is, BCFtools docs don't mention.

For VarScan:
module load bcftools-1.2
bcftools merge -O v -m none --force-samples -i ADP:sum *.gz > merge.vcf
module remove bcftools-1.2

For BCFTools:
module load bcftools-1.2
bcftools merge -O v -m none -i DP:sum *.gz > merge.vcf
module remove bcftools-1.2


6.) Fix headers as necessary in text editor. 


7.) Concatenate the merged BCFTools and Varscan files. Be sure the column order is the same for both files.
For duplicates, it will print the record with more samples called. If the variant is called in the same number of samples between both files, the variant from the
BCFTools file will be printed. Adds an INFO field (BOTH) that specifies which file the variant was found in (BCF, VS, or BOTH).
Python script (cat_merged_coding_vcfs.py):
python3 cat_merged_coding_vcfs.py <merged_BCF.vcf> <merged_VarScan.vcf> <output.vcf>


8.) Sort output file.
vcf-sort coding_VS_BCF_final.vcf > coding_VS_BCF_final.sorted.vcf


9.) Annotate with VEP. This is essentially impossible to get working on the cluster due to how perl is set up on it, so install and run locally. 
Be sure to use the GrCH37 cache (--port 3337) for hg19, not GrCH38. Motif info is pulled from JASPAR mainly, it seems. 

perl ~/bin/ensembl-tools-release-82/scripts/variant_effect_predictor/variant_effect_predictor.pl --everything --vcf --format vcf --fork 2 --symbol --cache --port 3337 -i coding_VS_BCF_final.vcf -o coding_VS_BCF_final_annotated.vcf


10.) Intersect with GM TF ChIP-Seq data.
bedtools intersect -wa -wb -a /scratch/jandrews/Data/Variant_Calling/Coding/Final_Results/coding_VS_BCF_final_annotated.vcf -b GM12878_TF151_names_final.bed > Coding_variants_GM_ChIP_TFs_isec.bed


# annotate with VEP
 - assigns variants to known proteins
 - predicts impact

# remove known SNPs

# filter for those predicted to have moderate to high impact


###-To figure out insert sizes for PINDEL-###
Our samples have an average/median insert size of 200 bp, though the histograms look more like majority are 150 bp.

Bash script(get_insert_sizes.sh):
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
