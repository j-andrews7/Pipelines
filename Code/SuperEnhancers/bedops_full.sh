#!/bin/sh
# give the job a name to help keep track of running jobs (optional)
#PBS -N BEDOPS_FULL
#PBS -m e
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=8gb
export PATH=/act/Anaconda3-2.3.0/bin:${PATH}
source activate anaconda

for fold in /scratch/jandrews/Data/ChIP_Seq/ROSE/ROSE_SEs_From_Ind_Sample_Peaks/SE_INTERSECTS/BEDOPS_Overlap_Method/SELECT_SAMPLES/; do

    cd "$fold"
    bedops --everything *.bed \
    | sort-bed - > all_files.bed ;
    bedmap --echo-map-range --ec --sweep-all --fraction-either 0.25 all_files.bed \
    | sort-bed - \
    | bedmap --echo-map-range --ec --sweep-all --fraction-either 0.25 - \
    | sort-bed - \
    | bedmap --echo-map-range --ec --sweep-all --fraction-either 0.25 - \
    | sort-bed - \
    | bedmap --echo-map-range --ec --sweep-all --fraction-either 0.25 - \
    | sort-bed - \
    | bedmap --echo-map-range --ec --sweep-all --fraction-either 0.25 - \
    | sort-bed - \
    | bedmap --echo-map-range --ec --sweep-all --fraction-either 0.25 - \
    | sort-bed - \
    | bedmap --echo-map-range --ec --sweep-all --fraction-either 0.25 - \
    | sort-bed - \
    | bedmap --echo-map-range --ec --sweep-all --fraction-either 0.25 - \
    | sort-bed - \
    | bedmap --echo-map-range --ec --sweep-all --fraction-either 0.25 - \
    | sort-bed - \
    | uniq - > ranges.bed
    bedmap --echo-map-id-uniq --ec --sweep-all --fraction-either 0.25 ranges.bed all_files.bed > samples.bed ;
    paste ranges.bed samples.bed > ./Results/All_SEs.bed
    python /scratch/jandrews/bin/get_recurrent_SEs.py ./Results/All_SEs.bed ./Results/Recurrent_All_SEs.bed ";" 

    sed '/CC\|CB\|VG\|TS/!d' ./Results/All_SEs.bed > ./Results/NORMAL_SEs.bed
    sed '/FL\|DL\|CLL/d' ./Results/All_SEs.bed > ./Results/NORMAL_ONLY_SEs.bed
    python /scratch/jandrews/bin/get_recurrent_SEs.py ./Results/NORMAL_ONLY_SEs.bed ./Results/Recurrent_NORMAL_ONLY_SEs.bed ";"
    sed '/FL\|DL\|CLL/!d' ./Results/All_SEs.bed > ./Results/TUMOR_SEs.bed
    sed '/CC\|CB\|VG\|TS/d' ./Results/All_SEs.bed > ./Results/TUMOR_ONLY_SEs.bed
    python /scratch/jandrews/bin/get_recurrent_SEs.py ./Results/TUMOR_ONLY_SEs.bed ./Results/Recurrent_TUMOR_ONLY_SEs.bed ";"

    grep DL ./Results/All_SEs.bed > ./Results/DL_SEs.bed
    sed '/FL\|CC\|CB\|TS\|CLL\|VG/d' ./Results/DL_SEs.bed > ./Results/Unique_DL_SEs.bed
    python /scratch/jandrews/bin/get_recurrent_SEs.py ./Results/Unique_DL_SEs.bed ./Results/Recurrent_Unique_DL_SEs.bed ";" 
    sed '/TS\|CC\|CB\|CLL\|VG/d' ./Results/DL_SEs.bed | sed '/FL/!d' - > ./Results/Unique_DL_FL_SEs.bed
    sed '/TS\|FL\|CLL\|VG/d' ./Results/DL_SEs.bed | sed '/CC\|CB/!d' - > ./Results/Unique_DL_CC_CB_SEs.bed
    sed '/FL\|CC\|CB\|CLL\|MEM\|VG/d' ./Results/DL_SEs.bed | sed '/NAIVE/!d' - > ./Results/Unique_DL_NAIVE_SEs.bed
    sed '/FL\|CC\|CB\|CLL\|NAIVE\|VG/d' ./Results/DL_SEs.bed | sed '/MEM/!d' - > ./Results/Unique_DL_MEM_SEs.bed
    sed '/FL\|CC\|CB\|CLL\|TS\|VGR/d' ./Results/DL_SEs.bed | sed '/VGA/!d' - > ./Results/Unique_DL_VGA_SEs.bed
    sed '/FL\|CC\|CB\|CLL\|TS\|VGA/d' ./Results/DL_SEs.bed | sed '/VGR/!d' - > ./Results/Unique_DL_VGR_SEs.bed
    sed '/FL\|CC\|CB\|VG\|TS/d' ./Results/DL_SEs.bed | sed '/CLL/!d' - > ./Results/Unique_DL_CLL_SEs.bed


    grep FL ./Results/All_SEs.bed > ./Results/FL_SEs.bed
    sed '/DL\|CC\|CB\|TS\|CLL\|VG/d' ./Results/FL_SEs.bed > ./Results/Unique_FL_SEs.bed
    python /scratch/jandrews/bin/get_recurrent_SEs.py ./Results/Unique_FL_SEs.bed ./Results/Recurrent_Unique_FL_SEs.bed ";"
    sed '/TS\|DL\|CLL\|VG/d' ./Results/FL_SEs.bed | sed '/CC\|CB/!d' - > ./Results/Unique_FL_CC_CB_SEs.bed
    sed '/DL\|CC\|CB\|CLL\|MEM\|VG/d' ./Results/FL_SEs.bed | sed '/NAIVE/!d' - > ./Results/Unique_FL_NAIVE_SEs.bed
    sed '/DL\|CC\|CB\|CLL\|NAIVE\|VG/d' ./Results/FL_SEs.bed | sed '/MEM/!d' - > ./Results/Unique_FL_MEM_SEs.bed
    sed '/DL\|CC\|CB\|CLL\|TS\|VGR/d' ./Results/FL_SEs.bed | sed '/VGA/!d' - > ./Results/Unique_FL_VGA_SEs.bed
    sed '/DL\|CC\|CB\|CLL\|TS\|VGA/d' ./Results/FL_SEs.bed | sed '/VGR/!d' - > ./Results/Unique_FL_VGR_SEs.bed
    sed '/DL\|CC\|CB\|VG\|TS/d' ./Results/FL_SEs.bed | sed '/CLL/!d' - > ./Results/Unique_FL_CLL_SEs.bed


    grep NAIVE ./Results/All_SEs.bed > ./Results/NAIVE_SEs.bed
    sed '/DL\|CC\|CB\|MEM\|CLL\|VG\|FL/d' ./Results/NAIVE_SEs.bed > ./Results/Unique_NAIVE_SEs.bed
    python /scratch/jandrews/bin/get_recurrent_SEs.py ./Results/Unique_NAIVE_SEs.bed ./Results/Recurrent_Unique_NAIVE_SEs.bed ";"
    sed '/FL\|DL\|CLL\|VG\|MEM/d' ./Results/NAIVE_SEs.bed | sed '/CC\|CB/!d' - > ./Results/Unique_NAIVE_CC_CB_SEs.bed
    sed '/DL\|CC\|CB\|CLL\|FL\|VG/d' ./Results/NAIVE_SEs.bed | sed '/MEM/!d' - > ./Results/Unique_NAIVE_MEM_SEs.bed
    sed '/DL\|CC\|CB\|CLL\|MEM\|VGR\|FL/d' ./Results/NAIVE_SEs.bed | sed '/VGA/!d' - > ./Results/Unique_NAIVE_VGA_SEs.bed
    sed '/DL\|CC\|CB\|CLL\|MEM\|VGA\|FL/d' ./Results/NAIVE_SEs.bed | sed '/VGR/!d' - > ./Results/Unique_NAIVE_VGR_SEs.bed
    sed '/DL\|CC\|CB\|VG\|MEM\|FL/d' ./Results/NAIVE_SEs.bed | sed '/CLL/!d' - > ./Results/Unique_NAIVE_CLL_SEs.bed


    grep MEM ./Results/All_SEs.bed > ./Results/MEM_SEs.bed
    sed '/DL\|CC\|CB\|NAIVE\|CLL\|VG\|FL/d' ./Results/MEM_SEs.bed > ./Results/Unique_MEM_SEs.bed
    python /scratch/jandrews/bin/get_recurrent_SEs.py ./Results/Unique_MEM_SEs.bed ./Results/Recurrent_Unique_MEM_SEs.bed ";"
    sed '/FL\|DL\|CLL\|VG\|NAIVE/d' ./Results/MEM_SEs.bed | sed '/CC\|CB/!d' - > ./Results/Unique_MEM_CC_CB_SEs.bed
    sed '/DL\|CC\|CB\|CLL\|NAIVE\|VGR\|FL/d' ./Results/MEM_SEs.bed | sed '/VGA/!d' - > ./Results/Unique_MEM_VGA_SEs.bed
    sed '/DL\|CC\|CB\|CLL\|NAIVE\|VGA\|FL/d' ./Results/MEM_SEs.bed | sed '/VGR/!d' - > ./Results/Unique_MEM_VGR_SEs.bed
    sed '/DL\|CC\|CB\|VG\|NAIVE\|FL/d' ./Results/MEM_SEs.bed | sed '/CLL/!d' - > ./Results/Unique_MEM_CLL_SEs.bed


    grep VGA ./Results/All_SEs.bed > ./Results/VGA_SEs.bed
    sed '/DL\|CC\|CB\|TS\|CLL\|VGR\|FL/d' ./Results/VGA_SEs.bed > ./Results/Unique_VGA_SEs.bed
    python /scratch/jandrews/bin/get_recurrent_SEs.py ./Results/Unique_VGA_SEs.bed ./Results/Recurrent_Unique_VGA_SEs.bed ";"
    sed '/FL\|DL\|CLL\|VGR\|TS/d' ./Results/VGA_SEs.bed | sed '/CC\|CB/!d' - > ./Results/Unique_VGA_CC_CB_SEs.bed
    sed '/DL\|CC\|CB\|CLL\|TS\|FL/d' ./Results/VGA_SEs.bed | sed '/VGR/!d' - > ./Results/Unique_VGA_VGR_SEs.bed
    sed '/DL\|CC\|CB\|VGR\|TS\|FL/d' ./Results/VGA_SEs.bed | sed '/CLL/!d' - > ./Results/Unique_VGA_CLL_SEs.bed


    grep VGR ./Results/All_SEs.bed > ./Results/VGR_SEs.bed
    sed '/DL\|CC\|CB\|TS\|CLL\|VGA\|FL/d' ./Results/VGR_SEs.bed > ./Results/Unique_VGR_SEs.bed
    python /scratch/jandrews/bin/get_recurrent_SEs.py ./Results/Unique_VGR_SEs.bed ./Results/Recurrent_Unique_VGR_SEs.bed ";"
    sed '/FL\|DL\|CLL\|VGA\|TS/d' ./Results/VGR_SEs.bed | sed '/CC\|CB/!d' - > ./Results/Unique_VGR_CC_CB_SEs.bed
    sed '/DL\|CC\|CB\|VGA\|TS\|FL/d' ./Results/VGR_SEs.bed | sed '/CLL/!d' - > ./Results/Unique_VGR_CLL_SEs.bed


    grep CLL ./Results/All_SEs.bed > ./Results/CLL_SEs.bed
    sed '/DL\|CC\|CB\|TS\|VG\|FL/d' ./Results/CLL_SEs.bed > ./Results/Unique_CLL_SEs.bed
    python /scratch/jandrews/bin/get_recurrent_SEs.py ./Results/Unique_CLL_SEs.bed ./Results/Recurrent_Unique_CLL_SEs.bed ";"
    sed '/FL\|DL\|VG\|TS/d' ./Results/CLL_SEs.bed | sed '/CC\|CB/!d' - > ./Results/Unique_CLL_CC_CB_SEs.bed


    grep CC ./Results/All_SEs.bed > ./Results/CCCB_SEs.bed
    grep CB ./Results/All_SEs.bed >> ./Results/CCCB_SEs.bed
    sort-bed ./Results/CCCB_SEs.bed | uniq - > ./Results/CC_CB_SEs.bed
    sed '/FL\|TS\|DL\|VG\|CLL/d' ./Results/CC_CB_SEs.bed > ./Results/Unique_CC_CB_SEs.bed
    python /scratch/jandrews/bin/get_recurrent_SEs.py ./Results/Unique_CC_CB_SEs.bed ./Results/Recurrent_Unique_CC_CB_SEs.bed ";"

    rm ranges.bed
    rm samples.bed
    rm all_files.bed

done
wait