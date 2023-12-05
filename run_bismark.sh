#!/bin/sh

## pass fail function
nd_se_test () 
{
        if [ $1 -eq 0 ]; then
                echo "[PASS] "$2
        else
                echo "[FAIL] "$2
                exit 1
        fi
}

## get hold of the latest Bismark
git clone https://github.com/FelixKrueger/Bismark.git

## extract the genome reference
cp -r bismark_example/TAIR10 .
gunzip ./TAIR10/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa.gz ./TAIR10/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa.gz ./TAIR10/TAIR10_chr_all.fa.gz

## run bismark
##zcat example_data/test_data_BSSeq_SE_nondirectional_100k.fastq.gz | awk '{if(NR%4==1) print $0"_"((NR-1)/4)+1; else print $0}' | gzip > example_data/test_data_BSSeq_SE_nondirectional_100k_bismark.fastq.gz
##Bismark/bismark --non_directional TAIR10 example_data/test_data_BSSeq_SE_nondirectional_100k_bismark.fastq.gz

zcat example_data/test_data_BSSeq_SE_directional_100k.fastq.gz | awk '{if(NR%4==1) print $0"_"((NR-1)/4); else print $0}' | gzip > example_data/test_data_BSSeq_SE_directional_100k_bismark.fastq.gz
Bismark/bismark TAIR10 example_data/test_data_BSSeq_SE_directional_100k_bismark.fastq.gz

## Create MSM index
#if [ ! -f example_data/TAIR10_chr_all.fa ]; then
#        echo "[TEST 1/5] Extracting and creating reference index"
#        cp bismark_example/TAIR10/TAIR10_chr_all.fa.gz example_data/
#        gunzip example_data/TAIR10_chr_all.fa.gz
#        nd_se_test $? "Reference copied and extracted"
#        ./build_reference -i example_data/TAIR10_chr_all.fa >> debug_index.txt 2>&1
#        nd_se_test $? "Reference index created";
#else
#        echo "[TEST 1/4] Reference index already exists"
#fi


## this is what the below should be when we get the XD tags working for directional SE
samtools view test_data_BSSeq_SE_directional_100k_bismark_bismark_bt2.bam | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$10"\t"$14}' | sort -k1 > bismark_reads_D_SE.dat
nd_se_test $? "Extract Directional SE Bismark Reads"
./miniature-sniffle-mapper -i example_data/test_data_BSSeq_SE_directional_100k.fastq.gz -r example_data/TAIR10_chr_all.fa --ambig ignore --debug -o test_D_SE.bam >> debug_D_SE.txt 2>&1
nd_se_test $? "Run MSM on Directional SE"
samtools view test_D_SE.bam | awk '{if(NF==20) print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$10"\t"$14; else print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$10"\t"$15 }' | sort -k1 > msm_reads_D_SE.dat
nd_se_test $? "Extract Directional SE MSM Reads"
diff -q bismark_reads_D_SE.dat msm_reads_D_SE.dat
nd_se_test $? "MSM and Bismark Directional SE results"

echo "[TEST 4/5] Running SE non-directional all overlapping chromosome edge"


### to run CX report

