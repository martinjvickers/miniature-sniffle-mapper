#!/bin/bash

nd_se_test () 
{
	if [ $1 -eq 0 ]; then
		echo "[PASS] "$2
	else
		echo "[FAIL] "$2
		exit 1
	fi
}


## full one when dealing with XD tags

## create an index
if [ ! -f example_data/TAIR10_chr_all.fa ]; then
	echo "[TEST 1/5] Extracting and creating reference index"
	cp bismark_example/TAIR10/TAIR10_chr_all.fa.gz example_data/
	gunzip example_data/TAIR10_chr_all.fa.gz
	nd_se_test $? "Reference copied and extracted"
	./build_reference -i example_data/TAIR10_chr_all.fa &> debug_index.txt
	nd_se_test $? "Reference index created";
else
	echo "[TEST 1/4] Reference index already exists"
fi

## Just getting the reads right for SE non-directional
#
#  Single End - Non-Directional
#
echo "[TEST 2/5] Running SE non-directional app test"

samtools view bismark_example/test_data_BSSeq_SE_nondirectional_100k_bismark_bismark_bt2.bam | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$10"\t"$14}' | sort -k1 > bismark_reads_ND_SE.dat
nd_se_test $? "Extract Non-directional SE Bismark Reads"
./miniature-sniffle-mapper -i example_data/test_data_BSSeq_SE_nondirectional_100k.fastq.gz -r example_data/TAIR10_chr_all.fa --ambig ignore --debug --non-directional -o test_ND_SE.bam &> debug_ND_SE.txt
nd_se_test $? "Run MSM on Non-directional SE"
samtools view test_ND_SE.bam| awk '{if(NF==20) print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$10"\t"$14; else print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$10"\t"$15 }' | sort -k1 > msm_reads_ND_SE.dat
nd_se_test $? "Extract Non-directional SE MSM Reads"
diff -q bismark_reads_ND_SE.dat msm_reads_ND_SE.dat
nd_se_test $? "MSM and Bismark Non-directional SE results"

## Just getting the reads right for SE directional

### to get the file 
## ./Bismark/bismark TAIR10 ../example_data/test_data_BSSeq_SE_directional_100k_bismark.fastq.gz

echo "[TEST 3/5] Running SE Directional app test"

## this is what the below should be when we get the XD tags working for directional SE
#
#  Single End - Directional
#
samtools view bismark_example/test_data_BSSeq_SE_directional_100k_bismark_bismark_bt2.bam | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$10"\t"$14}' | sort -k1 > bismark_reads_D_SE.dat
nd_se_test $? "Extract Directional SE Bismark Reads"
./miniature-sniffle-mapper -i example_data/test_data_BSSeq_SE_directional_100k.fastq.gz -r example_data/TAIR10_chr_all.fa --ambig ignore --debug -o test_D_SE.bam &> debug_D_SE.txt
nd_se_test $? "Run MSM on Directional SE"
samtools view test_D_SE.bam | awk '{if(NF==20) print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$10"\t"$14; else print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$10"\t"$15 }' | sort -k1 > msm_reads_D_SE.dat
nd_se_test $? "Extract Directional SE MSM Reads"
diff -q bismark_reads_D_SE.dat msm_reads_D_SE.dat
nd_se_test $? "MSM and Bismark Directional SE results"

echo "[TEST 4/5] Running SE non-directional all overlapping chromosome edge"

## Bismark BAM file created using;
## ./Bismark/bismark TAIR10 --non_directional ../example_data/test_data_BSSeq_SE_nondirectional_100k_alloverlapchr_bismark.fastq.gz

samtools view bismark_example/test_data_BSSeq_SE_nondirectional_100k_alloverlapchr_bismark_bismark_bt2.bam | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$10"\t"$14}' | sort -k1 > bismark_alloverlap_ND_SE.dat
nd_se_test $? "Extract Directional SE Bismark overlapping chromosome edge Reads"
./miniature-sniffle-mapper -i example_data/test_data_BSSeq_SE_nondirectional_100k_alloverlapchr.fastq.gz -r example_data/TAIR10_chr_all.fa --ambig ignore --non-directional --debug -o test_ND_SE_alloverlap.bam &> debug_ND_SE_alloverlap.txt
nd_se_test $? "Run MSM on Non-Directional SE overlapping chromosome edge "
samtools view test_ND_SE_alloverlap.bam| awk '{if(NF==20) print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$10"\t"$14; else print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$10"\t"$15 }' | sort -k1 > MSM_alloverlap_ND_SE.dat
nd_se_test $? "Extract Non-Directional SE MSM overlapping chromosome edge Reads"
diff -q bismark_alloverlap_ND_SE.dat MSM_alloverlap_ND_SE.dat
nd_se_test $? "MSM and Bismark Non-Directional SE overlapping chromosome edge results"

## let's set one up for paired end then
echo "[TEST 5/5] Running PE directional app test"
## Bismark BAM file created using;
## ./Bismark/bismark TAIR10 -1 ../example_data/test_data_BSSeq_PE_directional_1_100k.fastq.gz -2 ../example_data/test_data_BSSeq_PE_directional_2_100k.fastq.gz
#nd_se_test $? "Run MSM on Directional PE"
#
#   Paired End - Directional
#
./miniature-sniffle-mapper \
-1 example_data/test_data_BSSeq_PE_directional_1_100k.fastq.gz \
-2 example_data/test_data_BSSeq_PE_directional_2_100k.fastq.gz \
-r example_data/TAIR10_chr_all.fa \
--ambig ignore --debug \
-o test_D_PE_alloverlap.bam &> debug_D_PE_alloverlap.txt

nd_se_test $? "Run MSM on Directional PE"
samtools view test_D_PE_alloverlap.bam | awk '{if(NF==21) print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$10"\t"$14; else print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$10"\t"$15 }' | sort -k1 --version-sort > MSM_D_PE.dat
nd_se_test $? "Extract Directional PE MSM Reads"
samtools view bismark_example/test_data_BSSeq_PE_directional_1_100k_bismark_bt2_pe.bam | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$10"\t"$14}' | sort -k1 --version-sort > bismark_D_PE.dat
nd_se_test $? "Extract Directional PE Bismark Reads"


### For current debugging purposes
cat MSM_D_PE.dat | awk '{print $1}' | uniq | sort > MSM.txt
cat bismark_D_PE.dat | awk '{print $1}' | uniq | sort > BIS.txt
wc -l MSM.txt BIS.txt
comm -2 -3 MSM.txt BIS.txt | wc -l
comm -1 -3 MSM.txt BIS.txt | wc -l








## Cleanup
# rm debug_D_PE_alloverlap.txt test_D_PE_alloverlap.bam MSM_alloverlap_ND_SE.dat debug_ND_SE_alloverlap.txt test_ND_SE_alloverlap.bam bismark_alloverlap_ND_SE.dat msm_reads_D_SE.dat debug_D_SE.txt test_D_SE.bam bismark_reads_D_SE.dat msm_reads_ND_SE.dat debug_ND_SE.txt test_ND_SE.bam bismark_reads_ND_SE.dat BIS.txt MSM.txt bismark_D_PE.dat MSM_D_PE.dat

