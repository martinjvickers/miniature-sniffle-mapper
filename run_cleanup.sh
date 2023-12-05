#!/bin/sh

## Clean up indexes
gzip example_data/TAIR10_chr4.fa
rm example_data/TAIR10_chr4.fa_*
rm -rf TAIR10_CHR4

## Remove Bismark
rm -rf Bismark

## Remove index logs
rm debug_index.txt

## Remove modified input fastq files
rm example_data/BSSeq_SE_D_mod.fq.gz
rm example_data/BSSeq_SE_ND_mod.fq.gz

## remove BAMs, mapped logs and dat files
rm BSSeq_SE_D_mod_bismark_bt2.bam
rm BSSeq_SE_ND_mod_bismark_bt2.bam
rm msm_reads_D_SE.dat
rm msm_reads_ND_SE.dat
rm bismark_reads_D_SE.dat
rm bismark_reads_ND_SE.dat
rm BSSeq_SE_D_mod_bismark_bt2_SE_report.txt
rm BSSeq_SE_ND_mod_bismark_bt2_SE_report.txt
rm MSM_D_SE.bam
rm MSM_ND_SE.bam
