#!/bin/sh

## Run MSM methylation extract
/usr/bin/time -v ./extract_methylation -CX -i MSM_D_SE.bam -o MSM_D_SE -r example_data/TAIR10_chr4.fa
/usr/bin/time -v ./Bismark/bismark_methylation_extractor --comprehensive --cytosine_report --bedGraph --CX BSSeq_SE_D_mod_bismark_bt2.bam --genome_folder TAIR10_CHR4 >> bismark_methyl_extract_SE_D.log 2>&1
diff -y MSM_D_SE_CX_report.txt BSSeq_SE_D_mod_bismark_bt2.CX_report.txt --suppress-common-lines

/usr/bin/time -v ./extract_methylation -CX -i MSM_ND_SE.bam -o MSM_ND_SE -r example_data/TAIR10_chr4.fa
/usr/bin/time -v ./Bismark/bismark_methylation_extractor --comprehensive --cytosine_report --bedGraph --CX BSSeq_SE_ND_mod_bismark_bt2.bam --genome_folder TAIR10_CHR4 >> bismark_methyl_extract_SE_ND.log 2>&1
diff -y MSM_ND_SE_CX_report.txt BSSeq_SE_ND_mod_bismark_bt2.CX_report.txt --suppress-common-lines


#cat test_data_BSSeq_SE_nondirectional_100k_bismark_bismark_bt2.CX_report.txt.CX_report.txt | awk '{if($1==3) print $0}' | awk '{if($4!=0 || $5!=0) print $2"\t"$4"\t"$5"\t"$6}' | sort > bismark_CX_chr1.txt

#./extract_methylation -i bismark_example/test_data_BSSeq_SE_nondirectional_100k_bismark_bismark_bt2.bam -o test -r example_data/TAIR10_chr_all.fa > test_CX.txt

#cat test_CX.txt | awk '{if($1==3) print $0}' | awk '{if($4!=0 || $5!=0) print $2"\t"$3"\t"$4"\t"$5}' | sort > test_CX_chr1.txt

#wc -l bismark_CX_chr1.txt test_CX_chr1.txt

#comm -2 -3 bismark_CX_chr1.txt test_CX_chr1.txt > only_in_bismark.txt
#comm -1 -3 bismark_CX_chr1.txt test_CX_chr1.txt > only_in_mine.txt
#wc -l only_in_bismark.txt only_in_mine.txt
#diff bismark_CX_chr1.txt test_CX_chr1.txt 
