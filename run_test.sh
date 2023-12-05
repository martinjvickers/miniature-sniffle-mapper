#!/bin/sh

### SE directional test

export PATH=~/bin:$PATH
export PATH

echo "STARTING: SE Directional Test"

cd bismark_example
rm -rf Bismark
git clone https://github.com/FelixKrueger/Bismark.git
find . -name "*.gz"| awk '{print "gunzip "$1}'| sh
/usr/bin/time -v Bismark/bismark TAIR10 ../example_data/test_data_BSSeq_SE_directional_100k_bismark.fastq.gz
cd ..

samtools view bismark_example/test_data_BSSeq_SE_directional_100k_bismark_bismark_bt2.bam | awk '{print $1}'| sort > bismark-res.txt

/usr/bin/time -v ./miniature-sniffle-mapper --debug -i example_data/test_data_BSSeq_SE_directional_100k.fastq.gz -r example_data/TAIR10_chr_all.fa > looking-se-d.txt
samtools view results.bam | awk '{print $1}'| sort > results-res.txt
mv results.bam msm-SE-results.bam

echo "Reads only visible by Bismark"
comm -2 -3 bismark-res.txt results-res.txt | wc -l
echo "Reads only visible by mine"
comm -1 -3 bismark-res.txt results-res.txt | wc -l
echo "How many reads mapped? These numbers should be the same"
wc -l bismark-res.txt results-res.txt
#diff bismark-res.txt results-res.txt

echo "COMPLETE: SE Directional Test"

### SE Non-directional test
echo "STARTING: SE non-directional Test"

cd bismark_example
/usr/bin/time -v Bismark/bismark --non_directional TAIR10 ../example_data/test_data_BSSeq_SE_nondirectional_100k_bismark.fastq.gz
cd ..

samtools view bismark_example/test_data_BSSeq_SE_nondirectional_100k_bismark_bismark_bt2.bam | awk '{print $1}'| sort > bismark-nd-res.txt

/usr/bin/time -v ./miniature-sniffle-mapper --debug --non-directional -i example_data/test_data_BSSeq_SE_nondirectional_100k.fastq.gz -r example_data/TAIR10_chr_all.fa > looking-se-nd.txt
samtools view results.bam | awk '{print $1}'| sort > results-nd-res.txt
mv results.bam msm-SE-nd-results.bam

echo "Reads only visible by Bismark"
comm -2 -3 bismark-nd-res.txt results-nd-res.txt | wc -l
echo "Reads only visible by mine"
comm -1 -3 bismark-nd-res.txt results-nd-res.txt | wc -l
echo "How many reads mapped? These numbers should be the same"
wc -l bismark-nd-res.txt results-nd-res.txt

echo "COMPLETE: SE non-directional Test"

echo "STARTING: PE Directional Test"

cd bismark_example
/usr/bin/time -v Bismark/bismark TAIR10 -1 ../example_data/test_data_BSSeq_PE_directional_1_100k.fastq.gz -2 ../example_data/test_data_BSSeq_PE_directional_2_100k.fastq.gz
cd ..

/usr/bin/time -v ./miniature-sniffle-mapper --debug -1 example_data/test_data_BSSeq_PE_directional_1_100k.fastq.gz -2 example_data/test_data_BSSeq_PE_directional_2_100k.fastq.gz -r example_data/TAIR10_chr_all.fa > looking-pe.txt

echo "COMPLETE: PE directional Test"

cd example_data
rm *C2T.fastq *C2T.fastq.gz *G2A.fastq.gz *G2A.fastq
cd ../
