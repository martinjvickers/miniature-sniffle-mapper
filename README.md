# miniature-sniffle-mapper
BS Mapper with a name suggested by github

## Install from source

### Dependencies for Ubuntu 

```
apt-get install cmake g++-4.9 cmake-data zlib1g-dev libbz2-dev libboost-dev libtbb-dev bowtie2
export CXX="g++-4.9" CC="gcc-4.9"
git clone https://github.com/seqan/seqan.git
git clone https://github.com/martinjvickers/miniature-sniffle-mapper.git
```

### Compile

```
cd miniature-sniffle-mapper
cmake ../miniature-sniffle-mapper -DCMAKE_MODULE_PATH=../seqan/util/cmake -DSEQAN_INCLUDE_PATH=../seqan/include -DCMAKE_CXX_FLAGS=-std=c++14 -DCMAKE_BUILD_TYPE=Release
```


#### My machine

```
~/cmake-3.6.0-rc2-Linux-x86_64/bin/cmake ../miniature-sniffle-mapper -DCMAKE_MODULE_PATH=../seqan/util/cmake -DSEQAN_INCLUDE_PATH=../seqan/include -DCMAKE_CXX_FLAGS=-std=c++14 -DCMAKE_BUILD_TYPE=Release
```

## Usage

This will be run on our example data.

### Build the reference

```
build_reference -i example_data/TAIR10_reference_all.fa
```

### Run Mapper

```
miniature-sniffle-mapper -i example_data/test_data_BSSeq_SE_directional_100k.fastq.gz -r example_data/TAIR10_reference_all.fa -o example_data/test_data_BSSeq_SE_directional_100k.bam -c 1
```

and if you wish to keep and randomly assign the ambigious reads, then add the `--ambig` flag.

### My example (later to be used for testing)

My version treats fastq IDs without mergins and stuff. Bowtie2 uses the readname as part of its random seed. So when comparing my results to bismark when considering the following read.

```
zcat example_data/rdr2_sperm_BSseq_rep1_L008_R1_100k.fastq.gz | head -4
@HS3:420:C3EHMACXX:8:1101:1230:2142 1:N:0:ACCTCA
AGTAAGAATTTAGAGAGTGATAGAGAATTTAATTGTGTATTATAGGTGGATAAAGTTATTGTGGTATTTTTTATGAGATTTAAAGATTTCGTAGATACGA
+
@@@DDFFDFH?DHIIIIB9EFACF@BFGGCFH<H9EFIH9FG>F*C1?D=?GGCFDFAGGFFCG=BBHEH);DD>@DCEHFGEECCDDEF;AC=@DDC1=

```

So the way to deal with this in order to compare what is happening, I'm changing OUR input for now to make them match 100%

```
zcat example_data/rdr2_sperm_BSseq_rep1_L008_R1_100k.fastq.gz | awk '{print $1}' | gzip > modified_10k.fastq.gz
```

This will allow for the bowie2 results generated will be identical between bismark/miniature-sniffle-mapper. 

## To make a file that bismark with bowtie will parse the same way as mine would

```
zcat test_data_BSSeq_SE_directional_100k.fastq.gz | awk '{if(NR%4==1) print $0"_"((NR-1)/4)+1; else print $0}' | gzip > test_data_BSSeq_SE_directional_100k-bismark.fastq.gz
```



```
cat BAXF003A_all_sub_mod.fastq | awk 'BEGIN{sum=0}{if(NR % 4 == 1){ print $0"_"sum; sum=sum+1; }else{ print $0}}' > BAXF003A_all_sub_mod_meh.fastq
```

Bismark version 0.19.0 and Bowtie2 2.2.6 was used to create the bismark test results so bowtie2 2.2.6 should be used to run the test harness for comparison.


## Methylation Simulator

Using the methylation simulator. NOTE: this will need updating with some detail

```
./methylation_simulator -i example_data/TAIR10_chr_all.fa -o meh.fa -r meh.fastq.gz -n 10
./miniature-sniffle-mapper -i meh.fastq.gz -r example_data/TAIR10_chr_all.fa --ambig ignore -o meh_D_SE_uniq.bam
~/software/samtools-1.5/samtools sort meh_D_SE_uniq.bam -o meh_D_SE_uniq_sorted.bam
~/software/samtools-1.5/samtools index meh_D_SE_uniq_sorted.bam
~/bin/bismark_methylation_extractor --comprehensive --cytosine_report --bedGraph --single-end --CX meh_D_SE_uniq_sorted.bam --genome_folder bismark_example/TAIR10
```

Odd bit with extracting methylation.

```
mvickers@n108379:~/development/miniature-sniffle-mapper$ diff huh_final.txt 70SC_1_val_1.CX_report.txt.CX_report_chr1.txt |head
127,128c127,128
< 1	287	-	11	44	CHH	CAC
< 1	288	-	6	48	CHH	CCA
---
> 1	287	-	11	43	CHH	CAC
> 1	288	-	6	49	CHH	CCA
788,789c788,789
< 1	1806	-	417	3482	CHH	CAT
< 1	1807	-	386	3545	CHH	CCA
---


```
