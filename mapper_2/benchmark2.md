## First real benchmark

Now that the results are the same between MSM and Bismark, we will do an initial benchmark. This will be followed by a real benchmark script to run the setup on AWS with a variety of parameters but for now this is exploritory.

* NOTE1: At this point I am not clearing caches between runs. I will get to that at some point for the official tests.
* NOTE2: The bowtie2 part of the test needs a little thinking about I think since only one of the two bowtie2 runs for a `directional` test

## Test Hardware Setup

Date: 18-10-2017

AWS 8 v-core machine (m4.2xlarge) with 64GB general purpose SSD.

```
vendor_id       : GenuineIntel
cpu family      : 6
model           : 79
model name      : Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz
```

## Software

Bismark version 	v0.19.0
`f20b2ecd3473e5181839f7c1c44d2d21361d595d`

MSM			alpha-0.0.4
`596577065dd5e6c4a33cfafe083c2fd87547d14c`

```
NAME="Ubuntu"
VERSION="14.04.5 LTS, Trusty Tahr"
ID=ubuntu
ID_LIKE=debian
PRETTY_NAME="Ubuntu 14.04.5 LTS"
VERSION_ID="14.04"
```

### Software install and compile

```
sudo apt-get install build-essential
sudo apt-get install -y g++-4.9 cmake cmake-data zlib1g-dev libbz2-dev libboost-dev bowtie2 libtbb-dev samtools git
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
export CXX="g++-4.9" CC="gcc-4.9"
git clone https://github.com/martinjvickers/miniature-sniffle-mapper.git
git clone https://github.com/seqan/seqan.git seqan
cd miniature-sniffle-mapper/
cmake ../miniature-sniffle-mapper -DCMAKE_MODULE_PATH=../seqan/util/cmake -DSEQAN_INCLUDE_PATH=../seqan/include -DCMAKE_CXX_FLAGS=-std=c++14 -DCMAKE_BUILD_TYPE=Release
```

## Test

These are the test commands.

### Getting the data

```
~/sratoolkit.2.8.2-1-ubuntu64/bin/prefetch SRR3301576
~/sratoolkit.2.8.2-1-ubuntu64/bin/fastq-dump --gzip SRR3301576.sra
```

### Running MSM

```
/usr/bin/time -v ./miniature-sniffle-mapper -i SRR3301576.fastq.gz -o SRR3301576.bam -r TAIR10_chr_all.fa --ambig ignore -c 4
```

### Running Bismark

```
$ /usr/bin/time -v Bismark/bismark TAIR10 ../SRR3301576.fastq.gz -p 4
```

### Running Bowtie2

```
/usr/bin/time -v bowtie2 -q --score-min L,0,-0.2 --ignore-quals --nofw -p 4 -x TAIR10_chr_all.fa_G2A -U tmp/SRR3301576.fastq.gz_ref.temp.C2T.fastq > /dev/null
```

## Results

|   | Bismark  | MSM  | Bowtie2 |
|---|----------|------|---------|
| 1 |  14:06   | 4:37 | 2:29    |
| 2 |  14:03   | 4:22 | 2:20    |
| 3 |  13:57   | 4:28 | 1:56    | <--- probably should clear caches TBH
