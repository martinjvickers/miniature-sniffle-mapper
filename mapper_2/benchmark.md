# Benchmarking bismark and miniature-sniffle-mapper

Results and information regarding performance.

## Setup

```
cd /dev/shm
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX218/SRX2182198/SRR4280515/SRR4280515.sra

wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar xvf sratoolkit.current-ubuntu64.tar.gz
./sratoolkit.2.8.2-1-ubuntu64/bin/fastq-dump --gzip --split-files -v SRR4280515.sra

git clone https://github.com/FelixKrueger/Bismark.git

git clone https://github.com/seqan/seqan.git

sudo apt-get install bowtie2 samtools zlib1g-dev g++ cmake libtbb-dev libboost-dev

cd miniature-sniffle-mapper/bismark_example
gunzip TAIR10/TAIR10_chr_all.fa.gz
/usr/bin/time  -v ../../Bismark/bismark -p 4 --multicore 8 TAIR10 ../../SRR4280515_1.fastq.gz
```

## Results

## CPU AWS 244GB RAM 32 core instance

processor	: 31
vendor_id	: GenuineIntel
cpu family	: 6
model		: 79
model name	: Intel(R) Xeon(R) CPU E5-2686 v4 @ 2.30GHz
stepping	: 1
microcode	: 0xb00001b
cpu MHz		: 1391.769
cache size	: 46080 KB
physical id	: 0
siblings	: 32
core id		: 15
cpu cores	: 16
apicid		: 31
initial apicid	: 31
fpu		: yes
fpu_exception	: yes
cpuid level	: 13
wp		: yes
flags		: fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush mmx fxsr sse sse2 ht syscall nx pdpe1gb rdtscp lm constant_tsc rep_good nopl xtopology nonstop_tsc aperfmperf eagerfpu pni pclmulqdq monitor est ssse3 fma cx16 pcid sse4_1 sse4_2 x2apic movbe popcnt tsc_deadline_timer aes xsave avx f16c rdrand hypervisor lahf_lm abm 3dnowprefetch fsgsbase bmi1 hle avx2 smep bmi2 erms invpcid rtm rdseed adx xsaveopt ida
bugs		:
bogomips	: 4600.08
clflush size	: 64
cache_alignment	: 64
address sizes	: 46 b

## Bismark bowtie2 4 core and split 8 times

	Command being timed: "../../Bismark/bismark -p 4 --multicore 8 TAIR10 ../../SRR4280515_1.fastq.gz"
	User time (seconds): 10479.66
	System time (seconds): 26427.22
	Percent of CPU this job got: 2666%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 23:03.86
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 445056
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 1
	Minor (reclaiming a frame) page faults: 4097893
	Voluntary context switches: 2325300
	Involuntary context switches: 8209813547
	Swaps: 0
	File system inputs: 160
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0

## Reducing the number of multicores (4 this time)

	Command being timed: "../../Bismark/bismark -p 4 --multicore 4 TAIR10 ../../SRR4280515_1.fastq.gz"
	User time (seconds): 8658.06
	System time (seconds): 19475.55
	Percent of CPU this job got: 2253%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 20:48.52
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 460588
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 3896986
	Voluntary context switches: 1929285
	Involuntary context switches: 63094403
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0

## This time it beats mine :-(

	Command being timed: "../../Bismark/bismark -p 2 --multicore 16 TAIR10 ../../SRR4280515_1.fastq.gz"
	User time (seconds): 8678.24
	System time (seconds): 8301.16
	Percent of CPU this job got: 2296%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 12:19.29
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 263072
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 4254012
	Voluntary context switches: 2958745
	Involuntary context switches: 540152623
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0

## Mine using fastq.gz output

	Command being timed: "./miniature-sniffle-mapper -i SRR4280515_1.fastq.gz -r TAIR10_chr_all.fa"
	User time (seconds): 5919.84
	System time (seconds): 17148.87
	Percent of CPU this job got: 1960%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 19:36.90
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 4021548
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 1686807
	Voluntary context switches: 2284050
	Involuntary context switches: 37614716
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0

## This time I will stop using compressed temp file. This really increases speed

	Command being timed: "./miniature-sniffle-mapper -i SRR4280515_1.fastq.gz -r TAIR10_chr_all.fa"
	User time (seconds): 5544.46
	System time (seconds): 16875.76
	Percent of CPU this job got: 2799%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 13:20.82
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 4013860
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 2096547
	Voluntary context switches: 2314546
	Involuntary context switches: 26791436
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0

## I tried changing the buffer size, didn't really do anything

	Command being timed: "./miniature-sniffle-mapper -i SRR4280515_1.fastq.gz -r TAIR10_chr_all.fa"
	User time (seconds): 5382.06
	System time (seconds): 16580.39
	Percent of CPU this job got: 2792%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 13:06.43
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 4013420
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 2108831
	Voluntary context switches: 2301334
	Involuntary context switches: 25955646
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0

