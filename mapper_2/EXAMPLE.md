# EXAMPLE

This file is basically so that I can preserve the reads I'm interested in that exhibit interesting behaviour between bismark and miniature-sniffle-mapper when I switch between my work and home machines. As I get more aspects working this file will not be needed.


## non-direction

Bismark gets;


100000 reads; of these:
  100000 (100.00%) were unpaired; of these:
    85559 (85.56%) aligned 0 times
    11200 (11.20%) aligned exactly 1 time
    3241 (3.24%) aligned >1 times
14.44% overall alignment rate
100000 reads; of these:
  100000 (100.00%) were unpaired; of these:
    85361 (85.36%) aligned 0 times
    10899 (10.90%) aligned exactly 1 time
    3740 (3.74%) aligned >1 times
14.64% overall alignment rate
100000 reads; of these:
  100000 (100.00%) were unpaired; of these:
    89796 (89.80%) aligned 0 times
    7700 (7.70%) aligned exactly 1 time
    2504 (2.50%) aligned >1 times
10.20% overall alignment rate
100000 reads; of these:
  100000 (100.00%) were unpaired; of these:
    89860 (89.86%) aligned 0 times
    7764 (7.76%) aligned exactly 1 time
    2376 (2.38%) aligned >1 times
10.14% overall alignment rate
Processed 100000 sequences in total


I get


100000 reads; of these:
  100000 (100.00%) were unpaired; of these:
    99995 (100.00%) aligned 0 times
    2 (0.00%) aligned exactly 1 time
    3 (0.00%) aligned >1 times
0.01% overall alignment rate
100000 reads; of these:
  100000 (100.00%) were unpaired; of these:
    99995 (100.00%) aligned 0 times
    2 (0.00%) aligned exactly 1 time
    3 (0.00%) aligned >1 times
0.01% overall alignment rate
100000 reads; of these:
  100000 (100.00%) were unpaired; of these:
    85361 (85.36%) aligned 0 times
    10899 (10.90%) aligned exactly 1 time
    3740 (3.74%) aligned >1 times
14.64% overall alignment rate
100000 reads; of these:
  100000 (100.00%) were unpaired; of these:
    85559 (85.56%) aligned 0 times
    11200 (11.20%) aligned exactly 1 time
    3241 (3.24%) aligned >1 times
14.44% overall alignment rate


For direction I run;

Running: bowtie2 -q --score-min L,0,-0.2 --ignore-quals --nofw -p 1 -x example_data/TAIR10_chr_all.fa_G2A -U ref.temp.C2T.fastq.gz
Running: bowtie2 -q --score-min L,0,-0.2 --ignore-quals --norc -p 1 -x example_data/TAIR10_chr_all.fa_C2T -U ref.temp.C2T.fastq.gz

For non directional I run;

Running: bowtie2 -q --score-min L,0,-0.2 --ignore-quals --nofw -p 1 -x example_data/TAIR10_chr_all.fa_G2A -U ref.temp.G2A.fastq.gz
Running: bowtie2 -q --score-min L,0,-0.2 --ignore-quals --norc -p 1 -x example_data/TAIR10_chr_all.fa_C2T -U ref.temp.G2A.fastq.gz
Running: bowtie2 -q --score-min L,0,-0.2 --ignore-quals --nofw -p 1 -x example_data/TAIR10_chr_all.fa_G2A -U ref.temp.C2T.fastq.gz
Running: bowtie2 -q --score-min L,0,-0.2 --ignore-quals --norc -p 1 -x example_data/TAIR10_chr_all.fa_C2T -U ref.temp.C2T.fastq.gz

So the first two are failing.

	 bowtie2 -q --score-min L,0,-0.2 --ignore-quals --nofw -p 1 -x bismark_example/TAIR10/Bisulfite_Genome/CT_conversion/BS_CT -U ref.temp.G2A.fastq.gz > meh1

This works though

	 bowtie2 -q --score-min L,0,-0.2 --ignore-quals --nofw -p 1 -x example_data/TAIR10_chr_all.fa_C2T -U ref.temp.G2A.fastq.gz

	 bowtie2 -q --score-min L,0,-0.2 --ignore-quals --nofw -p 1 -x example_data/TAIR10_chr_all.fa_G2A -U ref.temp.G2A.fastq.gz
HUH!?!

      ### In the above example the number of transliterations required to transform the actual sequence
      ### to the C->T version would be TAGTTATGTGTGTGTG -> TAGTTATGTGTGTGTG = 0; (assuming this gives the correct alignment)
      ### in the G->A case it would be TAGTTATGTGTGTGTG -> TAATTATATATATATA = 6; (assuming this gives multiple wrong alignments)
      ### if the sequence giving a unique best alignment required a lower number of transliterations than the second best sequence yielding alignments
      ### while requiring a much higher number of transliterations, we are going to accept the unique best alignment with the lowest number of performed
      ### transliterations. As a threshold which does scale we will start with the number of tranliterations of the lowest best match x 2 must still be
      ### smaller than the number of tranliterations of the second best sequence. Everything will be flagged with $sequence_fails = 1 and discarded.
     
	if ($mismatches{$mismatch_number}->{$composite_location}->{index} == 0 or $mismatches{$mismatch_number}->{$composite_location}->{index} == 1){
	  $transliterations_performed = determine_number_of_transliterations_performed($sequence,'CT');
	}
	elsif ($mismatches{$mismatch_number}->{$composite_location}->{index} == 2 or $mismatches{$mismatch_number}->{$composite_location}->{index} == 3){
	  $transliterations_performed = determine_number_of_transliterations_performed($sequence,'GA');
	}




## Mismatches and transliterals

This is the example we are looking for

```
HS3:420:C3EHMACXX:8:1101:1936:2217 1:N:0:ACCTCA
Added 0 12814776 2 24 CT GA
Added 1 12857495 2 24 CT GA
Added 2 3056443 0 60 CT CT
Added 3 3056910 0 60 CT CT
```

The CT CT conversion has zero mismatches, but 60 transliterals. Meaning either of the others should be the best one. 


## Kicking things out

```
cat looking.txt | grep HS3:420:C3EHMACXX:8:1101:6433:2038
Adding HS3:420:C3EHMACXX:8:1101:6433:2038 1:N:0:ACCTTA 2 11724410 27 0 -7
Adding HS3:420:C3EHMACXX:8:1101:6433:2038 1:N:0:ACCTTA 4 3954339 27 1 -19
Adding HS3:420:C3EHMACXX:8:1101:6433:2038 1:N:0:ACCTTA 4 3624454 43 2 -19
Adding HS3:420:C3EHMACXX:8:1101:6433:2038 1:N:0:ACCTTA 4 3625340 43 3 -19
```

This read is one that mine keeps but bismark boots it out.

```
HS3:420:C3EHMACXX:8:1101:6433:2038	16	Chr5	11724411	24	100M	*	0	0	NGGATTTTGGTTATATTATGAAAGTTTTGAGAAGTAAGAAGAAGGTTGGTTAGTGTTTTGGTGTTGAATATGATTTGATGTTATGTGTATGATTGAGTATDCCAC@DCA;DFEDD@EEEEBFFCDFHHCHHHFHEIHGHAHIJIIGHHIHHHFHJIJIIHIHFGGFFHHGGCIJHHHJJIHJJJJIIFFFFFFFDDD=4#	AS:i:-7	XS:i:-19	XN:i:0	XM:Z:..................h.......h........z.............................h.................hh.h..h.....h....	XO:i:0	XG:Z:GA	NM:i:2	MD:Z:75T23A0	YT:Z:UU	XR:Z:CT	RF:Z:CGGACTTTGGCTACACCATGAAAGATTTGAGAAGCAAGAAGAAGGTTGGTTAGTGTTTTGGTGTCGAATATGACTTGATGTCATGTGTATGATTGAGTAT
```

Hmmm, bismark seems to be getting different information

```
    983 sequence NGGATTTTGGTTATATTATGAAAGTTTTGAGAAGTAAGAAGAAGGTTGGTTAGTGTTTTGGTGTTGAATATGATTTGATGTTATGTGTATGATTGAGTAT
    984 id HS3:420:C3EHMACXX:8:1101:6433:2038_1:N:0:ACCTTA
    985 quality: '#4=DDDFFFFFFFIIJJJJHIJJHHHJICGGHHFFGGFHIHIIJIJHFHHHIHHGIIJIHAHGHIEHFHHHCHHFDCFFBEEEE@DDEFD;ACD@CACCD'
    986 Index: 0
    987 HS3:420:C3EHMACXX:8:1101:6433:2038_1:N:0:ACCTTA 0       Chr3_CT_converted       13589300        0       100M    *       0       0       NGGATTTTGGTTATATTATGAAAGTTTTGAGAAGTAAGAAGAAGGTTGGTTAGTGTTTTGGTGTTGAAT    987 ATGATTTGATGTTATGTGTATGATTGAGTAT    #4=DDDFFFFFFFIIJJJJHIJJHHHJICGGHHFFGGFHIHIIJIJHFHHHIHHGIIJIHAHGHIEHFHHHCHHFDCFFBEEEE@DDEFD;ACD@CACCD    AS:i:-7 XS:i:-13        XN:i:0  XM:i:2  XO:i:0  XG:i:0  NM:i:2  MD    987 :Z:0T60A38    YT:Z:UU
    988 HS3:420:C3EHMACXX:8:1101:6433:2038_1:N:0:ACCTTA
    989 last seq id: HS3:420:C3EHMACXX:8:1101:6433:2038_1:N:0:ACCTTA and identifier: HS3:420:C3EHMACXX:8:1101:6433:2038_1:N:0:ACCTTA
    990 Index: 1
    991 HS3:420:C3EHMACXX:8:1101:6433:2038_1:N:0:ACCTTA 16      Chr2_GA_converted       3607693 0       100M    *       0       0       ATACTCAATCATACACATAACATCAAATCATATTCAACACCAAAACACTAACCAACCTTCTTCTTACTTCTCAAAAC    991 TTTCATAATATAACCAAAATCCN    DCCAC@DCA;DFEDD@EEEEBFFCDFHHCHHHFHEIHGHAHIJIIGHHIHHHFHJIJIIHIHFGGFFHHGGCIJHHHJJIHJJJJIIFFFFFFFDDD=4#    AS:i:-7 XS:i:-7 XN:i:0  XM:i:2  XO:i:0  XG:i:0  NM:i:2  MD:Z:38T60A0    YT    991 :Z:UU
    992 HS3:420:C3EHMACXX:8:1101:6433:2038_1:N:0:ACCTTA
    993 last seq id: HS3:420:C3EHMACXX:8:1101:6433:2038_1:N:0:ACCTTA and identifier: HS3:420:C3EHMACXX:8:1101:6433:2038_1:N:0:ACCTTA
    994 $alignment_ambiguous now: 1
    995 HS3:420:C3EHMACXX:8:1101:6433:2038_1:N:0:ACCTTA 256     *       0       0       *       *       0       0       NGGATTTTGGTTATATTATGAAAGTTTTGAGAAGTAAGAAGAAGGTTGGTTAGTGTTTTGGTGTTGAATATGATTTGATGTTATGTGTATGAT    995 TGAGTAT    #4=DDDFFFFFFFIIJJJJHIJJHHHJICGGHHFFGGFHIHIIJIJHFHHHIHHGIIJIHAHGHIEHFHHHCHHFDCFFBEEEE@DDEFD;ACD@CACCD
    996 sequence GGGAAGAAAATGTAGTTGTGAAAATATAATGGATAATTTTTTTATATTTTGATTATTATAGTAGTAGTGTTTGTGGTTAGATTATTTGGAGAAATTGATA
```

```
Mine
TTTTAAATTTTAAATTTTAAATTTTAAATTTTTGAATTTTTAATTTTTAAATTTTTAAATTTTTAAATTT
TATATTTATGAATTTTTAAATATTTAATTTTTTAAATTTGAAATTGGTTTTTTTGGTTGAAAATTATTGT

Bismark
TTTTAAATTTTAAATTTTAAATTTTAAATTTTTGAATTTTTAATTTTTAAATTTTTAAATTTTTAAATTT
TATATTTATGAATTTTTAAATATTTAATTTTTTAAATTTGAAATTGGTTTTTTTGGTTGAAAATTATTGT
```

```
GA_converted
CCCTAAACCCTAAACCCTAAACCCTAAACCTCTAAATCCTTAATCCCTAAATCCCTAAATCTTTAAATCCTACATCCAT
AAATCCCTAAATACCTAATTCCCTAAACCCAAAACCAATTTCTCTAATTAAAAATCATTATATATATAATAATAATTTT

>Chr1 
CCCTAAACCCTAAACCCTAAACCCTAAACCTCTAAATCCTTAATCCCTAAATCCCTAAATCTTTAAATCCTACATCCAT
AAATCCCTAAATACCTAATTCCCTAAACCCAAAACCAATTTCTCTAATTAAAAATCATTAT
```

Solved, the issue here was that I had the -k 2 option set which is only used in bowtie v1.

### New thing

I get the following, bismark doesnot get it

```
HS3:420:C3EHMACXX:8:1101:5227:2035      0       Chr2    1518    42      100M    *       0       0       NAAAAAAAATGTTGTTGAATAATTTTTGAAAATTATTGGATATGATGTAATGTTTTGTGATTGAATTTTTTAAAATATATTAATAAAGAGTTTAGGATGT    #4=AAAAAAB3?A+@B):ABBBBBBBB0=ABBBBBBBB9*==?A4>)=>7>ABBBB26)=>A)=AAAAA?=@@@@;??>?>BBBB??#############    AS:i:-1 XS:i:-1 XN:i:0  XM:Z:..............xz......h...z......h.............h.............z....h.h........h..h...................       XO:i:0  XG:Z:CT NM:i:1  MD:Z:0A99       YT:Z:UU XR:Z:CT RF:Z:AAAAAAAAATGTTGCCGAATAACTTTCGAAAATCATTGGATATGATGCAATGTTTTGTGATCGAATCTCTTAAAATACATCAATAAAGAGTTTAGGATGT
```


```
   1128 id HS3:420:C3EHMACXX:8:1101:5227:2035_1:Y:0:ACCTCA
   1129 quality: '#4=AAAAAAB3?A+@B):ABBBBBBBB0=ABBBBBBBB9*==?A4>)=>7>ABBBB26)=>A)=AAAAA?=@@@@;??>?>BBBB??#############'
   1130 HASH(0x1194ae8)
   1131 Index: 0
   1132 HS3:420:C3EHMACXX:8:1101:5227:2035_1:Y:0:ACCTCA 0       Chr2_CT_converted       1518    1       100M    *       0       0       NAAAAAAAATGTTGTTGAATAATTTTTGAAAATTATTGGATATGATGTAATGTTTTGTGATTGAATTTTTTAAAATA
   1132 TATTAATAAAGAGTTTAGGATGT    #4=AAAAAAB3?A+@B):ABBBBBBBB0=ABBBBBBBB9*==?A4>)=>7>ABBBB26)=>A)=AAAAA?=@@@@;??>?>BBBB??#############    AS:i:-1 XS:i:-1 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:0A99       YT
   1132 :Z:UU
   1133 HS3:420:C3EHMACXX:8:1101:5227:2035_1:Y:0:ACCTCA
   1134 last seq id: HS3:420:C3EHMACXX:8:1101:5227:2035_1:Y:0:ACCTCA and identifier: HS3:420:C3EHMACXX:8:1101:5227:2035_1:Y:0:ACCTCA
   1135 First alignment score, setting $best_AS_so_far to -1
   1136 First  best alignment_score is: '-1'
   1137 MD tag is: '0A99'
   1138 second best alignment_score is: '-1'
   1139 
   1140 Read is ambiguous within the same thread, or otherwise as good as the best one so far. Setting $amb_same_thread to 1 for currently best AS: -1
   1141 Index: 0        The current Seq-ID is HS3:420:C3EHMACXX:8:1101:5227:2035_1:Y:0:ACCTCA, skipped all ambiguous sequences until the next ID which is: HS3:420:C3EHMACXX:8:1101:5078:2129_1:N:0:ACCTCA
   1142 HASH(0x11babd8)
   1143 Index: 1
   1144 HS3:420:C3EHMACXX:8:1101:5227:2035_1:Y:0:ACCTCA 4       *       0       0       *       *       0       0       NAAAAAAAATGTTGTTGAATAATTTTTGAAAATTATTGGATATGATGTAATGTTTTGTGATTGAATTTTTTAAAATATATTAATAAAGAGTTT
   1144 AGGATGT    #4=AAAAAAB3?A+@B):ABBBBBBBB0=ABBBBBBBB9*==?A4>)=>7>ABBBB26)=>A)=AAAAA?=@@@@;??>?>BBBB??#############    YT:Z:UU
   1145 HS3:420:C3EHMACXX:8:1101:5227:2035_1:Y:0:ACCTCA
   1146 last seq id: HS3:420:C3EHMACXX:8:1101:5227:2035_1:Y:0:ACCTCA and identifier: HS3:420:C3EHMACXX:8:1101:5227:2035_1:Y:0:ACCTCA
   1147 $alignment_ambiguous now: 1
   1148 HS3:420:C3EHMACXX:8:1101:5227:2035_1:Y:0:ACCTCA 256     *       0       0       *       *       0       0       NAAAAAAAATGTTGTTGAATAATTTTTGAAAATTATTGGATATGATGTAATGTTTTGTGATTGAATTTTTTAAAATATATTAATAAAGAGTTT
   1148 AGGATGT    #4=AAAAAAB3?A+@B):ABBBBBBBB0=ABBBBBBBB9*==?A4>)=>7>ABBBB26)=>A)=AAAAA?=@@@@;??>?>BBBB??#############
```

Running through bowtie2 on there own.

Hmmm, but for bowtie2 on its own I get;


```
mvickers@n108379:~/development/miniature-sniffle-mapper/bismark_example_2$ bowtie2 -q --score-min L,0,-0.2 --ignore-quals --norc -x TAIR10/Bisulfite_Genome/CT_conversion/BS_CT -U test2_CT.fastq1 reads; of these:
  1 (100.00%) were unpaired; of these:
    0 (0.00%) aligned 0 times
    0 (0.00%) aligned exactly 1 time
    1 (100.00%) aligned >1 times
100.00% overall alignment rate
@HD	VN:1.0	SO:unsorted
@SQ	SN:Chr1_CT_converted	LN:30427671
@SQ	SN:Chr2_CT_converted	LN:19698289
@SQ	SN:Chr3_CT_converted	LN:23459830
@SQ	SN:Chr4_CT_converted	LN:18585056
@SQ	SN:Chr5_CT_converted	LN:26975502
@SQ	SN:ChrC_CT_converted	LN:154478
@SQ	SN:ChrM_CT_converted	LN:366924
@PG	ID:bowtie2	PN:bowtie2	VN:2.1.0
HS3:420:C3EHMACXX:8:1101:5227:2035_1:Y:0:ACCTCA	0	Chr2_CT_converted	1518	1	100M	*	0	0	NAAAAAAAATGTTGTTGAATAATTTTTGAAAATTATTGGATATGATGTAATGTTTTGTGATTGAATTTTTTAAAATATATTAATAAAGAGTTTAGGATGT	#4=AAAAAAB3?A+@B):ABBBBBBBB0=ABBBBBBBB9*==?A4>)=>7>ABBBB26)=>A)=AAAAA?=@@@@;??>?>BBBB??#############	AS:i:-1	XS:i:-1	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:0A99	YT:Z:UU
mvickers@n108379:~/development/miniature-sniffle-mapper/bismark_example_2$ bowtie2 -q --score-min L,0,-0.2 --ignore-quals --nofw -x TAIR10/Bisulfite_Genome/GA_conversion/BS_GA -U test2_CT.fastq
1 reads; of these:
  1 (100.00%) were unpaired; of these:
    1 (100.00%) aligned 0 times
    0 (0.00%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
0.00% overall alignment rate
@HD	VN:1.0	SO:unsorted
@SQ	SN:Chr1_GA_converted	LN:30427671
@SQ	SN:Chr2_GA_converted	LN:19698289
@SQ	SN:Chr3_GA_converted	LN:23459830
@SQ	SN:Chr4_GA_converted	LN:18585056
@SQ	SN:Chr5_GA_converted	LN:26975502
@SQ	SN:ChrC_GA_converted	LN:154478
@SQ	SN:ChrM_GA_converted	LN:366924
@PG	ID:bowtie2	PN:bowtie2	VN:2.1.0
HS3:420:C3EHMACXX:8:1101:5227:2035_1:Y:0:ACCTCA	4	*	0	0	*	*	0	0	NAAAAAAAATGTTGTTGAATAATTTTTGAAAATTATTGGATATGATGTAATGTTTTGTGATTGAATTTTTTAAAATATATTAATAAAGAGTTTAGGATGT#4=AAAAAAB3?A+@B):ABBBBBBBB0=ABBBBBBBB9*==?A4>)=>7>ABBBB26)=>A)=AAAAA?=@@@@;??>?>BBBB??#############	YT:Z:UU
```


CORRECT

ST-E00144:285:HVT3GCCXX:1:1101:20567:1678_2	16	3	9081562	40	100M	*	0	0	TTTTACAAGCTAATCACAAACATAATAATATAATTCTCAAAATTCAACAAATATATATATAAATTCTAATTTAAAAAAAACAAAAAAAACCATATAATCAAAJFJJA7-JAAFAAJJJFAFFJAJJJJFFJFJ-AJJJJJJJJFJJJJJJJJJJJJJFJFFFJJFJJJJJJJJJJJJJJJJJJJJFFFFJJJJJFJFFJF	AS:i:-6	XN:i:0	XM:Z:....h..x....h..................hh......x........z.h.h...h...................h....zx...h........h....	XO:i:0	XG:Z:GA	NM:i:1	MD:Z:8T91	YT:Z:UU	XR:Z:CT	MV:Z:TTTTACAAACTAATCACAAACATAATAATATAATTCTCAAAATTCAACAAATATATATATAAATTCTAATTTAAAAAAAACAAAAAAAACCATATAATCA

ST-E00144:285:HVT3GCCXX:1:1101:20567:1678_1:N:0:ATCACGTT	16	3	9081562	40	100M	*	0	0	TTTTACAAGCTAATCACAAACATAATAATATAATTCTCAAAATTCAACAAATATATATATAAATTCTAATTTAAAAAAAACAAAAAAAACCATATAATCA	AAJFJJA7-JAAFAAJJJFAFFJAJJJJFFJFJ-AJJJJJJJJFJJJJJJJJJJJJJFJFFFJJFJJJJJJJJJJJJJJJJJJJJFFFFJJJJJFJFFJF	NM:i:16	MD:Z:4G2G0T3G18G0G6G8G1G1G3G19G4G0G3G8G4	XM:Z:....h..x....h..................hh......x........z.h.h...h...................h....zx...h........h....	XR:Z:CT	XG:Z:GA

XM:Z:....h..x....h..................hh......x........z.h.h...h...................h....zx...h........h....
XM:Z:....h..x....h..................hh......x........z.h.h...h...................h....zx...h........h....

XR:Z:CT 
XG:Z:GA

TTTTACAAGCTAATCACAAACATAATAATATAATTCTCAAAATTCAACAAATATATATATAAATTCTAATTTAAAAAAAACAAAAAAAACCATATAATCA
TTTTACAAGCTAATCACAAACATAATAATATAATTCTCAAAATTCAACAAATATATATATAAATTCTAATTTAAAAAAAACAAAAAAAACCATATAATCA
....h..x....h..................hh......x........z.h.h...h...................h....zx...h........h....
ACTAATATACCAAAAAAAACAAAAAAAATTTAATCTTAAATATATATATAAACAACTTAAAACTCTTAATATAATAATACAAACACTAATCGAACATTTT
TTTTACAAACTAATCACAAACATAATAATATAATTCTCAAAATTCAACAAATATATATATAAATTCTAATTTAAAAAAAACAAAAAAAACCATATAATCA
....h..x....h..................hh......x........z.h.h...h...................h....zx...h........h....

NOT CORRECT

ST-E00144:285:HVT3GCCXX:1:1101:21745:1713_23    0       1       1191063 42      100M    *       0       0       ATTAATCAACATTAAATCCAACAAAACAAAAATTAACCTATATTTCATTAATCTATAAAACAACCAAATACTATACAAAAATCTAACTAAATTCAACAAAJJJJJJJJJJJFJJJJJJJJJJJJJJ<JJJJJFFJJJJJJJJJJFJJJFFJAJJJJJJJJJJJJJJJJJJJFJFJJJJJJJFJFJJJJJJJJJJJJJJJJ        AS:i:0  XN:i:0  XM:Z:......X..Z.......HH..Z....H.........HH.......H......H.......H..HX.....H....X......X...H......X..H...       XO:i:0  XG:Z:GA NM:i:0  MD:Z:100        YT:Z:UU XR:Z:GA MV:Z:ATTAATCAACATTAAATCCAACAAAACAAAAATTAACCTATATTTCATTAATCTATAAAACAACCAAATACTATACAAAAATCTAACTAAATTCAACAAA       SQ:Z:GTTGGTCAGCGTTGAATCCAACGGAACAAGAGTTGACCTATATTTCATTGATCTATGGAACAACCAGATACTATGCAGAGATCTGACTAGATTCAGCAAG

ST-E00144:285:HVT3GCCXX:1:1101:21745:1713_1:N:0:ATCACGTT	16	1	1191063	42	100M	*	0	0	ATTAATCAACATTAAATCCAACAAAACAAAAATTAACCTATATTTCATTAATCTATAAAACAACCAAATACTATACAAAAATCTAACTAAATTCAACAAA	JJJJJJJJJJJFJJJJJJJJJJJJJJ<JJJJJFFJJJJJJJJJJFJJJFFJAJJJJJJJJJJJJJJJJJJJFJFJJJJJJJFJFJJJJJJJJJJJJJJJJ	NM:i:22	MD:Z:0G2G0G3G1G2G8G0G5G1G2G14G6G0G8G7G2G1G4G4G5G3G0	XM:Z:z..hh...x.z..h........zx.....h.h..h..............h......hh........x.......h..x.h....x....h.....x...h	XR:Z:GA	XG:Z:GA

XG:Z:GA
XR:Z:GA

REF	GTTGGTCAGCGTTGAATCCAACGGAACAAGAGTTGACCTATATTTCATTGATCTATGGAACAACCAGATACTATGCAGAGATCTGACTAGATTCAGCAAG
QRY	ATTAATCAACATTAAATCCAACAAAACAAAAATTAACCTATATTTCATTAATCTATAAAACAACCAAATACTATACAAAAATCTAACTAAATTCAACAAA
MET	z..hh...x.z..h........zx.....h.h..h..............h......hh........x.......h..x.h....x....h.....x...h
NEW	h..hx...z.h..h........hh.....h.h..x..............h......hh........h.......z..h.h....x....h.....z....
	z..hh...x.z..h........zx.....h.h..h..............h......hh........x.......h..x.h....x....h.....x...h

	........z...........z......h.h........................h.................h....h...................h.z
	z.h...................h....h.................h........................h.h......z...........z........

	z..hh...x.z..h........zx.....h.h..h..............h......hh........x.......h..x.h....x....h.....x...h

NOT CORRECT

MINE

ST-E00144:285:HVT3GCCXX:1:1101:19948:1696_14	16	1	8501397	42	100M	*	0	0	TAGAATTTTTTAGTTGGAGATTTAATTTTAATTTTTATATGAAGAATATTATTTAATTTTATTTGTTTTTAATTAAAAAAAAAATTTAATTATTTTGTATFJFFAJJFJJFFJFJJF<JJFFJFJJJFJFJJJJFJJFJFJJJJJJJFAJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJF	AS:i:0	XN:i:0	XM:Z:..X.........X..ZX.H.....................H..H....................H...............................H...	XO:i:0	XG:Z:CT	NM:i:0	MD:Z:100	YT:Z:UU	XR:Z:GA	MV:Z:TAGAATTTTTTAGTTGGAGATTTAATTTTAATTTTTATATGAAGAATATTATTTAATTTTATTTGTTTTTAATTAAAAAAAAAATTTAATTATTTTGTAT

BISMARK

ST-E00144:285:HVT3GCCXX:1:1101:19948:1696_1:N:0:ATCACGTT	0	1	8501397	42	100M	*	0	0	TAGAATTTTTTAGTTGGAGATTTAATTTTAATTTTTATATGAAGAATATTATTTAATTTTATTTGTTTTTAATTAAAAAAAAAATTTAATTATTTTGTAT	FJFFAJJFJJFFJFJJF<JJFFJFJJJFJFJJJJFJJFJFJJJJJJJFAJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJF	NM:i:26	MD:Z:0C4C4C2C0C5C4C0C1C2C0C1C0C1C8C1C0C1C1C2C0C0C0C8C15C1C13	XM:Z:x....h....x..xz.....h....hh.h..hh.hh.h........h.hh.h.h..hhhh........h...............h.h.............	XR:Z:GA	XG:Z:CT

TAGAATTTTTTAGTTGGAGATTTAATTTTAATTTTTATATGAAGAATATTATTTAATTTTATTTGTTTTTAATTAAAAAAAAAATTTAATTATTTTGTAT
TAGAATTTTTTAGTTGGAGATTTAATTTTAATTTTTATATGAAGAATATTATTTAATTTTATTTGTTTTTAATTAAAAAAAAAATTTAATTATTTTGTAT
..X.........X..ZX.H.....................H..H....................H...............................H...
x....h....x..xz.....h....hh.h..hh.hh.h........h.hh.h.h..hhhh........h...............h.h.............

XR:Z:GA
XG:Z:CT

TAGAATTTTTTAGTTGGAGATTTAATTTTAATTTTTATATGAAGAATATTATTTAATTTTATTTGTTTTTAATTAAAAAAAAAATTTAATTATTTTGTAT





ST-E00144:285:HVT3GCCXX:1:1101:20567:1678_3	16	3	9081562	40	100M	*	0	0	TTTTACAAGCTAATCACAAACATAATAATATAATTCTCAAAATTCAACAAATATATATATAAATTCTAATTTAAAAAAAACAAAAAAAACCATATAATCAAAJFJJA7-JAAFAAJJJFAFFJAJJJJFFJFJ-AJJJJJJJJFJJJJJJJJJJJJJFJFFFJJFJJJJJJJJJJJJJJJJJJJJFFFFJJJJJFJFFJF	NM:i:16	MD:Z:4G2G0T3G18G0G6G8G1G1G3G19G4G0G3G8G4	XM:Z:....h..x....h..................hh......x........z.h.h...h...................h....zx...h........h....	XR:Z:CT	XG:Z:GA


....h..x....h..................hh......x........z.h.h...h...................h....zx...h........h....
....h..x....h..................hh......x........z.h.h...h...................h....zx...h........h....

....h..x....h..................hh......x........z.h.h...h...................h....zx...h........h....




## Meh


#### Insertion and deletion issue



71M1D4M





Parsing CIGAR
NS500422:333:H3TJVBGXY:1:11101:10255:2527_1523                         *
TGACTCCTGGTAAGGAACCGTCTTTACATCTCTAACTTGTTGGCTCGTTAGGAAGGGAAGTTAAAGCTTGCCTGTAAGG         orig
TGACTCCTGGTAAGGAACCGTCTTTACATCTCTAACTTGTTGGCTCGTTAGGAAGGGAAGTTAAAGCTTGCTGTAAGG          corrected
  ATTTTTGGTAAGGAATTGTTTTTATATTTTTAATTTGTTGGTTTGTTAGGAAGGGAAGTTAAAGTTTGTTTTAAG           read
  .h.hx..........xz..h....h..h.h...h.......h.z....................h...h......           mine

  .h.hx..........xz..h....h..h.h...h.......h.z....................h...hh.....           bismark




a




ANOTHER

<-----------

Parsing CIGAR
NS500422:333:H3TJVBGXY:1:11101:10016:14444_14207
                                                        *
AGGTTTCCAATGAGTAAAGATCATCCCCAACCACTGCGATGAGCGGTGGAGATTGGACAGAAGAGAAAGACATGTCT   orig
AGGTTTCCAATGAGTAAAGATCATCCCCAACCACTGCGATGAGCGGTGGAGATTGGCAGAAGAGAAAGACATGTCT    corrected
AAATTTCCAATAAATAAAAATCATCCCCAACCACTACAATAAACAATAAAAATTAAACAAAAAAAAAAACATATT     read
.hh........h.h....h................x.z..h.h.zx.hh.h...hh..z..h.h...h....x..     mine

.hh........h.h....h................x.z..h.h.zx.hh.h...hh..z..h.h...h....h..     bismark

< NS500422:333:H3TJVBGXY:1:11101:10016:14444    16      3       10009985        58M1D17M        AAATTTCCAATAAATAAAAATCATCCCCAACCACTACAATAAACAATAAAAATTAAACAAAAAAAAAAACATATT     XM:Z:.hh........h.h....h................x.z..h.h.zx.hh.h...hh..z..h.h...h....h..
---
> NS500422:333:H3TJVBGXY:1:11101:10016:14444    16      3       10009985        58M1D17M        AAATTTCCAATAAATAAAAATCATCCCCAACCACTACAATAAACAATAAAAATTAAACAAAAAAAAAAACATATT     XM:Z:.hh........h.h....h................x.z..h.h.zx.hh.h...hh..z..h.h...h....x..





