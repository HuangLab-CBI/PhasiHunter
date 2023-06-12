# phasihunter
welcome to phasihunter 😉
a multithreaded program for mining phasirna regulation pathway based on multiple reference sequence

![](https://sandbox-1314381151.cos.ap-nanjing.myqcloud.com/pic/202306122137212.png)

# table of contents
- dependencies 
- installation
- use case
- copyright

# dependencies
phasihunter is a cli program runing on linux platform. the correction runing of phasihunter depends on some existing software.
- bowtie
- biopython
- bedtools
- dnapi
- trim_galore
- seqkit

# installation
1. clone phasihunter
`git clone https://github.com/huanglab-cbi/phasihunter.git .`

2. setting enviroment variable in ~/.bashrc
`export path=$path:<phasihunter path>`

now type `phasihunter -h` to check phasihunter whether installation correct

# use case
parament in < > means necessary; parament in [ ] means optional

1. Data pre-process
```bash
phasiHunter preprocess -m <r> -i <SRR5049781.fastq.gz> -r <oryza_sativa_cdna.fa> -o [SRR5049781_cdna.map]

phasiHunter preprocess -m <r> -i <SRR5049781.fastq.gz> -r <oryza_sativa_gdna.fa> -o [SRR5049781_gdna.map]
```
2. phasiRNA and PHAS Loci prediction
```bash 
phasiHunter phase -cm <SRR5049781_cdna.map> -c <oryza_sativa_cdna.fa> -gm <SRR5049781_gdna.map> -g <oryza_sativa_gdna.fa> -fa <SRR7851621_trimmed_format_filter.fa> -a [SRR5049781_allsiRNA.txt] -o [SRR5049781_phasiRNA.txt] -pl [21] -j [10] -pv [0.0001] -ps [15] -pr [0.4] 
```
3. phasiRNA and PHAS Loci result integration
```bash
phasiHunter integration -io <SRR5049781_phasiRNA.txt> -ia <SRR5049781_allsiRNA.txt> -an <oryza_sativa_gdna.gff3> -o [SRR5049781_phasiRNA_dup.txt] -a [SRR5049781_allsiRNA_dup.txt] -s [SRR5049781_summary.txt] -po [SRR5049781_phas.txt] -g <y>
```
4. print phasiRNA_cluster plot, phasiRNA.fa, PHAS.fa
```bash
phasiHunter visulization -io <SRR5049781_phasiRNA_dup.txt> -ia <SRR5049781_allsiRNA_dup.txt> -ip <SRR5049781_phas.txt> -a [SRR5049781_alignment.txt] -o [SRR5049781.phasiRNA.fa] -p [SRR5049781.PHAS.fa] -c [oryza_sativa_cdna.fa] -g [oryza_sativa_gdna.fa] -pc [y] -pg [y]
```
5. initiator prediction and verification
```bash 
phasiHunter target -q <osa_miRNA.fa> -t <SRR5049781_PHAS.fa> -o <SRR5049781_miR.txt>

phasiHunter initiator -i <SRR5049781_phasiRNA_dup.txt> -j <SRR5049781_miR.txt> -ip <SRR5049781_phas.txt> -o <SRR5049781_initiator.txt>

phasiHunter deg -i <degradome_PHAS.map> -q <osa_miRNA.fa> -j <SRR5049781_initiator.txt> -t <SRR5049781_PHAS.fa> -o <SRR5049781_initiator_verified.txt> -in <y>
```

6. phasiRNA target prediction and verification
```bash
phasiHunter target -q <SRR5049781_phasiRNA.fa> -t <oryza_sativa_cdna.fa> -o <SRR5049781_phasiRNA_target.txt>

phasiHunter deg -i <degradome_cdna.map> -q <SRR5049781_phasiRNA.fa> -j <SRR5049781_phasiRNA_target.txt> -t <oryza_sativa_cdna.fa> -o <SRR5049781_phasiRNA_target_verified.txt> -in <n>
```

# copyright
copyright © crop bioinformatics group (cbi), college of agricultural, nanjing agricultural university
free for academic use. for commercial use, please contact us (huangji@njau.edu.cn)