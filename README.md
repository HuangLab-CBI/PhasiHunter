# phasiHunter
welcome to phasiHunter 😉
A multithreaded program for mining phasiRNA regulation pathway based on multiple reference sequence

![](https://sandbox-1314381151.cos.ap-nanjing.myqcloud.com/pic/202303200044049.png)

# Table of contents
- Dependencies 
- Installation
- Usage

# Dependencies
phasiHunter is a CLI program runing on Linux platform. The correction runing of phasiHunter depends on some existing software.
- Bowtie
- Biopython
- bedtools
- dnapi
- trim_galore
- seqkit

# Installation
1. clone phasiHunter
`git clone https://github.com/HuangLab-CBI/phasiHunter.git .`

2. setting enviroment variable in ~/.bashrc
`export PATH=$PATH:<phasiHunter path>`

now type `phasiHunter -h` to check phasiHunter whether installation correct

# Usage
type `phasiHunter -h` for Usage information
![](https://sandbox-1314381151.cos.ap-nanjing.myqcloud.com/pic/202303200117824.png)

parament in < > means necessary; parament in [ ] means optional
1. Data pre-process
```bash
phasiHunter preprocess -m r -i <SRR5049781.fastq.gz> -r <oryza_sativa_cdna.fa> -o <SRR5049781_cdna.map>

phasiHunter preprocess -m r -i <SRR5049781.fastq.gz> -r <oryza_sativa_gdna.fa> -o <SRR5049781_gdna.map>
```
2. phasiRNA and PHAS Loci prediction
```bash 
phasiHunter phase -cm <SRR5049781_cdna.map> -c <oryza_sativa_cdna.fa> -gm <SRR5049781_gdna.map> -g <oryza_sativa_gdna.fa> -fa <SRR7851621_trimmed_format_filter.fa> -a [SRR5049781_allsiRNA.txt] -o [SRR5049781_phasiRNA.txt] -pl [21] -j [10] -pv [0.0001] -ps [15] -pr [0.4] 
```
3. phasiRNA and PHAS Loci result integration
```bash
phasiHunter integration -io <SRR5049781_phasiRNA.txt> -ia <SRR5049781_allsiRNA.txt> -an <oryza_sativa_gdna.gff3> -o [SRR5049781_phasiRNA_dup.txt] -a [SRR5049781_allsiRNA_dup.txt] -s [SRR5049781_summary.txt] 
```
4. print phasiRNA_cluster plot, phasiRNA.fa, PHAS.fa
```bash
phasiHunter visulization -io <SRR5049781_phasiRNA_dup.txt> -ia <SRR5049781_allsiRNA_dup.txt> -a [SRR5049781_alignment.txt] -o [SRR5049781.phasiRNA.fa] -p [SRR5049781.PHAS.fa] -c [oryza_sativa_cdna.ga] -g [oryza_sativa_gdna.fa] -pc y -pg n
```
5. initiator prediction and verification
```bash 
phasiHunter target -q <osa_miRNA.fa> -t <SRR5049781_PHAS.fa> -o <SRR5049781_miR.txt>

phasiHunter initiator -i <SRR5049781_phasiRNA_dup.txt> -j <SRR5049781_miR.txt> -o <SRR5049781_initiator.txt> -a <SRR5049781_initiator.fa>

phasiHunter target -q <SRR5049781_initiator.fa> -t <SRR5049781_PHAS.fa> -o <SRR5049781_initiator_target.txt>

phasiHunter deg -i <degradome.map> -q <osa_miRNA.fa> -j <SRR5049781_initiator_target.txt> -t <oryza_sativa_cdna.fa> -o [SRR5049781_initiator_verified.txt]
```

6. phasiRNA target prediction and verification
```bash
phasiHunter target -q <SRR5049781_phasiRNA.fa> -t <oryza_sativa_cdna.fa> -o <SRR5049781_phasiRNA_target.txt>

phasiHunter deg -i <degradome.map> -q <SRR5049781_phasiRNA.fa> -j <SRR5049781_phasiRNA_target.txt> -t <oryza_sativa_cdna.fa> -o [SRR5049781_phasiRNA_target_verified.txt]
```
