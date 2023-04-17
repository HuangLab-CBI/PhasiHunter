# phasihunter
welcome to phasihunter 😉
a multithreaded program for mining phasirna regulation pathway based on multiple reference sequence

![](https://sandbox-1314381151.cos.ap-nanjing.myqcloud.com/pic/202303200044049.png)

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
1. data pre-process
```bash
phasihunter preprocess -m r -i <srr5049781.fastq.gz> -r <oryza_sativa_cdna.fa> -o <srr5049781_cdna.map>

phasihunter preprocess -m r -i <srr5049781.fastq.gz> -r <oryza_sativa_gdna.fa> -o <srr5049781_gdna.map>
```
2. phasirna and phas loci prediction
```bash 
phasihunter phase -cm <srr5049781_cdna.map> -c <oryza_sativa_cdna.fa> -gm <srr5049781_gdna.map> -g <oryza_sativa_gdna.fa> -fa <srr7851621_trimmed_format_filter.fa> -a [srr5049781_allsirna.txt] -o [srr5049781_phasirna.txt] -pl [21] -j [10] -pv [0.0001] -ps [15] -pr [0.4] 
```
3. phasirna and phas loci result integration
```bash
phasihunter integration -io <srr5049781_phasirna.txt> -ia <srr5049781_allsirna.txt> -an <oryza_sativa_gdna.gff3> -o [srr5049781_phasirna_dup.txt] -a [srr5049781_allsirna_dup.txt] -s [srr5049781_summary.txt] 
```
4. print phasirna_cluster plot, phasirna.fa, phas.fa
```bash
phasihunter visulization -io <srr5049781_phasirna_dup.txt> -ia <srr5049781_allsirna_dup.txt> -a [srr5049781_alignment.txt] -o [srr5049781.phasirna.fa] -p [srr5049781.phas.fa] -c [oryza_sativa_cdna.ga] -g [oryza_sativa_gdna.fa] -pc y -pg n
```
5. initiator prediction and verification
```bash 
phasihunter target -q <osa_mirna.fa> -t <srr5049781_phas.fa> -o <srr5049781_mir.txt>

phasihunter initiator -i <srr5049781_phasirna_dup.txt> -j <srr5049781_mir.txt> -o <srr5049781_initiator.txt> -a <srr5049781_initiator.fa>

phasihunter target -q <srr5049781_initiator.fa> -t <srr5049781_phas.fa> -o <srr5049781_initiator_target.txt>

phasihunter deg -i <degradome.map> -q <osa_mirna.fa> -j <srr5049781_initiator_target.txt> -t <oryza_sativa_cdna.fa> -o [srr5049781_initiator_verified.txt]
```

6. phasirna target prediction and verification
```bash
phasihunter target -q <srr5049781_phasirna.fa> -t <oryza_sativa_cdna.fa> -o <srr5049781_phasirna_target.txt>

phasihunter deg -i <degradome.map> -q <srr5049781_phasirna.fa> -j <srr5049781_phasirna_target.txt> -t <oryza_sativa_cdna.fa> -o [srr5049781_phasirna_target_verified.txt]
```

# copyright
copyright © crop bioinformatics group (cbi), college of agricultural, nanjing agricultural university

free for academic use. for commercial use, please contact us (huangji@njau.edu.cn)