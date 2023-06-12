# 提取用于构建数据库的PHAS Loci(直接抛弃intron和Intergenic区域)
cat all.24.integration.s | grep -v cluster_number | grep Intergenic > Intergenic.24
cat all.24.integration.s | grep -v cluster_number | grep -v Intergenic | grep -v intron > transcript.24
head -n 1 all.24.integration.s > header.24
cat header.24 >> PHAS_Loci_database.24
cat transcript.24 >> PHAS_Loci_database.24
rm Intergenic.24 transcript.24 header.24

## 获取phasiRNA.fa
[ -e database.phasiRNA.fa.24 ] && rm database.phasiRNA.fa.24
[ -e database_phasiRNA.tab.24 ] && rm database_phasiRNA.tab.24

seqkit fx2tab all.24.phasiRNA.fa | grep -P "HC\d+" > pool.24
# 使用HC的query从HC的pool中搜索
cat all.24.integration.s | grep -v intron| phasiHunter read 15 | grep -P "\sHC\s" | awk '{print $1}' | sort | uniq | while read id 
do
fgrep -w $id pool.24
done >> database_phasiRNA.tab.24

seqkit fx2tab all.24.phasiRNA.fa | grep -P "PC\d+" > pool.24
cat all.24.integration.s | grep -v intron| phasiHunter read 15 | grep -P "\sPC\s" | awk '{print $1}' | sort | uniq | while read id 
do
fgrep -w $id pool.24
done >> database_phasiRNA.tab.24

seqkit fx2tab all.24.phasiRNA.fa | grep -P "HF\d+" > pool.24
cat all.24.integration.s | grep -v intron|  phasiHunter read 15 | grep -P "\sHF\s" | awk '{print $1}' | sort | uniq | while read id 
do
fgrep -w $id pool.24
done >> database_phasiRNA.tab.24

seqkit fx2tab all.24.phasiRNA.fa | grep -P "PF\d+" > pool.24
cat all.24.integration.s | grep -v intron | phasiHunter read 15 | grep -P "\sPF\s" | awk '{print $1}' | sort | uniq | while read id 
do
fgrep -w $id pool.24
done >> database_phasiRNA.tab.24

seqkit fx2tab all.24.phasiRNA.fa | grep -P "HG\d+" > pool.24
cat all.24.integration.s | grep -v intron|  phasiHunter read 15 | grep -P "\sHG\s" | awk '{print $1}' | sort | uniq | while read id 
do
fgrep -w $id pool.24
done >> database_phasiRNA.tab.24

seqkit fx2tab all.24.phasiRNA.fa | grep -P "PG\d+" > pool.24
cat all.24.integration.s | grep -v intron | phasiHunter read 15 | grep -P "\sPG\s" | awk '{print $1}' | sort | uniq | while read id 
do
fgrep -w $id pool.24
done >> database_phasiRNA.tab.24

sort database_phasiRNA.tab.24| uniq | seqkit tab2fx > database.phasiRNA.fa.24
rm pool.24 database_phasiRNA.tab.24

## 获取PHAS.fa
[ -e database_PHAS.fa.24 ] && rm database_PHAS.fa.24
[ -e database_PHAS.tab.24 ] && rm database_PHAS.tab.24

seqkit fx2tab all.24.PHAS.fa | grep -P "HC\d+" | awk '{print $1,$2,$3,$4,$5,$6"\t"$8}' > pool.24
cat all.24.integration.s | grep -v intron| phasiHunter read 15 | grep -P "\sHC\s" | awk '{print $1}' | sort | uniq | while read id 
do
fgrep $id pool.24
done >> database_PHAS.tab.24

seqkit fx2tab all.24.PHAS.fa | grep -P "PC\d+" | awk '{print $1,$2,$3,$4,$5,$6"\t"$8}' > pool.24
cat all.24.integration.s | grep -v intron| phasiHunter read 15 | grep -P "\sPC\s" | awk '{print $1}' | sort | uniq | while read id 
do
fgrep $id pool.24
done >> database_PHAS.tab.24

seqkit fx2tab all.24.PHAS.fa | grep -P "HF\d+" | awk '{print $1,$2,$3,$4,$5,$6"\t"$8}' > pool.24
cat all.24.integration.s | grep -v intron|  phasiHunter read 15 | grep -P "\sHF\s" | awk '{print $1}' | sort | uniq | while read id 
do
fgrep $id pool.24
done >> database_PHAS.tab.24

seqkit fx2tab all.24.PHAS.fa | grep -P "PF\d+" | awk '{print $1,$2,$3,$4,$5,$6"\t"$8}' > pool.24
cat all.24.integration.s | grep -v intron | phasiHunter read 15 | grep -P "\sPF\s" | awk '{print $1}' | sort | uniq | while read id 
do
fgrep $id pool.24
done >> database_PHAS.tab.24

seqkit fx2tab all.24.PHAS.fa | grep -P "HG\d+" | awk '{print $1,$2,$3,$4,$5,$6"\t"$8}' > pool.24
cat all.24.integration.s | grep -v intron | grep -v Intergenic |  phasiHunter read 15 | grep -P "\sHG\s" | awk '{print $2}' | sort | uniq | while read id 
do
fgrep $id pool.24
done >> database_PHAS.tab.24

seqkit fx2tab all.24.PHAS.fa | grep -P "PG\d+" | awk '{print $1,$2,$3,$4,$5,$6"\t"$8}' > pool.24
cat all.24.integration.s | grep -v intron | grep -v Intergenic | phasiHunter read 15 | grep -P "\sPG\s" | awk '{print $2}' | sort | uniq | while read id 
do
fgrep $id pool.24
done >> database_PHAS.tab.24

sort database_PHAS.tab.24| uniq | seqkit tab2fx > database_PHAS.fa.24
rm pool.24 database_PHAS.tab.24