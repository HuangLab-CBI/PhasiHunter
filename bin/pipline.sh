# 提取用于构建数据库的PHAS Loci(直接抛弃intron和Intergenic区域)
cat all.integration.s | grep -v cluster_number | grep Intergenic > Intergenic
cat all.integration.s | grep -v cluster_number | grep -v Intergenic | grep -v intron > transcript
head -n 1 all.integration.s > header
cat header >> PHAS_Loci_database
cat transcript >> PHAS_Loci_database
rm Intergenic transcript header

## 获取phasiRNA.fa
[ -e database.phasiRNA.fa ] && rm database.phasiRNA.fa
[ -e database_phasiRNA.tab ] && rm database_phasiRNA.tab

seqkit fx2tab all.phasiRNA.fa | grep -P "HC\d+" > pool
# 使用HC的query从HC的pool中搜索
cat all.integration.s | grep -v intron| phasiHunter read 15 | grep -P "\sHC\s" | awk '{print $1}' | sort | uniq | while read id 
do
fgrep -w $id pool
done >> database_phasiRNA.tab

seqkit fx2tab all.phasiRNA.fa | grep -P "PC\d+" > pool
cat all.integration.s | grep -v intron| phasiHunter read 15 | grep -P "\sPC\s" | awk '{print $1}' | sort | uniq | while read id 
do
fgrep -w $id pool
done >> database_phasiRNA.tab

seqkit fx2tab all.phasiRNA.fa | grep -P "HF\d+" > pool
cat all.integration.s | grep -v intron|  phasiHunter read 15 | grep -P "\sHF\s" | awk '{print $1}' | sort | uniq | while read id 
do
fgrep -w $id pool
done >> database_phasiRNA.tab

seqkit fx2tab all.phasiRNA.fa | grep -P "PF\d+" > pool
cat all.integration.s | grep -v intron | phasiHunter read 15 | grep -P "\sPF\s" | awk '{print $1}' | sort | uniq | while read id 
do
fgrep -w $id pool
done >> database_phasiRNA.tab

seqkit fx2tab all.phasiRNA.fa | grep -P "HG\d+" > pool
cat all.integration.s | grep -v intron|  phasiHunter read 15 | grep -P "\sHG\s" | awk '{print $1}' | sort | uniq | while read id 
do
fgrep -w $id pool
done >> database_phasiRNA.tab

seqkit fx2tab all.phasiRNA.fa | grep -P "PG\d+" > pool
cat all.integration.s | grep -v intron | phasiHunter read 15 | grep -P "\sPG\s" | awk '{print $1}' | sort | uniq | while read id 
do
fgrep -w $id pool
done >> database_phasiRNA.tab

sort database_phasiRNA.tab | uniq | seqkit tab2fx > database.phasiRNA.fa
rm pool database_phasiRNA.tab

## 获取PHAS.fa
[ -e database_PHAS.fa ] && rm database_PHAS.fa
[ -e database_PHAS.tab ] && rm database_PHAS.tab

seqkit fx2tab all.PHAS.fa | grep -P "HC\d+" | awk '{print $1,$2,$3,$4,$5,$6"\t"$8}' > pool
cat all.integration.s | grep -v intron| phasiHunter read 15 | grep -P "\sHC\s" | awk '{print $1}' | sort | uniq | while read id 
do
fgrep $id pool
done >> database_PHAS.tab

seqkit fx2tab all.PHAS.fa | grep -P "PC\d+" | awk '{print $1,$2,$3,$4,$5,$6"\t"$8}' > pool
cat all.integration.s | grep -v intron| phasiHunter read 15 | grep -P "\sPC\s" | awk '{print $1}' | sort | uniq | while read id 
do
fgrep $id pool
done >> database_PHAS.tab

seqkit fx2tab all.PHAS.fa | grep -P "HF\d+" | awk '{print $1,$2,$3,$4,$5,$6"\t"$8}' > pool
cat all.integration.s | grep -v intron|  phasiHunter read 15 | grep -P "\sHF\s" | awk '{print $1}' | sort | uniq | while read id 
do
fgrep $id pool
done >> database_PHAS.tab

seqkit fx2tab all.PHAS.fa | grep -P "PF\d+" | awk '{print $1,$2,$3,$4,$5,$6"\t"$8}' > pool
cat all.integration.s | grep -v intron | phasiHunter read 15 | grep -P "\sPF\s" | awk '{print $1}' | sort | uniq | while read id 
do
fgrep $id pool
done >> database_PHAS.tab

seqkit fx2tab all.PHAS.fa | grep -P "HG\d+" | awk '{print $1,$2,$3,$4,$5,$6"\t"$8}' > pool
cat all.integration.s | grep -v intron | grep -v Intergenic |  phasiHunter read 15 | grep -P "\sHG\s" | awk '{print $2}' | sort | uniq | while read id 
do
fgrep $id pool
done >> database_PHAS.tab

seqkit fx2tab all.PHAS.fa | grep -P "PG\d+" | awk '{print $1,$2,$3,$4,$5,$6"\t"$8}' > pool
cat all.integration.s | grep -v intron | grep -v Intergenic | phasiHunter read 15 | grep -P "\sPG\s" | awk '{print $2}' | sort | uniq | while read id 
do
fgrep $id pool
done >> database_PHAS.tab

sort database_PHAS.tab | uniq | seqkit tab2fx > database_PHAS.fa
rm pool database_PHAS.tab