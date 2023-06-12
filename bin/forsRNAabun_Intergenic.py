import sys
integration = sys.argv[1]
phasiRNA = sys.argv[2]
from customeDic import *
from Bio import SeqIO
with open(integration, 'r') as fn:
    fucku_dic = nestedDic()
    for line in fn:
        line1 = line.strip()
        l = line.strip().split("\t")
        if line1.startswith('#') or line1.startswith('>'):
            continue
        target, strand, loc, cnm, srna, seq, length, a, b, c, d, pvalue, anno, record = \
        l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8],l[9],l[10],l[11],l[12],l[13]
        abun = srna.split('@')[1]
        key = f'{target}\t{seq}\t{record}'
        fucku_dic[key] = (srna, abun)

for query in SeqIO.parse(phasiRNA, 'fasta'):
    desc = query.description
    seq = str(query.seq)
    target = query.id    
    record = desc.split(' ')[2]
    anno = desc.split(' ')[1]
    if anno == 'Other:Intergenic':
        query_key = f'{target}\t{seq}\t{record}'
        print(f'{fucku_dic[query_key][0]}\t{seq}\t{fucku_dic[query_key][1]}')
