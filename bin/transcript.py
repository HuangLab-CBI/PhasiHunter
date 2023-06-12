# 统计转录本对应的Loci的信息
import sys
from customeDic import *

help = '<gene_Loci_info> <all_de_transcript>'
version = ''
for i in range(1, len(sys.argv)):
    if sys.argv[i] == '-v':
        print(version)
        sys.exit()
    elif sys.argv[i] == '-h':
        print(help)
        sys.exit()

inp = sys.argv[1]
# 将基因座信息存为字典
dic = OneDepDic()
with open(inp, 'r') as fn:
    for line in fn:
        line1 = line.strip()
        l = line.strip().split("\t")
        length = len(l)
        if length == 3:
            loc = l[2]
            for j in l[:2]:
                if j != '':
                    if j.startswith('XP_') or j.startswith('NP_') or j.startswith('YP_'):
                        continue
                    else:
                        dic[loc].append(j)
            dic[loc] = list(set(dic[loc]))
        elif length ==1:
            continue

# 过滤得到多转录本基因座位信息
num_dic = nestedDic()
multiple_dic = OneDepDic()
for k in dic:
    number = len(dic[k])
    if number > 1:
        multiple_dic[k] = dic[k]
        num_dic[k] = number

# 分析PHAS 转录本信息, 如果一个转录本在基因组中存在，那么这个转录本将从此基因座的列表中弹出
inp2 = sys.argv[2]
with open(inp2, 'r') as fn:
    for line in fn:
        line1 = line.strip()
        l = line.strip().split("\t")
        for k in list(multiple_dic.keys()):
            if l[0] in multiple_dic[k]:
                multiple_dic[k].remove(l[0])

# 输出结果信息
foo = open('./多转录本phasiRNA产生能力分析.txt', 'w')
foo.write(f'Locus\ttotalTranscriptNumber\tcannotGeneratingphasiRNANumber\tdetail\n')
print(f'Locus\ttotalTranscriptNumber\tcannotGeneratingphasiRNANumber\tgenerationgPHASLociNumber\tgeneratingPHASability\tdetail\n')
for i in multiple_dic:
    total_num = num_dic[i]
    cannotGeneratingPHASNumber = len(multiple_dic[i])
    transcript_cangeneratingPHASNumber = total_num - cannotGeneratingPHASNumber
    gengratingAbility = transcript_cangeneratingPHASNumber/total_num
    foo.write(f'{i}\t{total_num}\t{transcript_cangeneratingPHASNumber}\t{transcript_cangeneratingPHASNumber}\t{gengratingAbility}\t{multiple_dic[i]}\n')
    print(f'{i}\t{total_num}\t{transcript_cangeneratingPHASNumber}\t{transcript_cangeneratingPHASNumber}\t{gengratingAbility}\t{multiple_dic[i]}')
foo.close()