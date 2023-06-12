
import sys
sys.path.append('/home/riceUsers/fzr/phasiHunter/bin')
from customeDic import *

# 序列一样直接去重
dic = nestedDic()
count = 0
count_dic = OneDepDic()
for i in range(1, len(sys.argv)):
    inp = sys.argv[i]
    count += 1
    with open(inp, 'r') as fn:
        for line in fn:
            line1 = line.strip()
            l = line.strip().split("\t")
            tmp,seq,abun = l[0],l[1],l[2]
            sRNA = tmp.split('#')[0]
            dic[seq][sRNA] = float(abun)

for seq in dic:
    sum = 0
    for srna in dic[seq]:
        sum += dic[seq][srna]
    print(f'{seq}\t{sum/count}')


# for i in count_dic:
#     abun = sum(count_dic[i])/count
#     print(f'{i}\t{abun}')
