# -*- encoding: utf-8 -*-
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
'''
@File    :   plot_venn.py
@Time    :   2022/11/11 10:16:22
@Author  :   zrf 
@Version :   0.0
@Contact :   zrfeng1@gmail.com
@License :   (C)Copyright 2022-, NJAU-CBI
@Desc    :   None
'''
help = '''
    python3 plot_venn.py -h for help
        options:
            -a: afile
            -b: bfile
            -o: outfile, venn digram; <venn.png>
            -t1: afile lable; <a>
            -t2: bfile lable; <b>
            -A:  afile - bfile; <a_specific>
            -B:  bfile - afile; <b_specifc>
            -C:  afile overlap bfile; <ab.txt>
'''
version = '''
'''
out = './venn'
tag1 = 'a'
tag2 = 'b'
a_b = 'a_specific.txt'
b_a = 'b_specific.txt'
ab = 'ab.txt'

for i in range(1, len(sys.argv)):
    if sys.argv[i] == '-a':
        inp1 = sys.argv[i+1]
    elif sys.argv[i] == '-b':
        inp2 = sys.argv[i+1]
    elif sys.argv[i] == '-o':
        out = sys.argv[i+1]
    elif sys.argv[i] == '-t1':
        tag1 = sys.argv[i+1]
    elif sys.argv[i] == '-t2':
        tag2 = sys.argv[i+1]
    elif sys.argv[i] == '-A':
        a_b = sys.argv[i+1]
    elif sys.argv[i] == '-B':
        b_a = sys.argv[i+1]
    elif sys.argv[i] == '-C':
        ab = sys.argv[i+1]
    elif sys.argv[i] == '-v':
        print(version)
        sys.exit()
    elif sys.argv[i] == '-h':
        print(help)
        sys.exit()
print('Runing')

alist = np.loadtxt(inp1, dtype=str, delimiter=None)
blist = np.loadtxt(inp2, dtype=str, delimiter=None)

venn2(
    subsets = [set(alist),set(blist)],
    set_labels=(tag1,tag2),
    set_colors=('r','g'))

plt.savefig(out)

with open(a_b, 'w') as fo1:
    tmp1 = set(alist)-set(blist)
    for i in tmp1:
        fo1.write(i + "\n")

with open(b_a, 'w') as fo2:
    tmp2 = set(blist) - set(alist)
    for i in tmp2:
        fo2.write(i + "\n")

with open(ab, 'w') as fo3:
    tmp3 = set(alist) & set(blist)
    for i in tmp3:
        fo3.write(i + "\n")

print('Plot finished')
print(f'there are four outfile, 1: {out}.png, 2: {a_b}, 3: {b_a}, 4: {ab}')