# usage \fd * | xargs phasiHunter addcol
import sys
import pandas as pd

fo = open('combine_tab.txt', 'w')

frames = []
for i in range(1, len(sys.argv)):
    inp = sys.argv[i]
    name = inp.replace('./', '').replace('/', '_').replace('_dedup', '').replace('.abun', '')
    data = pd.read_table(inp, sep="\t", header=None, names=['sRNA', 'seq', name])[name]
    # data = pd.read_table(inp, sep="\t", header=None)[列名或者索引]
    frames.append(data)

frame1 = pd.concat(frames, axis=1)
# print(frame1)
frame1.to_csv('./combine_tab.txt', sep="\t", index=False, header=True)