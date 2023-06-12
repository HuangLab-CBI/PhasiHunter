# usage \fd * | xargs phasiHunter addcol
import sys
import pandas as pd

fo = open('combine_tab.txt', 'w')

frames = []
for i in range(1, len(sys.argv)):
    inp = sys.argv[i]
    name = inp.replace('./', '').replace('/', '_').replace('dedup', '').replace('.abun', '').replace('__','_')
    data = pd.read_table(inp, sep="\t", header=None, names=['sRNA', 'seq', name])[['seq', name]]
    # data = pd.read_table(inp, sep="\t", header=None)[列名或者索引]
    frames.append(data)


for i in range(0,len(frames)):
    if i == 0:
        df = frames[0]
    else:
        df = pd.merge(df, frames[i], on='seq', how='outer')

df.to_csv('./combine_tab.txt', sep="\t", index=False, header=True, na_rep=0)