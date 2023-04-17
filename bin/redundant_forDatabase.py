import sys
from customeDic import *
help = '''
'''
version = '''
'''
for i in range(1, len(sys.argv)):
    if sys.argv[i] == '-i':
        inp = sys.argv[i+1]
    elif sys.argv[i] == '-o':
        out = sys.argv[i+1]
    elif sys.argv[i] == '-v':
        print(version)
        sys.exit()
    elif sys.argv[i] == '-h':
        print(help)
        sys.exit()


inp = '/home/riceUsers/fzr/fzr_data/NewDownload/zma/trimmed_format_filter_cp10m/phasiHunter/all/all.integration.s'

with open(inp, 'r') as fn:
    for line in fn:
        line1 = line.strip()
        l = line.strip().split("\t")
        geneid, method, ref, hcc, hgc, hfc, pcc, pgc, pfc, hcp, hgp, hfp, pcr, pgr, pfr, pcs, pgs, pfs = \
        l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10], l[11], l[12], l[13], l[14], l[15], l[16], l[17]
