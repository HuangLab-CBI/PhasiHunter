import sys
from Bio import SeqIO
import gzip
'''
@File         :   filter.py
@Time         :   2022/11/30 19:59:14
@Author       :   Ji Huang
@Contributors :   Ji Hunag, Baoyi Zhang, Zerong Feng
@Version      :   1.0
@Contact      :   zrfeng1@gmail.com
@License      :   (C)Copyright 2022-, NJAU-CBI
@Desc         :   None
'''
# ----------> variable declaration <------------
inputfile = ""
outputfilename = ""
min = 16
max = 26
min_read_count = 1.0
compression=0
help = '''
Usage:
    -i     [str]   inputfile name, forced
    -o     [str]   filename of output file, forced 
    -min   [int]   Minimum length of read sequence, default:16
    -max   [int]   the maximum length of the read sequence, default:26
    -count [float] Minimum expression threshold, default:1.0
    -c     the output file is compressed using gzip 
    -h     this information//
'''
version = '''
'''

def try_gzip_open(file, type): # Zerong Feng append
    """_summary_

    Arguments:
        file {file} -- _description_
        type {mode} -- file open mode, e.g. 'rt', 'wt'

    Returns:
        object -- File handle object
    """
    if ".gz" in file:
        fn = gzip.open(file, type)
        return fn
    else:
        fn = open(file, type)
        return fn
# ----------> end <------------

for i in range(1, len(sys.argv)):
    if sys.argv[i] == '-i':
        inputfile = sys.argv[i+1]
    elif sys.argv[i] == '-o':
        outputfilename = sys.argv[i+1]
    elif sys.argv[i]=='-min':
        min=int(sys.argv[i+1])
    elif sys.argv[i]=='-max':
        max=int(sys.argv[i+1])
    elif sys.argv[i]=='-count':
        min_read_count=float(sys.argv[i+1])
    elif sys.argv[i]=='-c':
        compression=1
    elif sys.argv[i] == '-v':
        print(version)
        sys.exit()
    elif sys.argv[i] == '-h':
        print(help)
        sys.exit()

# ----------> for enhancing robustness <------------
if inputfile=="" or outputfilename=="":
    print("python filter.py -h")


if compression == 0:
    o=open(outputfilename, 'w')
elif compression == 1:
    if ".gz" in outputfilename:
        o = gzip.open(outputfilename, 'wt')
    else:
        o = gzip.open(outputfilename+".gz", 'wt')

for query in SeqIO.parse(try_gzip_open(inputfile, 'rt'), 'fasta'):
    name = query.id
    seq = str(query.seq)
    try:
        exp = float(name.split("@")[1])
    except IndexError:
        exp = float(name.split("_")[1])
    seqlen = len(seq)
    if exp>=min_read_count and seqlen>=min and seqlen<=max:
        o.write(">"+name+"\n"+seq+"\n")
o.close()