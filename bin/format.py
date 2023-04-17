#!/usr/bin/env python
# -*- encoding: utf-8 -*-
import sys
from Bio import SeqIO
import time
import gzip
'''
@File         :   format.py
@Time         :   2022/11/30 19:59:14
@Author       :   Ji Huang
@Contributors :   Ji Hunag, Baoyi Zhang, Zerong Feng
@Version      :   1.0
@Contact      :   zrfeng1@gmail.com
@License      :   (C)Copyright 2022-, NJAU-CBI
@Desc         :   None
'''
help = '''
//by Ji Huang. Current version 1.0
python3 format.py [-i inputfilename] [-o outputfilename] [-it inputfile type] [-ot outputfile type] [-of fasta_title name]

    -i:    [str] inputfilename
    -o:    [str] outputfilename
    -min:  [int] min length of reads, optional, default=1
    -max:  [int] max length of reads, optional, default=100
    -n:    [int] base of normalization, default=1000000
    -p:    poly(A) or poly(T) reads will be removed if -p is specified
    -c:    [int] min count of each read, default=0
    -b:    [int] startline for read processing, only for -ot g or p. default=1
    -left: [int] the bases removed, from left. default=0
    -it:   [str] inputfile-type:
            g, GEO type, non-redundant, like GAGAGATGATGTGATGAT 12
            p, pure sequence, redundant, like GAGAGAGGAGGGGG
            f, fasta format, redundant, like >root
            f1, fasta format, non-redundant, like >root1#@1
    -ot:   [str] outputfile-type:
            g, GEO type
            p, pure sequence
            f, fasta format, redundant 
            f1, with count infor,  non-redundant, like >root#1@12
            fn, non-redundant and normalized, like >root#1@0.38
            m, like >root0001_x1
    -of:   [str] the Fasta title name in outputfile, optional, default=reads
    -C:    the output file is compressed using gzip 
    -v:    print version of format.py
    -h:    this information
//
'''
version = '''
    Version:1.0 2022/11/30 , by Ji Huang"
'''

# ----------> variable and functino declaration <------------
start_time = time.time()
inputfilename=''
outputfilename=''
min_length=1
max_length=100
formatin=''
formatout=''
fastatitle=''
nor_base=1000000
polyA=0
min_read_count=0
quzhen=0
startline=1
left=0
compression=0

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
for i in range(1,len(sys.argv)):
    if sys.argv[i]=='-i':
        inputfilename=sys.argv[i+1]
    elif sys.argv[i]=='-o':
        outputfilename=sys.argv[i+1]
    elif sys.argv[i]=='-n':
        nor_base=int(sys.argv[i+1])
    elif sys.argv[i]=='-it':
        formatin=sys.argv[i+1]
    elif sys.argv[i]=='-ot':
        formatout=sys.argv[i+1]
    elif sys.argv[i]=='-of':
        fastatitle=sys.argv[i+1]
    elif sys.argv[i]=='-min':
        min_length=int(sys.argv[i+1])
    elif sys.argv[i]=='-c':
        min_read_count=int(sys.argv[i+1])
    elif sys.argv[i]=='-p':
        polyA=1
    elif sys.argv[i]=='-left':
        left=int(sys.argv[i+1])
    elif sys.argv[i]=='-b':
        startline=int(sys.argv[i+1])
    elif sys.argv[i]=='-C':
        compression=1
    elif sys.argv[i]=='-z':
        quzhen=1
    elif sys.argv[i]=='-max':
        max_length=int(sys.argv[i+1])
    elif sys.argv[i] == '-v':
        print(version)
        sys.exit()
    elif sys.argv[i] == '-h':
        print(help)
        sys.exit()

# for enhancing robustness
if len(sys.argv)==1:
    print("type python3 format.py -h for details")
    sys.exit()

if formatin==formatout:
    print("//")
    print("WARNING: The inputfile format is identical to outputfile format!")

if fastatitle=='':
    fastatitle='reads'

if startline>1:
    print("First line be processed:", startline)
    if formatin=='f' or formatin=='f1':
        print('Warning: -b not for -it f or f1; -b was ignored.')

if inputfilename=='' or outputfilename=='':
    print("Filename input error!")
    sys.exit()

if quzhen==1 and formatout!='fn':
    print('-z must be used together with -ot fn')
    sys.exit()

# ----------> analysis <------------
string_f=[]
string_f1=[]
string_fn=[]
nf=0
nf1=0
nfn=0
if compression == 0:
    w=open(outputfilename, 'w')
elif compression == 1:
    if ".gz" in outputfilename:
        w = gzip.open(outputfilename, 'wt')
    else:
        w = gzip.open(outputfilename+".gz", 'wt')

total_seq=0

# ----------> count total number <------------
print("Counting total reads...")
ll=0 
filee=try_gzip_open(inputfilename, 'rt')
heads = ["A","G","C","T"]

if formatin=='g':
    for line in filee:
        ll+=1
        if line[0].upper() in heads:
            line1=line.strip().split('\t')
            if len(line1)>=2: # zby append
                seq=line1[0]
                seq_count=int(float(line1[1]))
                seq_count_num=str(line1[1]).replace("\n","")
                total_seq+=seq_count
elif formatin=='f':
    for ss in filee:
        if '>' in ss:
            total_seq+=1
elif formatin=='p':
    for sp in filee:
        ll+=1
        if ll>=startline:
            if sp!='' and sp!='\n':
                total_seq+=1
elif formatin=='f1':
    for s1 in SeqIO.parse(try_gzip_open(inputfilename,'rt'),'fasta'):
        s11=float(s1.id.split('@')[1])
        total_seq+=s11
print("Total sequences:", total_seq)

print("Formating......")

def poly(seq):
    if polyA==0: return 0
    polya=float(seq.count ("A"))
    polyt=float(seq.count ("T"))
    if polya/len(seq)>=0.82 or polyt/len(seq)>=0.82:
        return 1
    else:
        return 0

ll=0
if formatin=='g':
    for line in try_gzip_open(inputfilename, 'rt'):
        ll+=1
        if line[0].upper() in heads:
            line1=line.strip().split('\t')
            if len(line1)>=2: # zby append
                seq=line1[0].upper()[left:]
                seq_count=int(float(line1[1]))
                seq_count_num=str(line1[1]).replace("\n","")
                if len(str(seq))>=min_length and len(str(seq))<=max_length and poly(seq)==0 and seq_count>=min_read_count:
                    if formatout=='f':
                        for j in range(seq_count):
                            string_f.append (">"+fastatitle+'#'+str(nf)+'\n')
                            string_f.append (seq+'\n')
                            nf+=1
                    elif formatout=='p':
                        for j in range(seq_count):
                            string_f.append (seq+'\n')
                            nf+=1
                    elif formatout=='f1':
                        string_f.append ('>'+fastatitle+'#'+str(nf)+'@'+str(seq_count)+'\n')
                        string_f.append (seq+'\n')
                        nf+=1
                    elif formatout=='fn':
                        if nor_base==0:
                            string_f.append ('>'+fastatitle+'#'+str(nf)+'@'+str(seq_count)+'\n')
                            string_f.append (seq+'\n')
                        else:
                            string_f.append ('>'+fastatitle+'#'+str(nf)+'@'+str(round(float(seq_count)*nor_base/total_seq,2))+'\n')
                            string_f.append (seq+'\n')
                        nf+=1
                    elif formatout=='m':
                        string_f.append ('>'+fastatitle+str(nf)+'_x'+str(seq_count)+'\n')
                        string_f.append (seq+'\n')
                        nf+=1
                    else:
                        print("Error in outputfile type")
                        sys.exit()
elif formatin=='p':
    if formatout=='p':
        print("Did nothing!")
        sys.exit()

    p_dic={}
    for line in try_gzip_open(inputfilename, 'rt'):
        ll+=1
        if ll>=startline:
            if str(line)!='' and line!='\n':
                seq_len=len(str(line[:-1]))
                if seq_len>=min_length and seq_len<=max_length and poly(str(line[:-1]))==0:
                    if line.upper() in p_dic:
                        p_dic[line.upper()]+=1
                    else:
                        p_dic[line.upper()]=1
        for key in p_dic.keys():
            if int(p_dic[key])>=min_read_count:
                if formatout=='f':
                    for j in range(int(p_dic[key])):
                        string_f.append(">"+fastatitle+'#'+str(nf)+'\n')
                        string_f.append (key[left:])
                        nf+=1
                elif formatout=='f1':
                    string_f.append(">"+fastatitle+'#'+str(nf)+'@'+str(p_dic[key])+'\n')
                    string_f.append (key[left:])
                    nf+=1
                elif formatout=='fn':
                    if nor_base==0:
                        string_f.append ('>'+fastatitle+'#'+str(nf)+'@'+str(p_dic[key])+'\n')
                        string_f.append (key[left:])
                    else:
                        string_f.append ('>'+fastatitle+'#'+str(nf)+'@'+str(round(float(p_dic[key])*nor_base/total_seq,2))+'\n')
                        string_f.append (key[left:])
                    nf+=1
                elif formatout=='m':
                    string_f.append(">"+fastatitle+str(nf)+'_x'+str(p_dic[key])+'\n')
                    string_f.append (key[left:])
                    nf+=1
                elif formatout=='g':
                    string_f.append (str(key[left:])+'\t'+str(p_dic[key])+'\n')
                    nf+=1
                else:
                    print("Error in outputfile type")
                    sys.exit()
elif formatin=='f':
    f_dic={}
    for se in SeqIO.parse (try_gzip_open(inputfilename, 'rt'),'fasta'):
        seq_len=len(str(se.seq))
        if seq_len>=min_length and seq_len<=max_length and poly(str(se.seq))==0:
            if str(se.seq).upper() in f_dic:
                f_dic[str(se.seq).upper()]+=1
            else:
                f_dic[str(se.seq).upper()]=1
    for key in f_dic.keys():
        if f_dic[key]>=min_read_count:
            if formatout=='f':
                for j in range(int(f_dic[key])):
                    string_f.append(">"+fastatitle+'#'+str(nf)+'\n')
                    string_f.append (key[left:]+'\n')
                    nf+=1
            elif formatout=='p':
                for j in range(int(f_dic[key])):
                    string_f.append (key[left:]+'\n')
                    nf+=1
            elif formatout=='f1':
                string_f.append(">"+fastatitle+'#'+str(nf)+'@'+str(f_dic[key])+'\n')
                string_f.append (key[left:]+'\n')
                nf+=1
            elif formatout=='m':
                string_f.append(">"+fastatitle+str(nf)+'_x'+str(f_dic[key])+'\n')
                string_f.append (key[left:]+'\n')
                nf+=1
            elif formatout=='fn':
                if nor_base==0:
                    string_f.append ('>'+fastatitle+'#'+str(nf)+'@'+str(f_dic[key])+'\n')
                    string_f.append (key[left:]+'\n')
                else:
                    string_f.append ('>'+fastatitle+'#'+str(nf)+'@'+str(round(float(f_dic[key])*nor_base/total_seq,2))+'\n')
                    # normalization 
                    # seq_abundance * nor_base / total_sRNA_number
                    string_f.append (key[left:]+'\n')
                nf+=1
            elif formatout=='g':
                string_f.append (str(key[left:])+'\t'+str(f_dic[key])+'\n')
                nf+=1
            else:
                print("Error in outputfile type")
                sys.exit()
elif formatin=='f1':
    for se in SeqIO.parse (try_gzip_open(inputfilename, 'rt'),'fasta'):
        seq_len=len(str(se.seq))
        seq_count=float(se.id.split('@')[1])
        if seq_len>=min_length and seq_len<=max_length and poly(str(se.seq))==0 and seq_count>=min_read_count:
            if formatout=='f':
                for j in range(int(seq_count)):
                    string_f.append(">"+'read'+'#'+str(nf)+'\n')
                    string_f.append (str(se.seq)[left:]+'\n')
                    nf+=1
            elif formatout=='p':
                for j in range(int(seq_count)):
                    string_f.append (str(se.seq)[left:]+'\n')
                    nf+=1
            elif formatout=='f1':
                string_f.append(">"+se.id+'\n')
                string_f.append (str(se.seq)[left:]+'\n')
                nf+=1
            elif formatout=='fn':
                if fastatitle=='':
                    if nor_base==0:
                        string_f.append ('>'+se.id.split('@')[0]+'@'+str(seq_count)+'\n')
                        string_f.append (str(se.seq)[left:]+'\n')
                    else:
                        string_f.append ('>'+se.id.split('@')[0]+'@'+str(round(float(seq_count)*nor_base/total_seq,2))+'\n')
                        string_f.append (str(se.seq)[left:]+'\n')
                    nf+=1
                else:
                    if nor_base==0:
                        string_f.append ('>'+fastatitle+'@'+str(seq_count)+'\n')
                        string_f.append (str(se.seq)[left:]+'\n')
                    else:
                        string_f.append ('>'+fastatitle+'@'+str(round(float(seq_count)*nor_base/total_seq,2))+'\n')
                        string_f.append (str(se.seq)[left:]+'\n')
                    nf+=1
            elif formatout=='m':
                if fastatitle=='reads':
                    string_f.append ('>'+se.id.split('@')[0]+'_x'+str(seq_count)+'\n')
                    string_f.append (str(se.seq)[left:]+'\n')
                    nf+=1
                else:
                    string_f.append ('>'+fastatitle+'_x'+str(seq_count)+'\n')
                    string_f.append (str(se.seq)[left:]+'\n')
                    nf+=1
            elif formatout=='g':
                string_f.append (str(se.seq)[left:]+'\t'+str(seq_count)+'\n')
                nf+=1
            else:
                print("Error in outputfile type")
                sys.exit()
else:
    print("Error in inputfile type")
    sys.exit()

print("Writing......")
w.writelines(string_f)
stop_time = time.time()
print("Time elasped (second):", stop_time - start_time)
