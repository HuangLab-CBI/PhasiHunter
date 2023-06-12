import sys
from customeDic import *
import os
help = '''
    phasiHunter database1 <initiator.txt> <phasiRNA.txt> <PHAS_tmp.txt>
'''
version = '''
'''
for i in range(1, len(sys.argv)):
    if sys.argv[i] == '-v':
        print(version)
        sys.exit()
    elif sys.argv[i] == '-h':
        print(help)
        sys.exit()

ini = sys.argv[1]  
phasiRNA = sys.argv[2]
PHAS_tmp = sys.argv[3]

TWD = os.getcwd()
PHAS = open(f'{TWD}/PHAS.txt', 'w')

PHAS.write(f'PHAS_Gene\tChromosome\tGenome_start\tGenome_End\tLocus_Gene\tTranscript_id\tPHAS_Gene_Type\tStart_in_transcript\tEnd_in_transcript\tP-value\tPhase_score\tAnnotation\tMethod\tReference_sequence\tRecorder\tPlot\tTarget_verified\tInitiator_num\tInitiator_verified\n')
with open(phasiRNA, 'r') as fn:
    phasiRNA_dic = nestedDic()
    for line in fn:
        line1 = line.strip()
        if line1.startswith('PhasiRNA'):
            continue
        l = line.strip().split("\t")
        PhasiRNA,Chromosome,Genome_Start,Genome_End,Strand,Sequence,Length,PHAS_Gene,Transcript_id,Start_in_Transcript,end_in_Transcript,T_verified,Verified_Target_num,Verified_Target,Degradome_Lib,Category,SmallRNA_loc,Deg_loc,Deg_Count,Shift,T_plot = \
            l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8],l[9],l[10],l[11],l[12],l[13],l[14],l[15],l[16],l[17],l[18],l[19],l[20]
        phasiRNA_dic[PHAS_Gene] = True

with open(ini, 'r') as fn:
    ini_dic_ver = nestedDic()
    ini_dic_num = OneDepDic()
    for line in fn:
        line1 = line.strip()
        if line1.startswith('PHAS_Gene'):
            continue
        l = line.strip().split("\t")
        PHAS_Gene,miRNA,GeneID,Cleavage_site,Target_pos_end,Posmax,Distance,Initiator_seq,Initiator_len,Recordmarker,Verified,Degradome,Category,Deg_Count,T_plot = \
            l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8],l[9],l[10],l[11],l[12],l[13],l[14]
        ini_dic_ver[PHAS_Gene] = True
        ini_dic_num[PHAS_Gene].append(miRNA)

with open(PHAS_tmp, 'r') as fn:
    for line in fn:
        line1 = line.strip()
        if line1.startswith('PHAS_Gene'):
            continue
        l = line.strip().split("\t")
        PHAS_Gene,Chromosome,Genome_start,Genome_End,Locus_Gene,Transcript_id,PHAS_Gene_Type,Start_in_transcript,End_in_transcript,P_value,Phase_score,Annotation,Method,Reference_sequence,Recorder,Plot = \
            l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8],l[9],l[10],l[11],l[12],l[13],l[14], l[15]
        if phasiRNA_dic[PHAS_Gene]:
            Target_verified = 'yes'
        else:
            Target_verified = 'no'
        
        if ini_dic_ver[PHAS_Gene]:
            Initiator_verified = 'yes'
        else:
            Initiator_verified = 'no'

        Initiator_num = len(set(ini_dic_num[PHAS_Gene]))

        out = f'{line1}\t{Target_verified}\t{Initiator_num}\t{Initiator_verified}'
        PHAS.write(out + "\n")
