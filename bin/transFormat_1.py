import sys
from Bio import SeqIO
import os
from customeDic import *

help = '<new_PHAS_Loci_file.txt> <deg_phasiRNA_target> <deg_initiator> <outdir> <new_initiator.txt> <new_phasiRNA.txt>'
version = ''
for i in range(1, len(sys.argv)):
    if sys.argv[i] == '-v':
        print(version)
        sys.exit()
    elif sys.argv[i] == '-h':
        print(help)
        sys.exit()

PHAS_Loci_file = sys.argv[1]
deg_phasiRNA_target = sys.argv[2]
deg_initiator = sys.argv[3]
outdir = sys.argv[4]
initiator_file = sys.argv[5]
phasiRNA_txt = sys.argv[6]

# output file
new_PHAS_Loci_file = outdir + '/database_PHAS_Loci.txt'
new_initiator_file = outdir + '/database_initiator.txt'
new_phasiRNA_file = outdir + '/database_phasiRNA.txt'
fo_new_PHAS_Loci_file = open(new_PHAS_Loci_file, 'w')
fo_new_initiator_file = open(new_initiator_file, 'w')
fo_new_phasiRNA_file = open(new_phasiRNA_file, 'w')


print('Loading deg_phasiRNA_target file')
with open(deg_phasiRNA_target, 'r') as fn:
    Target_verified_dic = nestedDic()
    phasiRNA_Target_dic = OneSetDic()
    phasiRNA_Target_number = OneSetDic()
    phasiRNA_target_verified = nestedDic()
    for line in fn:
        line1 = line.strip()
        if line1.startswith('Category'):
            continue
        l = line.strip().split("\t")
        Category,Small_RNA,Target_gene,sRNA_loc,Deg_loc,Deg_count,sRNA_seq,Shift,Gene_annotation,Library = \
        l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8],l[9]
        PHAS_Loci = '.'.join(Small_RNA.split('.')[:-1])
        Target_verified_dic[PHAS_Loci] = True
        phasiRNA_Target_dic[Small_RNA].add((Category,Target_gene,sRNA_loc,Deg_loc,Deg_count,sRNA_seq,Shift,Library))
        phasiRNA_Target_number[Small_RNA].add(Target_gene)
        phasiRNA_target_verified[Small_RNA] = True

print('Loading deg initiator file')
with open(deg_initiator, 'r') as fn:
    initiator_verified_dic = nestedDic()
    initiator_number_dic = OneSetDic()
    deg_MTI = nestedDic()
    deg_MTI_info = OneSetDic()
    for line in fn:
        line1 = line.strip()
        if line1.startswith('Category'):
            continue
        l = line.strip().split("\t")
        Category,Small_RNA,Target_gene,sRNA_loc,Deg_loc,Deg_count,sRNA_seq,Shift,Gene_annotation,Library = \
        l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8],l[9]
        PHAS_Loci = Target_gene
        initiator_verified_dic[PHAS_Loci] = True
        initiator_number_dic[PHAS_Loci].add(Small_RNA)
        MTI = Small_RNA + '\t' + Target_gene
        deg_MTI[MTI] = True
        deg_MTI_info[MTI].add((Deg_count,Shift,Library,Category))


print('Loading PHAS Loci file')
with open(PHAS_Loci_file, 'r') as fn:
    PHAS_Loci2Transcript = nestedDic()
    for line in fn:
        line1 = line.strip()
        l = line.strip().split("\t")
        PHAS_Loci,Chromosome,Genome_start,Genome_end,Locus_Gene,Transcript,PHAS_Gene_Type,Transcript_start,Transcript_end,Pvalue,PhaseScore,PhaseRatio,Annotation,Method_Ref,Tag,recorder = \
        l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8],l[9],l[10],l[11],l[12],l[13],l[14],l[15]
        if line1.startswith('PHAS_Loci'):
            fo_new_PHAS_Loci_file.write(line1 + '\tTarget_verified\tInitiator_Verified\tInitiator_Number\n')
            continue
        if Target_verified_dic[PHAS_Loci]:
            Target_Verified = 'yes'
        else:
            Target_Verified = 'no'
        if initiator_verified_dic[PHAS_Loci]:
            initiator_Verified = 'yes'
        else:
            initiator_Verified = 'no'
        PHAS_Loci2Transcript[PHAS_Loci] = Transcript
        initiator_number = len(initiator_number_dic[PHAS_Loci])
        fo_new_PHAS_Loci_file.write(line1+ f'\t{Target_Verified}\t{initiator_Verified}\t{initiator_number}\n')

print('Loading initiator file')
with open(initiator_file, 'r') as fn:
    for line in fn:
        line1 = line.strip()
        l = line.strip().split("\t")
        miRNA,GeneID,Cleavage_site,Target_pos_end,Posmax_min,Distance,Initiator_seq,Initiator_len,Recordmarker,Hit_type = \
            l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8],l[9]
        if line1.startswith('miRNA'):
            fo_new_initiator_file.write(f'PHAS_Gene\tmiRNA\tGeneID\tCleavege_site\tTarget_pos_end\tPosmax_min\tDistance\tInitiator_seq\tInitiator_len\tRecordmarker\tVerified\tDegradome\tCategory\tDeg_Count\tT-plot\n')
            continue
        gene = PHAS_Loci2Transcript[GeneID]
        MTI = miRNA + '\t' + GeneID
        if deg_MTI[MTI]:
            Verified = 'yes'
        else:
            Verified = 'no'
        if len(deg_MTI_info[MTI]) == 0:
            Degradome = '-'
            Category = '-'
            Deg_count = '-'
            T_plot = '-'
            fo_new_initiator_file.write(f'{GeneID}\t{miRNA}\t{gene}\t{Cleavage_site}\t{Target_pos_end}\t{Posmax_min}\t{Distance}\t{Initiator_seq}\t{Initiator_len}\t{Recordmarker}\t{Verified}\t{Degradome}\t{Category}\t{Deg_count}\t{T_plot}\n')
        elif len(deg_MTI_info[MTI]) >= 1:
            for i in deg_MTI_info[MTI]:
                Degradome = i[2]
                Category = i[3]
                Deg_count = i[0]
                Shift = i[1]
                T_plot = Degradome + '_' + miRNA + '-' + GeneID + '.png' 
                fo_new_initiator_file.write(f'{GeneID}\t{miRNA}\t{gene}\t{Cleavage_site}\t{Target_pos_end}\t{Posmax_min}\t{Distance}\t{Initiator_seq}\t{Initiator_len}\t{Recordmarker}\t{Verified}\t{Degradome}\t{Category}\t{Deg_count}\t{T_plot}\n')



print('Loading phasiRNA txt file')
with open(phasiRNA_txt, 'r') as fn:
    for line in fn:
        line1 = line.strip()
        l = line.strip().split("\t")
        PhasiRNA,Chromosome,Genome_Start,Genome_End,Strand,Sequence,Length,PHAS_Loci,Transcript_id,Start_in_Transcript,End_in_Transcript = \
            l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8],l[9],l[10]
        PHAS_Loci = '.'.join(PhasiRNA.split('.')[:-1])
        Verified_Target_num = len(phasiRNA_Target_number[PhasiRNA])
        if line1.startswith('PhasiRNA'):
            fo_new_phasiRNA_file.write(f'{line1}\tT_Verified\tVerified_Target_Num\tVerified_Target\tDegradome_Lib\tCategory\tSmallRNA_Loc\tDeg_Loc\tDeg_Count\tShift\tT_plot\n')
            continue
        if phasiRNA_target_verified[PhasiRNA]:
            Verified_Target = 'yes'
        else:
            Verified_Target = 'no'
        if Verified_Target_num >= 1:
            for i in phasiRNA_Target_dic[PhasiRNA]:
                Target_gene = i[1]
                Degradome_Lib = i[7]
                Category = i[0]
                Small_RNA_loc = i[2]
                Deg_loc = i[3]
                Deg_count = i[4]
                Shift = i[6]
                T_plot = Degradome_Lib + '_' + PhasiRNA + '-' + Target_gene + '.png' 
                fo_new_phasiRNA_file.write(f'{line1}\t{Verified_Target}\t{Verified_Target_num}\t{Target_gene}\t{Degradome_Lib}\t{Category}\t{Small_RNA_loc}\t{Deg_loc}\t{Deg_count}\t{Shift}\t{T_plot}\n')
        elif Verified_Target_num == 0:
            Target_gene = '-'
            Degradome_Lib = '-'
            Category = '-'
            Small_RNA_loc = '-'
            Deg_loc = '-'
            Deg_count = '-'
            Shift = '-'
            T_plot = '-'
            fo_new_phasiRNA_file.write(f'{line1}\t{Verified_Target}\t{Verified_Target_num}\t{Target_gene}\t{Degradome_Lib}\t{Category}\t{Small_RNA_loc}\t{Deg_loc}\t{Deg_count}\t{Shift}\t{T_plot}\n')



fo_new_PHAS_Loci_file.close()
fo_new_initiator_file.close()
fo_new_phasiRNA_file.close()