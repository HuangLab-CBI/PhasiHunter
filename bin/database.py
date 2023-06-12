import sys
from Bio import SeqIO
import os
from customeDic import *
spe = 'Zm'
phase_length = 21
help = '''
    phasiHunter database <PHAS_Loci_database> <feature_table> <database_phasiRNA.fa> <all.integtation.o> <database_PHAS.fa> <refcdna.fa> <initiator_prediction> <initiator_verified> <zma.phasiRNA.deg> <1|2> <0|1> <spe>
    <1|2> means print Cat_1 only or Cat_1 and Cat_2
    <0|1> means print Shift 0 or Shif 1
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

PHAS_Loci = sys.argv[1]
feature_table = sys.argv[2]
phasiRNA = sys.argv[3]
integration_o = sys.argv[4]
PHAS = sys.argv[5]
ref_cdna = sys.argv[6]
initiator_prediction = sys.argv[7]
initiator_verified = sys.argv[8]
phasiRNA_target = sys.argv[9]
Cat_cutoff = int(sys.argv[10])
Shift_cutoff = int(sys.argv[11])
spe = sys.argv[12]
phase_length = int(sys.argv[13])

# PHAS_Gene = sys.argv[4]
# initiator =sys.argv[5]
# target = sys.argv[6]

TWD = os.getcwd()

PHAS_txt = open(f'{TWD}/PHAS_tmp.txt', 'w')
phasiRNA_txt = open(f'{TWD}/phasiRNA.txt', 'w')
phasiRNA_fa = open(f'{TWD}/phasiRNA.fa', 'w')
PHAS_fa = open(f'{TWD}/PHAS.fa', 'w')
PHAS_fa_for_initiator = open(f'{TWD}/forInitiatorPHAS.fa', 'w')
ini_out = open(f'{TWD}/initiator.txt', 'w')

def phasiRNAFormat():
    print(f'Loading phasiRNA_target...')
    phasiRNA_deg_dic1 = OneDepDic()
    phasiRNA_deg_dic2 = OneDepDic()
    with open(phasiRNA_target, 'r') as fn:
        for line in fn:
            line1 = line.strip()
            l = line.strip().split("\t")
            if line1.startswith('phase-target'):
                library = line1.split('/')[1].split('_')[0].split('.')[0]
                continue
            if line1.startswith('Category'):
                continue
            cat, sRNA, target, sRNA_loc, Deg_loc, Deg_count, seq, Shift, Gene_annotation = \
                l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8]
            if Cat_cutoff == 1 and Shift_cutoff == 0:
                if cat == 'Cat_1':
                    if Shift == '0':
                        T_plot = f'{sRNA}-{target}.png'
                        phasiRNA_deg_dic1[sRNA].append(target)
                        phasiRNA_deg_dic2[sRNA + '\t' + target].append([library, cat, sRNA_loc, Deg_loc, Deg_count, Shift, T_plot])
            elif Cat_cutoff == 2 and Shift_cutoff == 0:
                if cat == 'Cat_1' or cat == 'Cat_2':
                    if Shift == '0':
                        T_plot = f'{sRNA}-{target}.png'
                        phasiRNA_deg_dic1[sRNA].append(target)
                        phasiRNA_deg_dic2[sRNA + '\t' + target].append([library, cat, sRNA_loc, Deg_loc, Deg_count, Shift, T_plot])
            elif Cat_cutoff == 1 and Shift_cutoff == 1:
                if cat == "Cat_1":
                    if Shift == '1' or Shift == '-1':
                        T_plot = f'{sRNA}-{target}.png'
                        phasiRNA_deg_dic1[sRNA].append(target)
                        phasiRNA_deg_dic2[sRNA + '\t' + target].append([library, cat, sRNA_loc, Deg_loc, Deg_count, Shift, T_plot])
            elif Cat_cutoff == 2 and Shift_cutoff == 1:
                if cat == 'Cat_1' or cat == 'Cat_2':
                    if Shift == '1' or Shift == '-1':
                        T_plot = f'{sRNA}-{target}.png'
                        phasiRNA_deg_dic1[sRNA].append(target)
                        phasiRNA_deg_dic2[sRNA + '\t' + target].append([library, cat, sRNA_loc, Deg_loc, Deg_count, Shift, T_plot])

    print(f'Loading initiator_verified')
    initiator_verified_dic = OneDepDic()
    with open(initiator_verified, 'r') as fn:
        for line in fn:
            line1 = line.strip()
            if line1.startswith('phase-ini'):
                library = line1.split('/')[1].split('_')[0].split('.')[0]
                continue
            l = line.strip().split("\t")
            cat, initiator, target = l[0],l[1],l[2]
            deg_count = l[5]
            T_plot = f'{initiator}-{target}.png'
            initiator_verified_dic[initiator + '@' + target].append([cat, library, deg_count, T_plot])

    print('Loading phasiRNA.fa')
    phasiRNA_dic = nestedDic()
    for query in SeqIO.parse(phasiRNA, 'fasta'):
        descr = query.description
        seq = str(query.seq)
        phasiRNA_dic[descr] = seq

    print('Loading PHAS.fa')
    PHAS_dic = nestedDic()
    for query in SeqIO.parse(PHAS, 'fasta'):
        descr = query.description
        seq = str(query.seq)
        PHAS_dic[descr] = seq
    
    print('Loading ref transcriptome')
    ref_dic = nestedDic()
    for query in SeqIO.parse(ref_cdna, 'fasta'):
        id = query.name
        seq = str(seq)
        ref_dic[id]=seq

    print('Loading feature_table table...')
    locus2transcript_dic = OneDepDic()
    transcript2_sympol = OneDepDic()
    transcript2chr_dic = nestedDic()
    transcript2anno_dic = nestedDic()
    transcript2_GeneID = nestedDic()
    PHAS_txt.write(f"PHAS_Gene\tChromosome\tGenome_start\tGenome_End\tLocus_Gene\tTranscript_id\tPHAS_Gene_Type\tStart_in_transcript\tEnd_in_transcript\tP-value\tPhase_score\tAnnotation\tMethod\tReference_sequence\tRecorder\tPlot" + "\n")
    phasiRNA_txt.write(f'PhasiRNA\tChromosome\tGenome_Start\tGenome_End\tStrand\tSequence\tLength\tPHAS_Gene\tTranscript_id\tStart_in_Transcript\tend_in_Transcript\tT_verified\tVerified_Target_num\tVerified_Target\tDegradome_Lib\tCategory\tSmallRNA_loc\tDeg_loc\tDeg_Count\tShift\tT_plot' + '\n')
    ini_out.write(f'PHAS_Gene\tmiRNA\tGeneID\tCleavage_site\tTarget_pos_end\tPosmax/min\tDistance\tInitiator_seq\tInitiator_len\tRecordmarker\tVerified\tDegradome\tCategory\tDeg_Count\tT-plot'+ "\n")

    print(f'Loading feature table...')
    with open(feature_table, 'r') as fn:
        for line in fn:
            line1 = line.strip()
            if line.startswith('#'):
                continue
            l = line.strip().split("\t")
            feature, class_, assembly, assembly_unit, seq_type, chromosome, genomic_accession, start,	end,	strand,	product_accession,\
                non_redundant_refseq, related_accession,	name,	symbol,	GeneID,	locus_tag,	feature_interval_length \
                        = l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10], l[11], l[12], l[13], l[14], l[15], l[16], l[17]
            if symbol != '':
                if product_accession != '':
                    locus2transcript_dic[symbol].append(product_accession)
                if related_accession != '':
                    locus2transcript_dic[symbol].append(related_accession)
            if product_accession != '':
                transcript2_sympol[product_accession] = symbol
                transcript2chr_dic[product_accession] = chromosome
                transcript2anno_dic[product_accession] = name
                transcript2_GeneID[product_accession] = GeneID
            if related_accession != '':
                transcript2chr_dic[related_accession] = chromosome
                transcript2anno_dic[related_accession] = name
                transcript2_GeneID[related_accession] = GeneID
                transcript2_sympol[related_accession] = symbol
    
    print(f'Loading integration.o')
    sRNA_dic = nestedDic()
    with open(integration_o, 'r') as fn:
        for line in fn:
            line1 = line.strip()
            if line1.startswith('#') or line1.startswith('>'):
                continue
            l = line.strip().split("\t")
            gene = l[0]
            pos = l[2]
            marker = l[-1]
            seq = l[5]
            sRNA_dic[marker+';'+seq] = pos
    
    # 处理输出文件
    with open(PHAS_Loci, 'r') as fn:
        for line in fn:
            line1 = line.strip()
            l = line.strip().split("\t")
            if l[0] == 'cluster_number':
                continue
            cluster_number,feature,geneid,method,reference_type,h_transcriptome_coordinate,h_genome_coordinate,h_fl_transcriptome_coordinate,p_transcriptome_coordinate,p_genome_coordinate,p_fl_transcriptome_coordinate,c_pvalue,g_pvalue,f_pvalue,c_phaseratio,g_phaseratio,f_phaseratio,c_phasescore,g_phasescore,f_phasescore = \
            int(l[0]), l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10], l[11], l[12], l[13], l[14], l[15], l[16], l[17], l[18], l[19]

            t_starts = []
            t_ends = []
            g_starts = []
            g_ends = []
            pvalues = []
            g_pvalues = []
            phaseScores = []
            g_phaseScores = []
            # 按照HC, PC, HF, PF, HF, PG的顺序使用坐标，pvalue和cutoff
            if 'C' in reference_type:
                if 'H' in method:
                    for coor in h_transcriptome_coordinate.split(';'):
                        coor1 = coor.replace('(', '').replace(')', '').split(':')
                        t_start = coor1[0]
                        t_end = coor1[-1]
                        pvalues.append(c_pvalue)
                        phaseScores.append(c_phasescore)
                        t_starts.append(t_start)
                        t_ends.append(t_end)
                        Record = "HC"
                elif 'P' in method:
                    for coor in p_transcriptome_coordinate.split(';'):
                        coor1 = coor.replace('(', '').replace(')', '').split(':')
                        pvalues.append(c_pvalue)
                        phaseScores.append(c_phasescore)
                        t_start = coor1[0]
                        t_end = coor1[-1]
                        t_starts.append(t_start)
                        t_ends.append(t_end)
                        Record = "PC"
            elif 'F' in reference_type:
                if 'H' in method:
                    for coor in h_fl_transcriptome_coordinate.split(';'):
                        coor1 = coor.replace('(', '').replace(')', '').split(':')
                        t_start = coor1[0]
                        t_end = coor1[-1]
                        pvalues.append(f_pvalue)
                        phaseScores.append(f_phasescore)
                        t_starts.append(t_start)
                        t_ends.append(t_end)
                        Record = "HF"
                elif 'P' in method:
                    for coor in p_fl_transcriptome_coordinate.split(';'):
                        coor1 = coor.replace('(', '').replace(')', '').split(':')
                        pvalues.append(f_pvalue)
                        phaseScores.append(f_phasescore)
                        t_start = coor1[0]
                        t_end = coor1[-1]
                        t_starts.append(t_start)
                        t_ends.append(t_end)
                        Record = "PF"

            elif 'G' in reference_type:
                if 'H' in method:
                    for coor in h_genome_coordinate.split(';'):
                        coor1 = coor.replace('(', '').replace(')', '').split(':')
                        g_start = coor1[0]
                        g_end = coor1[-1]
                        if len(pvalues) == 0:
                            g_pvalues.append(g_pvalue)
                        if len(phaseScores) == 0:
                            g_phaseScores.append(g_phasescore)
                        g_starts.append(g_start)
                        g_ends.append(g_end)
                        Record = "HG"
                elif 'P' in method:
                    for coor in p_genome_coordinate.split(';'):
                        coor1 = coor.replace('(', '').replace(')', '').split(':')
                        g_start = coor1[0]
                        g_end = coor1[-1]
                        g_starts.append(g_start)
                        g_ends.append(g_end)
                        if len(pvalues) == 0:
                            g_pvalues.append(g_pvalue)
                        if len(phaseScores) == 0:
                            g_phaseScores.append(g_phasescore)
                        Record = "PG"

            chromosome = transcript2chr_dic[geneid]
            if chromosome == {}:
                chromosome = "NULL"
            Locus_Gene = transcript2_sympol[geneid]
            transcript_id = transcript2_GeneID[geneid]
            if transcript_id == {}:
                transcript_id = geneid


            if len(pvalues) == 0 and len(g_pvalues) != 0:
                out_pvalues = g_pvalues
            else:
                out_pvalues = pvalues
            if len(phaseScores) == 0 and len(g_phaseScores) != 0:
                out_phaseScores = g_phaseScores
            else:
                out_phaseScores = phaseScores 

            out_pvalues = set(out_pvalues)
            out_phaseScores = set(out_phaseScores)

            Anno = transcript2anno_dic[geneid]
            if Anno == {}:
                Anno = "NULL"
            PHAS_fa_query = feature + ':' + geneid
            Plot = f'{geneid}_{feature}:{geneid}_{Record}_plot.jpg'
            PHAS_Gene = f'{spe}PHAS{chromosome}g{transcript_id}@{geneid}@{str(phase_length)}\t{chromosome}\t{",".join(g_starts)}\t{",".join(g_ends)}\t{Locus_Gene}\t{geneid}\
                \t{feature}\t{",".join(t_starts)}\t{",".join(t_ends)}\t{",".join(out_pvalues)}\t{",".join(out_phaseScores)}\t{Anno}\t{method}\t{reference_type}\t{Record}\t{Plot}'
            PHAS_txt.write(PHAS_Gene + "\n")
            out_PHAS_Gene = f'{spe}PHAS{chromosome}g{transcript_id}@{geneid}@{str(phase_length)}'
            if geneid in ref_dic:
                PHAS_fa_for_initiator.write(f'>{geneid}\n{ref_dic[geneid]}\n')
            for descr in PHAS_dic:
                cnm = descr.split(' ')
                feature_gene = cnm[1]    
                if PHAS_fa_query in feature_gene:
                    seq = PHAS_dic[descr]
                    PHAS_fa.write(f'>{spe}PHAS{chromosome}g{transcript_id}@{geneid}@{str(phase_length)}\n{seq}' + "\n")

            count = 0
            # 处理phasiRNA
            for descr in phasiRNA_dic:
                if geneid in descr:
                    cnm = descr.split(' ')
                    strand = cnm[-1]
                    marker = cnm[-4]
                    seq = phasiRNA_dic[descr]
                    length = len(phasiRNA_dic[descr])
                    count += 1
                    pos = sRNA_dic[marker + ";" + seq]
                    # 直接使用marker和序列来识别
                    pos_s = pos
                    g_start = 'NA'
                    g_end = 'NA'
                    pos_e = int(pos) + phase_length
                    phasiRNA_fa_out = f'>{spe}PHAS{chromosome}g{transcript_id}@{geneid}@{str(phase_length)}.{str(count)}\n{seq}'
                    phasiRNA_id = f'{spe}PHAS{chromosome}g{transcript_id}@{geneid}@{str(phase_length)}.{str(count)}'
                    phasiRNA_fa.write(phasiRNA_fa_out + "\n")
                    if phasiRNA_id in phasiRNA_deg_dic1:
                        verified = 'yes'
                        verified_Target_num = len(set(phasiRNA_deg_dic1[phasiRNA_id]))
                        for target in set(phasiRNA_deg_dic1[phasiRNA_id]):
                            for j in phasiRNA_deg_dic2[phasiRNA_id + '\t' + target]:
                                verified_Target = target
                                degradome = j[0]
                                category = j[1]
                                sRNA_loc = j[2]
                                Deg_loc = j[3]
                                deg_count = j[4]
                                shift = j[5]
                                T_plot = degradome + '_' + j[6]
                                phasiRNA_out = f'{spe}PHAS{chromosome}g{transcript_id}@{geneid}@{str(phase_length)}.{str(count)}\t{chromosome}\t{g_start}\t{g_end}\t{strand}\t{seq}\t{length}\t{out_PHAS_Gene}\t{geneid}\t{pos_s}\t{pos_e}\t{verified}\t{verified_Target_num}\t{verified_Target}\t{degradome}\t{category}\t{sRNA_loc}\t{Deg_loc}\t{deg_count}\t{shift}\t{T_plot}'
                                phasiRNA_txt.write(phasiRNA_out + "\n")
                    else:
                        verified = 'no'
                        verified_Target = '-'
                        verified_Target_num = 0
                        degradome = '-'
                        category = '-'
                        sRNA_loc = '-'
                        Deg_loc = '-'
                        deg_count = '-'
                        shift = '-'
                        T_plot = '-'
                        phasiRNA_out = f'{spe}PHAS{chromosome}g{transcript_id}@{geneid}@{str(phase_length)}.{str(count)}\t{chromosome}\t{g_start}\t{g_end}\t{strand}\t{seq}\t{length}\t{out_PHAS_Gene}\t{geneid}\t{pos_s}\t{pos_e}\t{verified}\t{verified_Target_num}\t{verified_Target}\t{degradome}\t{category}\t{sRNA_loc}\t{Deg_loc}\t{deg_count}\t{shift}\t{T_plot}'
                        phasiRNA_txt.write(phasiRNA_out + "\n")


            # 处理initiator
            with open(initiator_prediction, 'r') as fn:
                for line in fn:
                    line1 = line.strip()
                    l = line.strip().split("\t")
                    miRNA, Gid, cleavage_site, target_pos_end, posmax, distance, initiator_seq, initiator_len, recordermaker = \
                    l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8]
                    if miRNA + '@' + Gid not in initiator_verified_dic:
                        verified = 'No'
                        degradome = '-'
                        category = '-'
                        deg_count = '-'
                        T_plot = '-'
                        if Gid in out_PHAS_Gene:
                            ini_out.write(f'{out_PHAS_Gene}\t{miRNA}\t{Gid}\t{cleavage_site}\t{target_pos_end}\t{posmax}\t{distance}\t{initiator_seq}\t{initiator_len}\t{recordermaker}\t{verified}\t{degradome}\t{category}\t{deg_count}\t{T_plot}'+ "\n")
                    else:
                        verified = 'Yes'
                        tmp = initiator_verified_dic[miRNA + '@' + Gid]
                        for i in tmp:
                            degradome = i[1]
                            category = i[0]
                            deg_count = i[2]
                            T_plot = degradome + '_' + i[3]
                            if Gid in out_PHAS_Gene:
                                ini_out.write(f'{out_PHAS_Gene}\t{miRNA}\t{Gid}\t{cleavage_site}\t{target_pos_end}\t{posmax}\t{distance}\t{initiator_seq}\t{initiator_len}\t{recordermaker}\t{verified}\t{degradome}\t{category}\t{deg_count}\t{T_plot}' + "\n")

def main():
    phasiRNAFormat()

if __name__ == "__main__":
    main()