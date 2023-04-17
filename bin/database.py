import sys
from Bio import SeqIO
import os
from customeDic import *
spe = 'Zm'
phase_length = 21
help = '''
    phasiHunter database <PHAS_Loci_database> <feature_table> <database_phasiRNA.fa> <all.integtation.o> <database_PHAS.fa> <refcdna.fa>
    # <PHAS_Gene.fa> <initiator.info> <target.info>
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
# PHAS_Gene = sys.argv[4]
# initiator =sys.argv[5]
# target = sys.argv[6]

TWD = os.getcwd()

PHAS_txt = open(f'{TWD}/PHAS.txt', 'w')
phasiRNA_txt = open(f'{TWD}/phasiRNA.txt', 'w')
phasiRNA_fa = open(f'{TWD}/phasiRNA.fa', 'w')
PHAS_fa = open(f'{TWD}/PHAS.fa', 'w')
PHAS_fa_for_initiator = open(f'{TWD}/forInitiatorPHAS.fa', 'w')

def phasiRNAFormat():
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
    PHAS_txt.write(f"PHAS_Gene\tChromosome\tStart\tEnd\tLocus_Gene\tTranscript_id\tPHAS_Gene_Type\tStart_in_transcript\tEnd_in_transcript\tP-value\tPhasiRNA_num\tTarget_verified\tInitiator_num\tInitiator_verified\tCluster\tAnnotation" + "\n")
    phasiRNA_txt.write(f'PhasiRNA\tChromosome\tStart\tEnd\tStrand\tSequence\tLength\tSource\tPHAS_Gene\tTranscript_id\tPosition_in_PHAS_Gene\tStrand_in_PHAS_Gene\tCluster\tVerified\tVerified_Target_num\tVerified_Target\tDegradome\tCategory\tSmallRNA_loc\tDeg_loc\tDeg_Count\tShift\tT_plot')
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
            if 'C' in reference_type:
                if 'H' in method:
                    for coor in h_transcriptome_coordinate.split(';'):
                        coor1 = coor.replace('(', '').replace(')', '').split(':')
                        t_start = coor1[0]
                        t_end = coor1[-1]
                        pvalues.append(c_pvalue)
                        t_starts.append(t_start)
                        t_ends.append(t_end)
                elif 'P' in method:
                    for coor in p_transcriptome_coordinate.split(';'):
                        coor1 = coor.replace('(', '').replace(')', '').split(':')
                        t_start = coor1[0]
                        t_end = coor1[-1]
                        t_starts.append(t_start)
                        t_ends.append(t_end)
            elif 'F' in reference_type:
                if 'H' in method:
                    for coor in h_fl_transcriptome_coordinate.split(';'):
                        coor1 = coor.replace('(', '').replace(')', '').split(':')
                        t_start = coor1[0]
                        t_end = coor1[-1]
                        pvalues.append(f_pvalue)
                        t_starts.append(t_start)
                        t_ends.append(t_end)
                elif 'P' in method:
                    for coor in p_fl_transcriptome_coordinate.split(';'):
                        coor1 = coor.replace('(', '').replace(')', '').split(':')
                        t_start = coor1[0]
                        t_end = coor1[-1]
                        t_starts.append(t_start)
                        t_ends.append(t_end)

            if 'G' in reference_type:
                if 'H' in method:
                    for coor in h_genome_coordinate.split(';'):
                        coor1 = coor.replace('(', '').replace(')', '').split(':')
                        g_start = coor1[0]
                        g_end = coor1[-1]
                        if len(pvalues) == 0:
                            g_pvalues.append(g_pvalue)
                        g_starts.append(g_start)
                        g_ends.append(g_end)
                elif 'P' in method:
                    for coor in p_genome_coordinate.split(';'):
                        coor1 = coor.replace('(', '').replace(')', '').split(':')
                        g_start = coor1[0]
                        g_end = coor1[-1]
                        g_starts.append(g_start)
                        g_ends.append(g_end)

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

            out_pvalues = set(out_pvalues)
            Target_verified = ''
            Initiator_num = ''
            Initiator_verified = ''
            Anno = transcript2anno_dic[geneid]
            if Anno == {}:
                Anno = "NULL"
            PHAS_fa_query = feature + ':' + geneid
            PHAS_Gene = f'{spe}PHAS{chromosome}g{transcript_id}@{geneid}@{str(phase_length)}\t{chromosome}\t{",".join(g_starts)}\t{",".join(g_ends)}\t{Locus_Gene}\t{geneid}\
                \t{",".join(t_starts)}\t{",".join(t_ends)}\t{",".join(out_pvalues)}\t{method}\t{reference_type}\t{Target_verified}\t{Initiator_num}\t{Initiator_verified}\t{Anno}'
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
            for descr in phasiRNA_dic:
                if geneid in descr:
                    cnm = descr.split(' ')
                    strand = cnm[-1]
                    marker = cnm[-4]
                    seq = phasiRNA_dic[descr]
                    length = len(phasiRNA_dic[descr])
                    count += 1
                    pos = sRNA_dic[marker + ";" + seq]
                    pos_s = pos
                    g_start = 'NA'
                    g_end = 'NA'
                    pos_e = int(pos) + phase_length
                    # phasiRNA_out = f'{spe}PHAS{chromosome}g{transcript_id}@{str(phase_length)}.{str(count)}\t{chromosome}\t{g_start}\t{g_end}\t{strand}\t{seq}\t{length}\t{out_PHAS_Gene}\t{geneid}\t{pos_s}\t{pos_e}\t{Verifed}\{V_Target_num}\t{V_Target}\t{Degradome}\t{Category}\t{SmallRNA_loc}\t{Deg_loc}\t{Deg_Count}\t{Shift}\t{T-plot}'
                    phasiRNA_out = f'{spe}PHAS{chromosome}g{transcript_id}@{geneid}@{str(phase_length)}.{str(count)}\t{chromosome}\t{g_start}\t{g_end}\t{strand}\t{seq}\t{length}\t{out_PHAS_Gene}\t{geneid}\t{pos_s}\t{pos_e}'
                    phasiRNA_txt.write(phasiRNA_out + "\n")
                    phasiRNA_fa_out = f'>{spe}PHAS{chromosome}g{transcript_id}@{geneid}@{str(phase_length)}.{str(count)}\n{seq}'
                    phasiRNA_fa.write(phasiRNA_fa_out + "\n")

def main():
    phasiRNAFormat()

if __name__ == "__main__":
    main()