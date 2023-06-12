# %%
import sys
from Bio import SeqIO
import os
from customeDic import *
# feature_table = sys.argv[2]
# feature_table = '/home/riceUsers/fzr/refseq/GCF_000001735.4_TAIR10.1_feature_table.txt'
# PHAS_Loci_file = '/data1/fzr_data/thesis/ath/all.integration.phas'
# phasiRNA_file = '/data1/fzr_data/thesis/ath/phasiRNA.fa'
# PHAS_file = '/data1/fzr_data/thesis/ath/PHAS.fa'
# Spe = 'Ath'
# phase_length = 21

help = '<Spe> <feature_table> <PHAS_Loci file> <phasiRNA file> <PHAS fa file> <integration_file> <phase_length> <out_dir> <all.integration.s>'
version = ''
for i in range(1, len(sys.argv)):
    if sys.argv[i] == '-v':
        print(version)
        sys.exit()
    elif sys.argv[i] == '-h':
        print(help)
        sys.exit()

Spe = sys.argv[1]
feature_table = sys.argv[2]
PHAS_Loci_file = sys.argv[3]
phasiRNA_file = sys.argv[4]
PHAS_file = sys.argv[5]
integration_file = sys.argv[6]
phase_length = int(sys.argv[7])
out_dir = sys.argv[8]
integration_summary = sys.argv[9]

# out_dir = './'
new_PHAS_Loci_file = out_dir + '/new_PHAS_Loci_file.txt'
new_phasiRNA_fa_file = out_dir + '/new_phasiRNA.fa'
new_PHAS_Gene_fa_file = out_dir + '/new_PHAS_Gene.fa'
new_integration_file = out_dir + '/new_integration.o'
new_integration_summary = out_dir + '/new_integration.s'
new_phasiRNA_txt = out_dir + '/new_phasiRNA.txt'
fo_integration_s = open(new_integration_summary, 'w')
fo_integration = open(new_integration_file, 'w')
fo_PHAS_Loci = open(new_PHAS_Loci_file, 'w')
fo_phasiRNA_fa = open(new_phasiRNA_fa_file, 'w')
fo_PHAS_GENE_fa = open(new_PHAS_Gene_fa_file, 'w')
fo_phasiRNA_txt = open(new_phasiRNA_txt, 'w')
# Spe = sys.argv[1]

print('Loading feature table')
locus2transcript_dic = OneDepDic()
transcript2_sympol = OneDepDic()
transcript2chr_dic = nestedDic()
transcript2anno_dic = nestedDic()
transcript2_GeneID = nestedDic()
chr2chr = nestedDic()
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
        if genomic_accession != '':
            chr2chr[genomic_accession] = chromosome


print('Loading PHAS_Loci information')
with open(PHAS_Loci_file, 'r') as fn:
    relist = []
    filter_list = []
    PHAS_gene2phasiRNA_dic = nestedDic()
    PHAS_gene2chromosome_dic = nestedDic()
    PHAS_gene2transcript = nestedDic()
    order_dic = OneDepDic()
    remember_order = nestedDic()
    cluster_order_dic = nestedDic()
    for line in fn:
        line1 = line.strip()
        l = line.strip().split("\t")
        feature,geneid,Genome_start,Genome_end,transcript_start,transcript_end,method_ref,pvalue,phase_score,phase_ratio,tag,recorder = l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8],l[9],l[10],l[11]
        if method_ref == 'PG' or method_ref == 'PC' or method_ref == "PF" or method_ref == "PC,PG" or method_ref == "PG,PC":
            relist.append(line)
            continue
        if line1.startswith('feature'):
            fo_PHAS_Loci.write('PHAS_Loci\tChromosome\tGenome_start\tGenome_end\tLocus_Gene\tTranscript\tPHAS_Gene_Type\tTranscript_start\tTranscript_end\tPvalue\tPhaseScore\tPhaseRatio\tAnnotation\tMethod_Ref\tTag\trecorder\n')
            continue
        filter_list.append(tag)
        chromosome = transcript2chr_dic[geneid]
        if len(chromosome) == 0:
            if feature == 'Intergenic':
                chromosome = chr2chr[geneid.split('#')[0]]
            else:
                chromosome = 'n'
        if type(chromosome) == int:
            if int(chromosome) < 10:
                chromosome = '0' + chromosome

        transcript_id = transcript2_GeneID[geneid]
        if len(transcript_id) == 0:
            if feature == "Intergenic":
                transcript_id = 'Intergenic'
                geneid = geneid.split('#')[0]
            elif 'novel_in_catalog' in feature:
                transcript_id = 'NIC_' + feature.split('_')[-2]
            elif 'novel_not_in_catalog' in feature:
                transcript_id = 'NNC_' + feature.split('_')[-2]
            else:
                transcript_id = feature + '_' + geneid
        
        order_dic[transcript_id].append(geneid)
        order = len(set(order_dic[transcript_id]))
        remember_order[geneid] = order
        
        
        if feature == 'Intergenic':
            cluster_order = cluster_order_dic[geneid.split('#')[0]]
        else:
            cluster_order = cluster_order_dic[geneid]
        if cluster_order == {}:
            cluster_order_dic[geneid.split('#')[0]] = 1
            cluster_order = cluster_order_dic[geneid.split('#')[0]]

        if feature == 'Intergenic':
            order = 1
        PHAS_Loci_name = f'{Spe}PHAS{chromosome}g{transcript_id}.{order}@{str(phase_length)}#{cluster_order}'
        cluster_order_dic[geneid.split('#')[0]] += 1
        if len(transcript2_GeneID[geneid]) == 0:
            pass
        else:
            transcript_id = 'LOC' + transcript_id
        PHAS_gene2phasiRNA_dic[recorder] = PHAS_Loci_name
        PHAS_gene2chromosome_dic[recorder] = chromosome
        PHAS_gene2transcript[recorder] = geneid
        annotation = transcript2anno_dic[geneid]
        if annotation == {}:
            annotation = 'None'
        # fo_PHAS_Loci.write(l[0] + '\t' + PHAS_Loci_name + '\t' + geneid + '\t' + annotation + '\t' + '\t'.join(l[2:]) + '\n')
        fo_PHAS_Loci.write(PHAS_Loci_name + '\t' + chromosome + '\t' + Genome_start + '\t' + Genome_end + '\t' + transcript_id + '\t' + geneid + '\t' + feature + '\t' + transcript_start + '\t' + transcript_end + '\t' + pvalue + '\t' + phase_score + '\t' + phase_ratio + '\t' + annotation + '\t' + method_ref +  '\t' + tag + '\t' + recorder + '\n')

    for line in relist:
        line1 = line.strip()
        l = line.strip().split("\t")
        feature,geneid,Genome_start,Genome_end,transcript_start,transcript_end,method_ref,pvalue,phase_score,phase_ratio,tag,recorder = l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8],l[9],l[10],l[11]
        filter_list.append(tag)
        chromosome = transcript2chr_dic[geneid]
        if len(chromosome) == 0:
            if feature == 'Intergenic':
                chromosome = chr2chr[geneid.split('#')[0]]
            else:
                chromosome = 'n'
        if type(chromosome) == int:
            if int(chromosome) < 10:
                chromosome = '0' + chromosome

        transcript_id = transcript2_GeneID[geneid]
        if len(transcript_id) == 0:
            if feature == "Intergenic":
                transcript_id = 'Intergenic'
                geneid = geneid.split('#')[0]
            elif 'novel_in_catalog' in feature:
                transcript_id = 'NIC_' + feature.split('_')[-2]
            elif 'novel_not_in_catalog' in feature:
                transcript_id = 'NNC_' + feature.split('_')[-2]
            else:
                transcript_id = feature + '_' + geneid
        
        order_dic[transcript_id].append(geneid)
        order = len(set(order_dic[transcript_id]))
        if geneid in remember_order:
            order = remember_order[geneid]
        
        if feature == 'Intergenic':
            cluster_order = cluster_order_dic[geneid.split('#')[0]]
        else:
            cluster_order = cluster_order_dic[geneid]
        if cluster_order == {}:
            cluster_order_dic[geneid.split('#')[0]] = 1
            cluster_order = cluster_order_dic[geneid.split('#')[0]]

        if feature == 'Intergenic':
            order = 1
        PHAS_Loci_name = f'{Spe}PHAS{chromosome}g{transcript_id}.{order}@{str(phase_length)}#{cluster_order}'
        if len(transcript2_GeneID[geneid]) == 0:
            pass
        else:
            transcript_id = 'LOC' + transcript_id
        cluster_order_dic[geneid.split('#')[0]] += 1
        PHAS_gene2phasiRNA_dic[recorder] = PHAS_Loci_name
        PHAS_gene2chromosome_dic[recorder] = chromosome
        PHAS_gene2transcript[recorder] = geneid
        annotation = transcript2anno_dic[geneid]
        if annotation == {}:
            annotation = 'None'
        # fo_PHAS_Loci.write(l[0] + '\t' + PHAS_Loci_name + '\t' + geneid + '\t' + annotation + '\t' + '\t'.join(l[2:]) + '\n')
        fo_PHAS_Loci.write(PHAS_Loci_name + '\t' + chromosome + '\t' + Genome_start + '\t' + Genome_end + '\t' + transcript_id + '\t' + geneid + '\t' + feature + '\t' + transcript_start + '\t' + transcript_end + '\t' + pvalue + '\t' + phase_score + '\t' + phase_ratio + '\t' + annotation + '\t' + method_ref +  '\t' + tag + '\t' + recorder + '\n')




for query in SeqIO.parse(phasiRNA_file, 'fasta'):
    id = query.description
    seq = str(query.seq)
    l = id.split('__')
    tag, order = l[0], l[-1]
    new_phasiRNA_name = PHAS_gene2phasiRNA_dic[tag] + '.' + order
    fo_phasiRNA_fa.write('>'+new_phasiRNA_name+'\n'+seq+'\n')

PHAS_GENE_name2integration = OneDepDic()
for query in SeqIO.parse(PHAS_file, 'fasta'):
    id = query.description
    seq = str(query.seq)
    l = id.split('__')
    tag, order = l[0], l[-1]
    gene = l[1]
    new_PHAS_Gene_name = PHAS_gene2phasiRNA_dic[tag]
    fo_PHAS_GENE_fa.write('>'+new_PHAS_Gene_name+'\n'+seq+'\n')
    PHAS_GENE_name2integration[gene + '\t' + order].append(new_PHAS_Gene_name)

with open(integration_file, 'r') as fn:
    for line in fn:
        line1 = line.strip()
        l = line.strip().split("\t")
        if line1.startswith('>'):
            continue
        if line1.startswith('#'):
            continue
        recorder = l[-1]
        gene = l[0]
        if recorder not in filter_list:
            continue
        tmps = PHAS_GENE_name2integration[gene + '\t' + recorder]
        for tmp in tmps:
            fo_integration.write(tmp + '\t' + '\t'.join(l[1:]) + '\n')


with open(integration_summary, 'r') as fn:
    for line in fn:
        line1 = line.strip()
        l = line.strip().split("\t")
        cluster_number,feature,geneid,method,reference_type,h_transcriptome_coordinate,h_genome_coordinate,h_fl_transcriptome_coordinate,p_transcriptome_coordinate,p_genome_coordinate,p_fl_transcriptome_coordinate,c_pvalue,g_pvalue,f_pvalue,c_phaseratio,g_phaseratio,f_phaseratio,c_phasescore,g_phasescore,f_phasescore = \
        l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8],l[9],l[10],l[11],l[12],l[13],l[14],l[15],l[16],l[17],l[18],l[19]
        if line1.startswith('cluster_number'):
            fo_integration_s.write(f'Locus\t{line1}\tAnnotation\n')
            continue
        annotation = transcript2anno_dic[geneid]
        transcript_id = transcript2_sympol[geneid]
        if len(transcript_id) == 0:
            if feature == "Intergenic":
                transcript_id = 'Intergenic'
                geneid = geneid.split('#')[0]
            elif 'novel_in_catalog' in feature:
                transcript_id = 'NIC_' + feature.split('_')[-2]
            elif 'novel_not_in_catalog' in feature:
                transcript_id = 'NNC_' + feature.split('_')[-2]
            else:
                transcript_id = feature + '_' + geneid
        fo_integration_s.write(f'{transcript_id}\t{line1}\t{annotation}\n')


fo_phasiRNA_txt.write(f'PhasiRNA\tChromosome\tGenome_Start\tGenome_End\tStrand\tSequence\tLength\tPHAS_Loci\tTranscript_id\tStart_in_Transcript\tEnd_in_Transcript\n')
for query in SeqIO.parse(phasiRNA_file, 'fasta'):
    id = query.description
    seq = str(query.seq)
    l = id.split('__')
    tag,chr_transcript,start,abun,strand,order = l[0],l[1],int(l[2]),l[3],l[4],l[5]
    phasiRNA = PHAS_gene2phasiRNA_dic[tag] + '.' + order
    chromosome = PHAS_gene2chromosome_dic[tag]
    seq_len = len(seq)
    PHAS_Loci = '.'.join(phasiRNA.split('.')[:-1])
    transcript = PHAS_gene2transcript[tag]
    if 'Intergenic' in phasiRNA:
        Genome_start = start
        Genome_end = start + len(seq)
        Transcript_start = '-'
        Transcript_end = '-'
    else:
        Genome_start = '-'
        Genome_end = '-'
        Transcript_start = start
        Transcript_end = start + len(seq)
    fo_phasiRNA_txt.write(f'{phasiRNA}\t{chromosome}\t{Genome_start}\t{Genome_end}\t{strand}\t{seq}\t{seq_len}\t{PHAS_Loci}\t{transcript}\t{Transcript_start}\t{Transcript_end}\n')

fo_PHAS_Loci.close()
fo_PHAS_GENE_fa.close()
fo_phasiRNA_fa.close()
fo_integration.close()
fo_integration_s.close()