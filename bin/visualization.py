# %%
from Bio import SeqIO
import matplotlib
from matplotlib import pyplot as plt
matplotlib.use('Agg')
from collections import defaultdict
import sys
import os
import time

def Load_PHAS_info(PHAS_Loci_file):
    feature_gene_dic = OneDepDic()
    list = []
    dic = OneDepDic()
    with open(PHAS_Loci_file, 'r') as fn:
        for line in fn:
            line1 = line.strip()
            if line1.startswith('feature'):
                continue
            l = line.strip().split("\t")
            feature,PHAS_Loci,Genome_start,Genome_end,transcript_start,transcript_end,method_ref,pvalue,phase_score,phase_ratio,tag, recorder = l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8],l[9],l[10],l[11]
            list.append(tag)
            dic[tag].append(recorder)
            if 'F' not in tag:
                feature_gene_dic[tag].append(feature + ':' + PHAS_Loci)
            else:
                feature_gene_dic[tag].append(PHAS_Loci)
    return (list, dic, feature_gene_dic)

def GetCurTime():
    """return formated current localtime 

    Returns:
        <str> -- formated current localtime
    """
    ctime = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
    return ctime

def CompleteStrand(x):
    x_string=''
    rc={"A":"T","T":"A","C":"G","G":"C","N":"N", "a":"t", "t":"a", "c":"g", "g":"c", "n":"n"}
    for i in range(len(x)):
        x_string+=rc[x[i]] 
    return x_string

def nestedDic():
    return defaultdict(nestedDic)

def OneDepDic():
    return defaultdict(list)

def TwoDepDic():
    return defaultdict(OneDepDic)

def ThreeDepDic():
    return defaultdict(TwoDepDic)

def FourDepDic():
    return defaultdict(ThreeDepDic)

def FiveDepDic():
    return defaultdict(FourDepDic)

def LoadsRNA(sRNA_file, filter_list):
    dic = ThreeDepDic()
    with open(sRNA_file, 'r') as fn:
        for line in fn:
            l = line.strip().split("\t")
            if line.startswith('>'):
                method = l[0]
                continue
            if line.startswith('#'):
                ref = "_".join(l)
                continue
            if method == '>Hypergeometric':
                geneid, strand, pos, abun, sRNAid, seq, seq_len, phaseRatio, phaseNumber, phaseAbun, phaseScore, pvalue, anno, reocrd_marker \
                    = l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10], l[11], l[12], l[13]
                if reocrd_marker in filter_list:
                    dic[method][ref][(geneid, pvalue, anno, reocrd_marker)].append((strand, int(pos), abun, seq))
            if method == '>PhaseScore':
                geneid, strand, pos, abun, sRNAid, seq, seq_len, phaseRatio, phaseNumber, phaseAbun, phaseScore, pvalue, anno, reocrd_marker \
                    = l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10], l[11], l[12], l[13]
                if reocrd_marker in filter_list:
                    dic[method][ref][(geneid, phaseScore, anno, reocrd_marker)].append((strand, int(pos), abun, seq))
    return dic

def LoadallsRNA(sRNA_file, filter_list):
    dic = FourDepDic()
    with open(sRNA_file, 'r') as fn:
        for line in fn:
            l = line.strip().split("\t")
            if line.startswith('>'):
                method = l[0]
                continue
            if line.startswith('#'):
                ref = "_".join(l)
                continue
            if method == '>Hypergeometric':
                geneid, strand, pos, abun, sRNAid, seq, seq_len, phaseRatio, phaseNumber, phaseAbun, phaseScore, pvalue, anno, reocrd_marker \
                    = l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10], l[11], l[12], l[13]
                if reocrd_marker in filter_list:
                    dic[method][ref][(geneid, pvalue, anno, reocrd_marker)][pos].append((strand, int(pos), abun, seq))
            if method == '>PhaseScore':
                geneid, strand, pos, abun, sRNAid, seq, seq_len, phaseRatio, phaseNumber, phaseAbun, phaseScore, pvalue, anno, reocrd_marker \
                    = l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10], l[11], l[12], l[13]
                if reocrd_marker in filter_list:
                    dic[method][ref][(geneid, phaseScore, anno, reocrd_marker)][pos].append((strand, int(pos), abun, seq))
    return dic

def LoadRef(ref_file):
    dic = nestedDic()
    if ref_file != '':
        for query in SeqIO.parse(ref_file, 'fasta'):
            name = str(query.id)
            seq = str(query.seq)
            dic[name] = seq
        return dic
    else:
        return {}

def Alignment(cdna_ref_dic, gdna_ref_dic, flnc_ref_dic,phasiRNA_dic, fo, phase_length, recorder_dic):
    transcript_PHAS_Gene = nestedDic()
    transcript_phasiRNA = TwoDepDic() 
    genome_PHAS_Gene = nestedDic()
    genome_phasiRNA = TwoDepDic()
    flnc_PHAS_Gene = nestedDic()
    flnc_phasiRNA = TwoDepDic()
    genome_PHAS_Gene
    for i in phasiRNA_dic:
        # fo.write(f'{i}\n')
        for j in phasiRNA_dic[i]:
            # fo.write(f'{j}\n')
            if 'cDNA' in j and cdna_ref_dic != {}:
                for k in phasiRNA_dic[i][j]:
                    tmp_recorder = k[3]
                    recorders = recorder_dic[tmp_recorder]
                    for recorder in recorders:
                        fo.write('>' + recorder + '\t' + "\t".join(k) + '\n')
                        start_end = []
                        for e in phasiRNA_dic[i][j][k]:
                            start_end.append(e[1])
                        region = sorted(start_end)
                        region1 = (region[0], region[-1] + phase_length)
                        # transcript_PHAS_Gene[k] = PHASGene_seq(cdna_ref_dic[k[0]], region1)
                        transcript_PHAS_Gene[k] = cdna_ref_dic[k[0]]
                        for e in phasiRNA_dic[i][j][k]:
                            if e[0] == "+":
                                fo.write(f"{'.'*(int(e[1])-int(region[0]))}{e[3]}{'.'*(int(region[-1])-int(e[1]))}" + "______" + k[0] + "\t" + str(e[1])+ "\n")
                                transcript_phasiRNA[k]['+'].append(e)
                        seq = cdna_ref_dic[k[0]][int(region[0]):int(region[-1])+phase_length]
                        fo.write(seq + '______excised_seq' "\n")
                        fo.write(CompleteStrand(seq) + '______anti_seq' + "\n")
                        for e in phasiRNA_dic[i][j][k]:
                            if e[0] == '-':
                                fo.write(f"{'.'*(int(e[1])-int(region[0]))}{CompleteStrand(e[3])}{'.'*(int(region[-1])-int(e[1]))}" + "______" + k[0] + "\t" + str(e[1])+ "\n")
                                transcript_phasiRNA[k]['-'].append(e)
            if 'gDNA' in j and gdna_ref_dic != {}:
                for k in phasiRNA_dic[i][j]:
                    tmp_recorder = k[3]
                    recorders = recorder_dic[tmp_recorder]
                    for recorder in recorders:
                        fo.write('>' + recorder + '\t' + "\t".join(k) + '\n')
                        start_end = []
                        for e in phasiRNA_dic[i][j][k]:
                            start_end.append(e[1])
                        region = sorted(start_end)
                        region1 = (region[0], region[-1] + phase_length)
                        genome_PHAS_Gene[k] = PHASGene_seq(gdna_ref_dic[k[0]], region1)
                        for e in phasiRNA_dic[i][j][k]:
                            if e[0] == "+":
                                fo.write(f"{'.'*(int(e[1])-int(region[0]))}{e[3]}{'.'*(int(region[-1])-int(e[1]))}" + "______" + k[0] + "\t" + str(e[1])+ "\n")
                                genome_phasiRNA[k]['+'].append(e)
                        seq = gdna_ref_dic[k[0]][int(region[0]):int(region[-1])+phase_length]
                        fo.write(seq + '______excised_seq' "\n")
                        fo.write(CompleteStrand(seq) + '______anti_seq' + "\n")
                        for e in phasiRNA_dic[i][j][k]:
                            if e[0] == '-':
                                fo.write(f"{'.'*(int(e[1])-int(region[0]))}{CompleteStrand(e[3])}{'.'*(int(region[-1])-int(e[1]))}" + "______" + k[0] + "\t" + str(e[1])+ "\n")
                                genome_phasiRNA[k]['-'].append(e)
            if 'FLNC' in j and flnc_ref_dic != {}:
                for k in phasiRNA_dic[i][j]:
                    tmp_recorder = k[3]
                    recorders = recorder_dic[tmp_recorder]
                    for recorder in recorders:
                        fo.write('>' + recorder + '\t' + "\t".join(k) + '\n')
                        start_end = []
                        for e in phasiRNA_dic[i][j][k]:
                            start_end.append(e[1])
                        region = sorted(start_end)
                        region1 = (region[0], region[-1] + phase_length)
                        # flnc_PHAS_Gene[k] = PHASGene_seq(flnc_ref_dic[k[0]], region1)
                        flnc_PHAS_Gene[k] = flnc_ref_dic[k[0]]
                        for e in phasiRNA_dic[i][j][k]:
                            if e[0] == "+":
                                fo.write(f"{'.'*(int(e[1])-int(region[0]))}{e[3]}{'.'*(int(region[-1])-int(e[1]))}" + "______" + k[0] + "\t" + str(e[1])+ "\n")
                                flnc_phasiRNA[k]['+'].append(e)
                        seq = flnc_ref_dic[k[0]][int(region[0]):int(region[-1])+phase_length]
                        fo.write(seq + '______excised_seq' "\n")
                        fo.write(CompleteStrand(seq) + '______anti_seq' + "\n")
                        for e in phasiRNA_dic[i][j][k]:
                            if e[0] == '-':
                                fo.write(f"{'.'*(int(e[1])-int(region[0]))}{CompleteStrand(e[3])}{'.'*(int(region[-1])-int(e[1]))}" + "______" + k[0] + "\t" + str(e[1])+ "\n")
                                flnc_phasiRNA[k]['-'].append(e)
    return transcript_PHAS_Gene, transcript_phasiRNA, genome_PHAS_Gene, genome_phasiRNA, flnc_PHAS_Gene, flnc_phasiRNA

def PHASGene_seq(seq, region1):
    start = region1[0]
    end = region1[-1]
    for num in [1000, 900, 800, 700, 600, 500, 400, 300, 200, 100, 50, 0]:
        if start - num <= 0:
            continue
        else:
            start = start - num
            break
    for num in [1000, 900, 800, 700, 600, 500, 400, 300, 200, 100, 50, 0]:
        if end + num <= len(seq):
            continue
        else:
            end = end + num
            break
    return [seq[start :end], region1, (start, end)]
    

def Plot(phasiRNA_dic, allsiRNA_dic, phase_length, outdir, max_abun, cdna_ref, gdna_ref, flnc_ref, plot_c, plot_g, plot_f, recorder_dic, feature_gene_dic):
    count = 1
    for method in phasiRNA_dic:
        for ref in phasiRNA_dic[method]:
            for feature in phasiRNA_dic[method][ref]:
                if cdna_ref == {} or plot_c == 'n':
                    if 'PC' in feature[-1] or 'HC' in feature[-1]:
                        continue
                if gdna_ref == {} or plot_g == 'n':
                    if 'PG' in feature[-1] or 'HG' in feature[-1]:
                        continue
                if flnc_ref == {} or plot_f == 'n':
                    if 'PF' in feature[-1] or "HF" in feature[-1]:
                        continue
                recorders = recorder_dic[feature[3]]
                for recorder in recorders:
                    plt.figure(figsize=(8,8), dpi=300)
                    xtck = []
                    x=[] # + all_pos
                    y=[] # + all_read
                    x1=[] # - phas_pos
                    y1=[] # - phas_read
                    y2=[] # - all_read
                    x2=[] # - all_pos
                    x3=[] # + phas_pos
                    y3=[] # + phas_read
                    x4 = []
                    x5 = []
                    y4 = []
                    y5 = []
                    exsicedStrand_max_abun = max_abun
                    antiStrand_max_abun = max_abun
                    for pos_ in phasiRNA_dic[method][ref][feature]:
                        strand = pos_[0]
                        pos = pos_[1]
                        abun = float(pos_[2])
                        if strand == '+':
                            if exsicedStrand_max_abun < abun: abun = exsicedStrand_max_abun
                            x.append(pos)
                            y.append(abun)
                        if strand == '-':
                            if antiStrand_max_abun < abun: abun = antiStrand_max_abun
                            x1.append(pos)
                            y1.append(0 - abun)
                    for pos_ in allsiRNA_dic[method][ref][feature]:
                        sum_abun_plus = 0
                        sum_abun_minus = 0
                        strand_plus = 0
                        strand_minus = 0
                        for pos__ in allsiRNA_dic[method][ref][feature][pos_]:
                            strand = pos__[0]
                            pos = pos__[1]
                            abun = float(pos__[2])
                            if strand == '+':
                                sum_abun_plus += abun
                                strand_plus = 1
                            elif strand == '-':
                                sum_abun_minus += abun
                                strand_minus = 1
                        if strand_minus == 1:
                            if antiStrand_max_abun < sum_abun_minus: sum_abun_minus = antiStrand_max_abun
                            if pos in x1:
                                x5.append(pos)
                                y5.append(0 - sum_abun_minus)
                            else:
                                x2.append(pos)
                                y2.append(0 - sum_abun_minus)
                        if strand_plus == 1:
                            if exsicedStrand_max_abun < sum_abun_plus: sum_abun_plus = exsicedStrand_max_abun
                            if pos in x:
                                x4.append(pos)
                                y4.append(sum_abun_plus)
                            else:
                                x3.append(pos)
                                y3.append(sum_abun_plus)
                    x7 = x2 + x3 + x4 + x5
                    try:
                        max_pos = sorted(x7)[-1] + phase_length
                    except IndexError:
                        count += 1
                        # print(str(count))
                        print(count)
                        plt.close()
                        continue
                    min_pos = sorted(x7)[0]
                    for i in range(min_pos, max_pos, phase_length):
                        plt.axvline(i, lw=0.8, color='grey', ls='--')
                        xtck.append(i)
                    plt.vlines(x3,0,y3, color='c',linestyles='-', label='Non-phased reads')
                    plt.vlines(x2,y2,0, color='c',linestyles='-')
                    plt.vlines(x5,y5,0, color='red',label='Phased reads')
                    plt.vlines(x4,0,y4, color='red')
                    plt.xticks(xtck, rotation=315)
                    plt.axhline(y=0, color='grey')
                    # if 'PG' in feature[-1] or 'HG' in feature[-1]:
                    #     # tmp_feature = feature[2].split(';')
                    #     # tmp_list = []
                    #     # for fuck in tmp_feature:
                    #     #     if 'exon' in fuck or 'CDS' in fuck:
                    #     #         pass
                    #     #     else:
                    #     #         tmp_list.append(fuck)
                    #     if 'Intergenic' in feature[2]:
                    #         out_name = feature[0]+'_'+feature[2]+'_'+feature[-1]
                    #     else:
                    #         out_name = feature[0]+'_'+'Transcript'+'_'+feature[-1]
                    # else:
                    #     out_name = feature[0]+'_'+feature[2]+'_'+feature[-1]
                    out_name = recorder
                    plt.xlabel ("Nucleotide position") 
                    plt.ylabel ("Read Abundance \nReverse strand(-)/Forward strand(+)")
                    plt.legend (bbox_to_anchor=(0,1.01),loc='lower left', borderaxespad=0)
                    title1_list = feature_gene_dic[feature[3]]
                    if len(title1_list) > 1:
                        title1 = feature[0]
                    else:
                        title1 = title1_list[0]
                    # pvalue = str("%.3e" % float(feature[1]))
                    # plt.title(title1+'\npvalue='+pvalue, loc='right')
                    plt.title(title1, loc='right')
                    plt.savefig("{}/{}.jpg".format(outdir, out_name))
                    plt.close()

def FaOut(PHAS_Gene_phasiRNA_dic, fo_PHAS_Gene, fo_phasiRNA, recorder_dic):
    print('writing PHAS Gene')
    for gene in PHAS_Gene_phasiRNA_dic[0]:
        seq = PHAS_Gene_phasiRNA_dic[0][gene]
        # start = PHAS_Gene_phasiRNA_dic[0][gene][1][0]
        # end = PHAS_Gene_phasiRNA_dic[0][gene][1][1]
        # start1 = PHAS_Gene_phasiRNA_dic[0][gene][2][0]
        # end1 = PHAS_Gene_phasiRNA_dic[0][gene][2][1]
        recorders = recorder_dic[gene[3]]
        for recorder in recorders:
            fo_PHAS_Gene.write(f'>{recorder}__{gene[0]}__{gene[3]}\n')
            # fo_PHAS_Gene.write(f'>{gene[0]}\n')
            # fo_PHAS_Gene.write('>{} {}\t{} {}\t{} {}\t{}'.format(gene[0], gene[2], start, end, start1, end1, gene[-1]) + "\n")
            fo_PHAS_Gene.write(seq + '\n')
    for gene in PHAS_Gene_phasiRNA_dic[2]:
        seq = PHAS_Gene_phasiRNA_dic[2][gene][0]
        start = PHAS_Gene_phasiRNA_dic[2][gene][1][0]
        end = PHAS_Gene_phasiRNA_dic[2][gene][1][1]
        start1 = PHAS_Gene_phasiRNA_dic[2][gene][2][0]
        end1 = PHAS_Gene_phasiRNA_dic[2][gene][2][1]
        recorders = recorder_dic[gene[3]]
        for recorder in recorders:
            fo_PHAS_Gene.write('>{}__{}__{}__{}__{}__{}__{}'.format(recorder, gene[0], start, end, start1, end1, gene[3])+ '\n')
            # fo_PHAS_Gene.write('>{}'.format(gene[0])+ '\n')
            fo_PHAS_Gene.write(seq+ '\n')
    for gene in PHAS_Gene_phasiRNA_dic[4]:
        seq = PHAS_Gene_phasiRNA_dic[4][gene]
        # start = PHAS_Gene_phasiRNA_dic[4][gene][1][0]
        # end = PHAS_Gene_phasiRNA_dic[4][gene][1][1]
        # start1 = PHAS_Gene_phasiRNA_dic[4][gene][2][0]
        # end1 = PHAS_Gene_phasiRNA_dic[4][gene][2][1]
        recorders = recorder_dic[gene[3]]
        for recorder in recorders:
            fo_PHAS_Gene.write(f'>{recorder}__{gene[0]}__{gene[3]}\n')
            # fo_PHAS_Gene.write(f'>{gene[0]}\n')
            # fo_PHAS_Gene.write('>{} {}\t{} {}\t{} {}\t{}'.format(gene[0], gene[2], start, end, start1, end1, gene[-1])+ '\n')
            fo_PHAS_Gene.write(seq+ '\n')
    print('writing phasiRNA')
    for gene in PHAS_Gene_phasiRNA_dic[1]:
        recorders = recorder_dic[gene[3]]
        for recorder in recorders:
            count = 0
            for strand in PHAS_Gene_phasiRNA_dic[1][gene]:
                for phasiRNA in PHAS_Gene_phasiRNA_dic[1][gene][strand]:
                    count += 1
                    coor = phasiRNA[1]
                    abun = phasiRNA[2]
                    seq = phasiRNA[3]
                    fo_phasiRNA.write(f">{recorder}__{gene[0]}__{coor}__{abun}__{strand}__{count}\n")
                    fo_phasiRNA.write(f'{seq}\n')
                    # fo_phasiRNA.write('>{} {} {} phasiRNA {} {}'.format(gene[0], gene[2], gene[-1], count, strand) + "\n")
                    # fo_phasiRNA.write(phasiRNA+ '\n')
    for gene in PHAS_Gene_phasiRNA_dic[3]:
        recorders = recorder_dic[gene[3]]
        for recorder in recorders:
            out_list = []
            count = 0
            for strand in PHAS_Gene_phasiRNA_dic[3][gene]:
                for phasiRNA in PHAS_Gene_phasiRNA_dic[3][gene][strand]:
                    if phasiRNA in out_list:
                        continue
                    else:
                        out_list.append(phasiRNA)
                    count += 1
                    coor = phasiRNA[1]
                    abun = phasiRNA[2]
                    seq = phasiRNA[3]
                    fo_phasiRNA.write(f">{recorder}__{gene[0]}__{coor}__{abun}__{strand}__{count}\n")
                    fo_phasiRNA.write(f'{seq}\n')
                    # fo_phasiRNA.write('>{} {} {} phasiRNA {} {}'.format(gene[0], gene[2], gene[-1], count, strand) + "\n")
                    # fo_phasiRNA.write(phasiRNA+ '\n')
    for gene in PHAS_Gene_phasiRNA_dic[5]:
        recorders = recorder_dic[gene[3]]
        for recorder in recorders:
            count = 0
            for strand in PHAS_Gene_phasiRNA_dic[5][gene]:
                for phasiRNA in PHAS_Gene_phasiRNA_dic[5][gene][strand]:
                    count += 1
                    coor = phasiRNA[1]
                    abun = phasiRNA[2]
                    seq = phasiRNA[3]
                    fo_phasiRNA.write(f">{recorder}__{gene[0]}__{coor}__{abun}__{strand}__{count}\n")
                    fo_phasiRNA.write(f'{seq}\n')
                    # fo_phasiRNA.write('>{} {} {} phasiRNA {} {}'.format(gene[0], gene[2], gene[-1], count, strand) + "\n")
                    # fo_phasiRNA.write(phasiRNA+ '\n')

def main():
    phase_length = 21
    phasiRNA_file = ''
    allsiRNA_file = ''
    cdna_file = ''
    gdna_file = ''
    flnc_file = ''
    outfile = 'alignment.txt'
    PHAS_fa_outfile = 'PHAS.fa'
    phasiRNA_outfile = 'phasiRNA.fa'
    max_abun = 10
    plot_c = 'y'
    plot_g = 'y'
    plot_f = 'y'
    help = '''
    phase usage:
        option:
            # necessary options:
            -io: file  --  integration -io outputfile
            -ia: file  --  integration -ia outputfile
            -ip: file  --  integration -po outputfile
            -a:  out   --  alignment file
            -o:  out   --  phasiRNA fasta file
            -p:  out   --  PHAS Gene fasta file; Format: >geneid/chr\\tphasiRNA_cluster_region(start end)\\tseq_region(start end)

            # options with default value
            -pl: int   --  phase length, 21 | 24, default=21
            -m:  float  --  the number for reducing the size of Y-axis. default=10

            # optional options
            -c:  file  --  reference transcritome sequence, fasta file, enable cdna based phasiRNA.fa, PHAS.fa, Alignmen, Plot output
            -g:  file  --  reference genome sequence, fasta file, enable gdna based phasiRNA.fa, PHAS.fa, Alignmen, Plot output
            -f:  file  --  full length transcriptome sequence, fasta file, enable flnc based phasiRNA.fa, PHAS.fa, Alignmen, Plot output
            -pc: str  --  plot cdna based phasiRNA cluster, y | n, defaut=y
            -pg: str  --  plot gdna based phasiRNA cluster, y | n, defaut=y
            -pf: str  --  plot flnc based phasiRNA cluster, y | n, defaut=y

            # other
            -v:        --  print version information
            -h:        --  print help information

    '''
    version = '''
    version v1.0
    '''

    for i in range(1, len(sys.argv)):
        if sys.argv[i] == '-io':
            phasiRNA_file = sys.argv[i+1]
        elif sys.argv[i] == '-ia':
            allsiRNA_file = sys.argv[i+1]
        elif sys.argv[i] == '-pl':
            phase_length = int(sys.argv[i+1])
        elif sys.argv[i] == '-c':
            cdna_file = sys.argv[i+1]
        elif sys.argv[i] == '-g':
            gdna_file = sys.argv[i+1]
        elif sys.argv[i] == '-f':
            flnc_file = sys.argv[i+1]
        elif sys.argv[i] == '-a':
            outfile = sys.argv[i+1]
        elif sys.argv[i] == '-pc':
            plot_c = sys.argv[i+1]
        elif sys.argv[i] == '-pg':
            plot_g = sys.argv[i+1]
        elif sys.argv[i] == '-ip':
            PHAS_Loci_file = sys.argv[i+1]
        elif sys.argv[i] == '-pf':
            plot_f = sys.argv[i+1]
        elif sys.argv[i] == '-p':
            PHAS_fa_outfile = sys.argv[i+1]
        elif sys.argv[i] == '-m':
            max_abun = float(sys.argv[i+1])
        elif sys.argv[i] == '-o':
            phasiRNA_outfile = sys.argv[i+1]
        elif sys.argv[i] == '-v':
            print(version)
            sys.exit()
        elif sys.argv[i] == '-h':
            print(help)
            sys.exit()

    fo = open(outfile, 'w')
    fo_PHAS_gene = open(PHAS_fa_outfile, 'w')
    fo_phasiRNA = open(phasiRNA_outfile, 'w')
    cdan_ref = LoadRef(cdna_file)
    flnc_ref = LoadRef(flnc_file)
    gdna_ref = LoadRef(gdna_file)
    tmp_filter_list = Load_PHAS_info(PHAS_Loci_file)
    filter_list = tmp_filter_list[0]
    recorder_dic = tmp_filter_list[1]
    feature_gene_dic = tmp_filter_list[2]
    phasiRNA_dic = LoadsRNA(phasiRNA_file, filter_list)
    allsiRNA_dic = LoadallsRNA(allsiRNA_file, filter_list)
    tmp = Alignment(cdan_ref, gdna_ref, flnc_ref, phasiRNA_dic, fo, phase_length, recorder_dic)
    FaOut(tmp, fo_PHAS_gene, fo_phasiRNA, recorder_dic)
    # transcript_PHAS_Gene = tmp[0]
    # transcript_phasiRNA = tmp[1]
    # genome_PHAS_Gene = tmp[2]
    # genome_phasiRNA = tmp[3]
    # flnc_PHAS_Gene = tmp[4]
    # flnc_phasiRNA = tmp[5]
    ctime = GetCurTime()
    TMP_WD = os.getcwd()
    outdir = TMP_WD + '/' + ctime + '_Fig'
    if plot_c == 'n' and plot_f == 'n' and plot_g == 'n':
        pass
    else:
        os.system(f'mkdir {outdir}')
    # Plot Fig
    Plot(phasiRNA_dic, allsiRNA_dic, phase_length, outdir, max_abun, cdan_ref, gdna_ref, flnc_ref, plot_c, plot_g, plot_f, recorder_dic, feature_gene_dic)

if __name__ == '__main__':
    main()