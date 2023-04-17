# %%
from Bio import SeqIO
import matplotlib
from matplotlib import pyplot as plt
matplotlib.use('Agg')
from collections import defaultdict
import sys
import os
import time

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

def LoadsRNA(sRNA_file):
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
                dic[method][ref][(geneid, pvalue, anno, reocrd_marker)].append((strand, int(pos), abun, seq))
            if method == '>PhaseScore':
                geneid, strand, pos, abun, sRNAid, seq, seq_len, phaseRatio, phaseNumber, phaseAbun, phaseScore, pvalue, anno, reocrd_marker \
                    = l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10], l[11], l[12], l[13]
                dic[method][ref][(geneid, phaseScore, anno, reocrd_marker)].append((strand, int(pos), abun, seq))
    return dic

def LoadallsRNA(sRNA_file):
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
                dic[method][ref][(geneid, pvalue, anno, reocrd_marker)][pos].append((strand, int(pos), abun, seq))
            if method == '>PhaseScore':
                geneid, strand, pos, abun, sRNAid, seq, seq_len, phaseRatio, phaseNumber, phaseAbun, phaseScore, pvalue, anno, reocrd_marker \
                    = l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10], l[11], l[12], l[13]
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

def Alignment(cdna_ref_dic, gdna_ref_dic, flnc_ref_dic,phasiRNA_dic, fo, phase_length):
    transcript_PHAS_Gene = nestedDic()
    transcript_phasiRNA = TwoDepDic() 
    genome_PHAS_Gene = nestedDic()
    genome_phasiRNA = TwoDepDic()
    flnc_PHAS_Gene = nestedDic()
    flnc_phasiRNA = TwoDepDic()
    genome_PHAS_Gene
    for i in phasiRNA_dic:
        fo.write(f'{i}\n')
        for j in phasiRNA_dic[i]:
            fo.write(f'{j}\n')
            if 'cDNA' in j and cdna_ref_dic != {}:
                for k in phasiRNA_dic[i][j]:
                    fo.write('#'+"\t".join(k)+'\n')    
                    start_end = []
                    for e in phasiRNA_dic[i][j][k]:
                        start_end.append(e[1])
                    region = sorted(start_end)
                    region1 = (region[0], region[-1] + phase_length)
                    transcript_PHAS_Gene[k] = PHASGene_seq(cdna_ref_dic[k[0]], region1)
                    for e in phasiRNA_dic[i][j][k]:
                        if e[0] == "+":
                            fo.write(f"{'.'*(int(e[1])-int(region[0]))}{e[3]}{'.'*(int(region[-1])-int(e[1]))}" + "______" + k[0] + "\t" + str(e[1])+ "\n")
                            transcript_phasiRNA[k]['+'].append(e[3])
                    seq = cdna_ref_dic[k[0]][int(region[0]):int(region[-1])+phase_length]
                    fo.write(seq + '______excised_seq' "\n")
                    fo.write(CompleteStrand(seq) + '______anti_seq' + "\n")
                    for e in phasiRNA_dic[i][j][k]:
                        if e[0] == '-':
                            fo.write(f"{'.'*(int(e[1])-int(region[0]))}{CompleteStrand(e[3])}{'.'*(int(region[-1])-int(e[1]))}" + "______" + k[0] + "\t" + str(e[1])+ "\n")
                            transcript_phasiRNA[k]['-'].append(e[3])
            if 'gDNA' in j and gdna_ref_dic != {}:
                for k in phasiRNA_dic[i][j]:
                    fo.write('#'+"\t".join(k)+'\n')    
                    start_end = []
                    for e in phasiRNA_dic[i][j][k]:
                        start_end.append(e[1])
                    region = sorted(start_end)
                    region1 = (region[0], region[-1] + phase_length)
                    genome_PHAS_Gene[k] = PHASGene_seq(gdna_ref_dic[k[0]], region1)
                    for e in phasiRNA_dic[i][j][k]:
                        if e[0] == "+":
                            fo.write(f"{'.'*(int(e[1])-int(region[0]))}{e[3]}{'.'*(int(region[-1])-int(e[1]))}" + "______" + k[0] + "\t" + str(e[1])+ "\n")
                            genome_phasiRNA[k]['+'].append(e[3])
                    seq = gdna_ref_dic[k[0]][int(region[0]):int(region[-1])+phase_length]
                    fo.write(seq + '______excised_seq' "\n")
                    fo.write(CompleteStrand(seq) + '______anti_seq' + "\n")
                    for e in phasiRNA_dic[i][j][k]:
                        if e[0] == '-':
                            fo.write(f"{'.'*(int(e[1])-int(region[0]))}{CompleteStrand(e[3])}{'.'*(int(region[-1])-int(e[1]))}" + "______" + k[0] + "\t" + str(e[1])+ "\n")
                            genome_phasiRNA[k]['-'].append(e[3])
            if 'FLNC' in j and flnc_ref_dic != {}:
                for k in phasiRNA_dic[i][j]:
                    fo.write('#'+"\t".join(k)+'\n')    
                    start_end = []
                    for e in phasiRNA_dic[i][j][k]:
                        start_end.append(e[1])
                    region = sorted(start_end)
                    region1 = (region[0], region[-1] + phase_length)
                    flnc_PHAS_Gene[k] = PHASGene_seq(flnc_ref_dic[k[0]], region1)
                    for e in phasiRNA_dic[i][j][k]:
                        if e[0] == "+":
                            fo.write(f"{'.'*(int(e[1])-int(region[0]))}{e[3]}{'.'*(int(region[-1])-int(e[1]))}" + "______" + k[0] + "\t" + str(e[1])+ "\n")
                            flnc_phasiRNA[k]['+'].append(e[3])
                    seq = flnc_ref_dic[k[0]][int(region[0]):int(region[-1])+phase_length]
                    fo.write(seq + '______excised_seq' "\n")
                    fo.write(CompleteStrand(seq) + '______anti_seq' + "\n")
                    for e in phasiRNA_dic[i][j][k]:
                        if e[0] == '-':
                            fo.write(f"{'.'*(int(e[1])-int(region[0]))}{CompleteStrand(e[3])}{'.'*(int(region[-1])-int(e[1]))}" + "______" + k[0] + "\t" + str(e[1])+ "\n")
                            flnc_phasiRNA[k]['-'].append(e[3])
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
    

def Plot(phasiRNA_dic, allsiRNA_dic, phase_length, outdir, max_abun, cdna_ref, gdna_ref, flnc_ref, plot_c, plot_g, plot_f):
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
                plt.figure(figsize=(10,10))
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
                    sum_abun = 0
                    for pos__ in allsiRNA_dic[method][ref][feature][pos_]:
                        strand = pos__[0]
                        pos = pos__[1]
                        abun = float(pos__[2])
                        sum_abun += abun
                    if strand == '-':
                        if antiStrand_max_abun < sum_abun: sum_abun = antiStrand_max_abun
                        if pos in x1:
                            x5.append(pos)
                            y5.append(0 - sum_abun)
                        else:
                            x2.append(pos)
                            y2.append(0 - sum_abun)
                    if strand == '+':
                        if exsicedStrand_max_abun < sum_abun: sum_abun = exsicedStrand_max_abun
                        if pos in x:
                            x4.append(pos)
                            y4.append(sum_abun)
                        else:
                            x3.append(pos)
                            y3.append(sum_abun)
                x7 = x2 + x3 + x4 + x5
                try:
                    max_pos = sorted(x7)[-1] + phase_length
                except IndexError:
                    count += 1
                    # print(str(count))
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
                if 'PG' in feature[-1] or 'HG' in feature[-1]:
                    # tmp_feature = feature[2].split(';')
                    # tmp_list = []
                    # for fuck in tmp_feature:
                    #     if 'exon' in fuck or 'CDS' in fuck:
                    #         pass
                    #     else:
                    #         tmp_list.append(fuck)
                    if 'Intergenic' in feature[2]:
                        out_name = feature[0]+'_'+feature[2]+'_'+feature[-1]
                    else:
                        out_name = feature[0]+'_'+'Transcript'+'_'+feature[-1]
                else:
                    out_name = feature[0]+'_'+feature[2]+'_'+feature[-1]
                plt.xlabel ("Nucleotide position") 
                plt.ylabel ("Read Abundance \nReverse strand(-)/Forward strand(+)")
                plt.legend (bbox_to_anchor=(0,1.01),loc='lower left', borderaxespad=0)
                plt.title(feature[2], loc='right')
                plt.savefig("{}/{}_plot.jpg".format(outdir, out_name))
                plt.close()

def FaOut(PHAS_Gene_phasiRNA_dic, fo_PHAS_Gene, fo_phasiRNA):
    print('writing PHAS Gene')
    for gene in PHAS_Gene_phasiRNA_dic[0]:
        seq = PHAS_Gene_phasiRNA_dic[0][gene][0]
        start = PHAS_Gene_phasiRNA_dic[0][gene][1][0]
        end = PHAS_Gene_phasiRNA_dic[0][gene][1][1]
        start1 = PHAS_Gene_phasiRNA_dic[0][gene][2][0]
        end1 = PHAS_Gene_phasiRNA_dic[0][gene][2][1]
        fo_PHAS_Gene.write('>{} {}\t{} {}\t{} {}\t{}'.format(gene[0], gene[2], start, end, start1, end1, gene[-1]) + "\n")
        fo_PHAS_Gene.write(seq + '\n')
    for gene in PHAS_Gene_phasiRNA_dic[2]:
        seq = PHAS_Gene_phasiRNA_dic[2][gene][0]
        start = PHAS_Gene_phasiRNA_dic[2][gene][1][0]
        end = PHAS_Gene_phasiRNA_dic[2][gene][1][1]
        start1 = PHAS_Gene_phasiRNA_dic[2][gene][2][0]
        end1 = PHAS_Gene_phasiRNA_dic[2][gene][2][1]
        fo_PHAS_Gene.write('>{} {}\t{} {}\t{} {}\t{}'.format(gene[0], gene[2], start, end, start1, end1, gene[-1])+ '\n')
        fo_PHAS_Gene.write(seq+ '\n')
    for gene in PHAS_Gene_phasiRNA_dic[4]:
        seq = PHAS_Gene_phasiRNA_dic[4][gene][0]
        start = PHAS_Gene_phasiRNA_dic[4][gene][1][0]
        end = PHAS_Gene_phasiRNA_dic[4][gene][1][1]
        start1 = PHAS_Gene_phasiRNA_dic[4][gene][2][0]
        end1 = PHAS_Gene_phasiRNA_dic[4][gene][2][1]
        fo_PHAS_Gene.write('>{} {}\t{} {}\t{} {}\t{}'.format(gene[0], gene[2], start, end, start1, end1, gene[-1])+ '\n')
        fo_PHAS_Gene.write(seq+ '\n')
    print('writing phasiRNA')
    for gene in PHAS_Gene_phasiRNA_dic[1]:
        count = 0
        for strand in PHAS_Gene_phasiRNA_dic[1][gene]:
            for phasiRNA in PHAS_Gene_phasiRNA_dic[1][gene][strand]:
                count += 1
                fo_phasiRNA.write('>{} {} {} phasiRNA {} {}'.format(gene[0], gene[2], gene[-1], count, strand) + "\n")
                fo_phasiRNA.write(phasiRNA+ '\n')
    for gene in PHAS_Gene_phasiRNA_dic[3]:
        count = 0
        for strand in PHAS_Gene_phasiRNA_dic[3][gene]:
            for phasiRNA in PHAS_Gene_phasiRNA_dic[3][gene][strand]:
                count += 1
                fo_phasiRNA.write('>{} {} {} phasiRNA {} {}'.format(gene[0], gene[2], gene[-1], count, strand) + "\n")
                fo_phasiRNA.write(phasiRNA+ '\n')
    for gene in PHAS_Gene_phasiRNA_dic[5]:
        count = 0
        for strand in PHAS_Gene_phasiRNA_dic[5][gene]:
            for phasiRNA in PHAS_Gene_phasiRNA_dic[5][gene][strand]:
                count += 1
                fo_phasiRNA.write('>{} {} {} phasiRNA {} {}'.format(gene[0], gene[2], gene[-1], count, strand) + "\n")
                fo_phasiRNA.write(phasiRNA+ '\n')

def main():
    phase_length = 21
    phasiRNA_file = ''
    allsiRNA_file = ''
    cdna_file = ''
    gdna_file = ''
    flnc_file = ''
    outfile = ''
    PHAS_fa_outfile = ''
    phasiRNA_outfile = ''
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
            -pc: str  --  plot cdna based phasiRNA cluster
            -pg: str  --  plot gdna based phasiRNA cluster
            -pf: str  --  plot flnc based phasiRNA cluster

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
    phasiRNA_dic = LoadsRNA(phasiRNA_file)
    allsiRNA_dic = LoadallsRNA(allsiRNA_file)
    tmp = Alignment(cdan_ref, gdna_ref, flnc_ref, phasiRNA_dic, fo, phase_length)
    FaOut(tmp, fo_PHAS_gene, fo_phasiRNA)
    # transcript_PHAS_Gene = tmp[0]
    # transcript_phasiRNA = tmp[1]
    # genome_PHAS_Gene = tmp[2]
    # genome_phasiRNA = tmp[3]
    # flnc_PHAS_Gene = tmp[4]
    # flnc_phasiRNA = tmp[5]
    ctime = GetCurTime()
    TMP_WD = os.getcwd()
    outdir = TMP_WD + '/' + ctime + '_Fig'
    os.system(f'mkdir {outdir}')
    # Plot Fig
    Plot(phasiRNA_dic, allsiRNA_dic, phase_length, outdir, max_abun, cdan_ref, gdna_ref, flnc_ref, plot_c, plot_g, plot_f)

if __name__ == '__main__':
    main()