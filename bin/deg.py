# %%
# -*- encoding: utf-8 -*-
import sys
from Bio import SeqIO
import time
import os
from matplotlib import pyplot as plt
from customeDic import *
'''
@File    :   deg.py
@Time    :   2022/11/05 23:50:16
@Author  :   zrf 
@Version :   0.0
@Contact :   zrfeng1@gmail.com
@License :   (C)Copyright 2022-, NJAU-CBI
@Desc    :   None
'''
# --------------------------------------------------------------> DEFAULT VARIABLE <---------------------------------------------------------------- #
start_time = time.time()
shifts = [-1, 0, 1]
min_num = 0
plot_function = 'n'
library = 'n'
plot12=1
less = 'n'
outfoldername='./'
initiator = 'n'
help = '''
// designed by Ji Huang, CBG, NJAU; revised by Yuhan Fei, Baoyi Zhang and Zerong Feng; 
    http://cbi.njau.edu.cn
// function: vertified the sRNA - Target interaction with degradome data

    options:
    -i: <inputfilename>     --    mapping file for degradome data mapping transcripts, by bowtie
    -q: <sRNA fasta>        --    small RNA sequences used for target prediction, fasta
    -j: <inputfilename>     --    from psRNATarget batch download file or initiator output
    -t: <inputfilename>     --    transcripts file, fasta
    -o: <outputfilename>    --    matched map file with only matched records
    -s: <shift_number>      --    if shifts=0 then cleaved exactly at pos.10, default=1
    -m: <minum deg_num>     --    minum number of degradome reads, int, default=0
    -p: <T-plot function>   --    enable the plot function, y | n, default='n' 
    -in: <bool>             --    y | n, use initiator output information
    -pl [int]               --    1,plot only category 1; 2, plot categories 1 and 2, default=1"
    -pf [str]               --    output folder name, for exporting t-plot images, optional, useful when enable plot function"
    --lib [str]             --    library, optional
    -less                   --    only output cat_1 and cat_2 information

    ***********************
    //About the categories:
    Cat #1, degradome read at the cleavage site is most abundant.
    Cat #2, the read is less than the most abudant one, but higher than the median.
    Cat #3, the read is less than the median, but high than 1
    Cat #4, the read is identical or less than 1 (if degradome data is normalized)
    "
'''
version = '''
    version 0.1
'''
filter = ''
# --------------------------------------------------------------> END <---------------------------------------------------------------- #
for i in range(1, len(sys.argv)):
    if sys.argv[i] == '-i':
        mapfile = sys.argv[i+1]
    elif sys.argv[i] == '-q':
        inp = sys.argv[i+1]
    elif sys.argv[i] == '-j':
        psRNAfile = sys.argv[i+1]
    elif sys.argv[i] == '-t':
        targetfa = sys.argv[i+1]
    elif sys.argv[i] == '-s':
        SHIFT = int(sys.argv[i+1])
    elif sys.argv[i] == '--lib':
        library = 'y'
        libraryinfo = sys.argv[i+1]
    elif sys.argv[i]=='-pl':
        plot12= int(sys.argv[i+1])
    elif sys.argv[i]=='-pf':
        outfoldername= sys.argv[i+1]
    elif sys.argv[i] == '-m':
        min_num = float(sys.argv[i+1])
    elif sys.argv[i] == '-in':
        initiator = sys.argv[i+1]
    elif sys.argv[i] == '-o':
        outputfile = sys.argv[i+1]
    elif sys.argv[i] == '-less':
        less = 'y'
    elif sys.argv[i] == '-p':
        plot_function = sys.argv[i+1]
    elif sys.argv[i] == '-f':
        filter = sys.argv[i+1]
    elif sys.argv[i] == '-v':
        print(version)
        sys.exit()
    elif sys.argv[i] == '-h':
        print(help)
        sys.exit()

if 'SHIFT' in dir():
    shifts = []
    for i in range(SHIFT + 1):
        shifts.append(-i)
        shifts.append(i)

#######################################--> test module <--###############################################
# test block, you can delete me anytime and anywhere or you can leave me for debug next time
#  inp = '/home/riceUsers/fzr/new_pipline_test/result/cdna_m10_1e-05_phasiRNA.fa'
#  mapfile = '/home/riceUsers/fzr/new_pipline_test/degradome_data/outdir/GSM1040650/GSM1040650_norm_filter.map'
#  psRNAfile = '/home/riceUsers/fzr/new_pipline_test/result/psRNATargetJob-1664868827125964.txt'
#  targetfa = '/home/riceUsers/fzr/new_pipline_test/result/osa_msu7_cdna.fa'
#  outputfile = '/home/riceUsers/fzr/deg_out'
#########################################################################################################
print('loading miRNA seq data ... ')
miRNA_dic = {}
with open(inp, 'r') as fn:
    for i in SeqIO.parse(fn, 'fasta'):
        if i.name not in miRNA_dic:
            miRNA_dic[i.name] = {str(i.seq): len(i.seq)} 

print('loading Prediction Target Gene fasta data ... ')
targetfaDic = {}
with open(targetfa, 'r') as fn:
    for i in SeqIO.parse(fn, 'fasta'):
        if i.name not in targetfaDic:
            targetfaDic[i.name] = len(i.seq)

# with open(targetfa, 'r') as fn:
#     last_line = "CBI"
#     for line_ in fn:
#         line = line_.strip()
#         if last_line.startswith('>') and ~line.startswith('>'):
#             if last_line.replace('>','') not in targetfaDic:
#                 targetfaDic[last_line.replace('>', '')] = len(line)
#         last_line = line
if filter != '':
    filter_dic = OneDepDic()   
    with open(filter, 'r') as fn:
        for line in fn:
            line1 = line.strip()
            l = line.strip().split("\t")
            miRNA = l[0]
            target = l[1]
            filter_dic[target].append(miRNA)

print('loading psRNATarget result data ... ')
psRNA_dic = {}
with open(psRNAfile, 'r') as fn:
    if initiator == 'n':
        for line in fn:
            line1 = line.strip()
            l = line.strip().split("\t")
            if line1.startswith('#'):
                continue
            if l[0] == 'miRNA_Acc.':
                continue
            target_gene = l[1]
            miRNA = l[0]
            miRNA_loc = int(l[7])
            try:
                gene_annotation = l[11]
            except IndexError:
                gene_annotation = 'NULL'
            if target_gene not in psRNA_dic:
                psRNA_dic[target_gene] = {miRNA: {miRNA_loc: gene_annotation}}
            else:
                pass
            if miRNA not in psRNA_dic[target_gene]:
                psRNA_dic[target_gene][miRNA] = {miRNA_loc: gene_annotation}
            else:
                pass
            if miRNA_loc not in psRNA_dic[target_gene][miRNA]:
                psRNA_dic[target_gene][miRNA][miRNA_loc] = gene_annotation
            else:
                pass
    elif initiator == 'y':
        for line in fn:
            line1 = line.strip()
            l = line.strip().split("\t")
            if line1.startswith('miRNA'):
                continue
            target_gene = l[1]
            miRNA = l[0]
            miRNA_loc = int(l[3]) + 1
            try:
                gene_annotation = l[11]
            except IndexError:
                gene_annotation = l[-2]
            if target_gene not in psRNA_dic:
                psRNA_dic[target_gene] = {miRNA: {miRNA_loc: gene_annotation}}
            else:
                pass
            if miRNA not in psRNA_dic[target_gene]:
                psRNA_dic[target_gene][miRNA] = {miRNA_loc: gene_annotation}
            else:
                pass
            if miRNA_loc not in psRNA_dic[target_gene][miRNA]:
                psRNA_dic[target_gene][miRNA][miRNA_loc] = gene_annotation
            else:
                pass


print('loading degra data and start analysis ... ')
# 保存在deg_target 中的
int_ = {}
annotation = {}
with open(mapfile, 'r') as fn:
    for line in fn:
        line1 = line.strip()
        l = line.strip().split("\t")
        deg_num = float((l[0].split('@')[1]))
        deg_loc = int(l[3]) + 1
        target_gene = l[2]
        if target_gene in psRNA_dic:
            miRNAs = psRNA_dic[target_gene]
            for miRNA in miRNAs:
                miRNA_locs = psRNA_dic[target_gene][miRNA]
                for miRNA_loc in miRNA_locs:
                    gene_annotation = psRNA_dic[target_gene][miRNA][miRNA_loc]
                    if target_gene not in int_:
                        int_[target_gene] = {targetfaDic[target_gene]: None}
                        annotation[target_gene] = gene_annotation
                    else:
                        pass
                    if 'sRNA' not in int_[target_gene]:
                        int_[target_gene]['sRNA'] = {}
                    if miRNA not in int_[target_gene]['sRNA']:
                            int_[target_gene]['sRNA'][miRNA] = {miRNA_loc: None}
                    if miRNA_loc not in int_[target_gene]['sRNA'][miRNA]:
                            int_[target_gene]['sRNA'][miRNA][miRNA_loc] = None
            if 'deg' not in int_[target_gene]:
                int_[target_gene]['deg'] = {}
            if deg_loc not in int_[target_gene]['deg']:
                int_[target_gene]['deg'][deg_loc] = deg_num
            else:
                int_[target_gene]['deg'][deg_loc] += deg_num

def get_median(data):
    data.sort()
    half = len(data)//2
    return (data[half] + data[~half])/2

for gene in list(int_):
    max_deg_num = max(int_[gene]['deg'].values())
    median_deg_num = get_median(list(int_[gene]['deg'].values()))
    if 'max_deg_num' not in int_:
        int_[gene]['max_deg_num'] = max_deg_num
    if 'median_deg_num' not in int_:
        int_[gene]['median_deg_num'] = median_deg_num
    for sRNA in int_[gene]['sRNA']:
        for sRNA_loc in list(int_[gene]['sRNA'][sRNA].keys()):
            for shift in shifts:
                if sRNA_loc -10 + 1 - shift in int_[gene]['deg']:
                    cut_site_num = int_[gene]['deg'][sRNA_loc - 10 + 1 - shift]
                    if cut_site_num >= int_[gene]['max_deg_num']:
                        if 'cat' not in int_[gene]:
                            int_[gene]['cat'] = {}
                        if 'Cat_1' not in int_[gene]['cat']:
                            int_[gene]['cat']['Cat_1'] = {}
                        if sRNA not in int_[gene]['cat']['Cat_1']:
                            int_[gene]['cat']['Cat_1'][sRNA] = {}
                        if sRNA_loc not in int_[gene]['cat']['Cat_1'][sRNA]:
                            int_[gene]['cat']['Cat_1'][sRNA][sRNA_loc] = {}
                        if shift not in int_[gene]['cat']['Cat_1'][sRNA][sRNA_loc]:
                            int_[gene]['cat']['Cat_1'][sRNA][sRNA_loc][shift] = cut_site_num
                    elif cut_site_num < int_[gene]['max_deg_num'] and cut_site_num >= int_[gene]['median_deg_num']:
                        if 'cat' not in int_[gene]:
                            int_[gene]['cat'] = {}
                        if 'Cat_2' not in int_[gene]['cat']:
                            int_[gene]['cat']['Cat_2'] = {}
                        if sRNA not in int_[gene]['cat']['Cat_2']:
                            int_[gene]['cat']['Cat_2'][sRNA] = {}
                        if sRNA_loc not in int_[gene]['cat']['Cat_2'][sRNA]:
                            int_[gene]['cat']['Cat_2'][sRNA][sRNA_loc] = {}
                        if shift not in int_[gene]['cat']['Cat_2'][sRNA][sRNA_loc]:
                            int_[gene]['cat']['Cat_2'][sRNA][sRNA_loc][shift] = cut_site_num
                    elif cut_site_num < int_[gene]['median_deg_num'] and cut_site_num > 1:
                        if 'cat' not in int_[gene]:
                            int_[gene]['cat'] = {}
                        if 'Cat_3' not in int_[gene]['cat']:
                            int_[gene]['cat']['Cat_3'] = {}
                        if sRNA not in int_[gene]['cat']['Cat_3']:
                            int_[gene]['cat']['Cat_3'][sRNA] = {}
                        if sRNA_loc not in int_[gene]['cat']['Cat_3'][sRNA]:
                            int_[gene]['cat']['Cat_3'][sRNA][sRNA_loc] = {}
                        if shift not in int_[gene]['cat']['Cat_3'][sRNA][sRNA_loc]:
                            int_[gene]['cat']['Cat_3'][sRNA][sRNA_loc][shift] = cut_site_num
                    elif cut_site_num < int_[gene]['median_deg_num'] and cut_site_num <= 1:
                        if 'cat' not in int_[gene]:
                            int_[gene]['cat'] = {}
                        if 'Cat_4' not in int_[gene]['cat']:
                            int_[gene]['cat']['Cat_4'] = {}
                        if sRNA not in int_[gene]['cat']['Cat_4']:
                            int_[gene]['cat']['Cat_4'][sRNA] = {}
                        if sRNA_loc not in int_[gene]['cat']['Cat_4'][sRNA]:
                            int_[gene]['cat']['Cat_4'][sRNA][sRNA_loc] = {}
                        if shift not in int_[gene]['cat']['Cat_4'][sRNA][sRNA_loc]:
                            int_[gene]['cat']['Cat_4'][sRNA][sRNA_loc][shift] = cut_site_num
                    else:
                        if 'cat' not in int_[gene]:
                            int_[gene]['cat'] = {}
                        if 'Cat_NA' not in int_[gene]['cat']:
                            int_[gene]['cat']['Cat_NA'] = {}
                        if sRNA not in int_[gene]['cat']['Cat_NA']:
                            int_[gene]['cat']['Cat_NA'][sRNA] = {}
                        if sRNA_loc not in int_[gene]['cat']['Cat_NA'][sRNA]:
                            int_[gene]['cat']['Cat_NA'][sRNA][sRNA_loc] = {}
                        if shift not in int_[gene]['cat']['Cat_NA'][sRNA][sRNA_loc]:
                            int_[gene]['cat']['Cat_NA'][sRNA][sRNA_loc][shift] = cut_site_num

print(f'analysis finished, start writing data to {outputfile}')

with open(outputfile, 'w') as fo:
    if library == 'y':
        string_list = ['Category', 'Small_RNA', 'Target_gene', 'sRNA_loc', 'Deg_loc', 'Deg_count', 'sRNA_seq', 'Shift', 'Gene_annotation', 'Library']
    elif library == 'n':
        string_list = ['Category', 'Small_RNA', 'Target_gene', 'sRNA_loc', 'Deg_loc', 'Deg_count', 'sRNA_seq', 'Shift', 'Gene_annotation']
    fo.write("\t".join(string_list) + "\n")
    for gene in int_:
        gene_annotation = annotation[gene]
        if 'cat' in int_[gene]:
            for category in int_[gene]['cat']:
                for vertified_sRNA in int_[gene]['cat'][category]:
                    sRNA_seq = list(miRNA_dic[vertified_sRNA].keys())[0]
                    for sRNA_loc in int_[gene]['cat'][category][vertified_sRNA]:
                        for shift in int_[gene]['cat'][category][vertified_sRNA][sRNA_loc]:
                            deg_loc = sRNA_loc - 10 + 1 - shift
                            deg_count = int_[gene]['cat'][category][vertified_sRNA][sRNA_loc][shift]
                            if deg_count <= min_num:
                                continue
                            if less == 'y' and category == 'Cat_3':
                                continue
                            if less == 'y' and category == 'Cat_4':
                                continue
                            if library == 'y':
                                string_list = [category, vertified_sRNA, gene, str(sRNA_loc), str(deg_loc), str(deg_count), str(sRNA_seq), str(shift), gene_annotation, libraryinfo]
                            elif library == 'n':
                                string_list = [category, vertified_sRNA, gene, str(sRNA_loc), str(deg_loc), str(deg_count), str(sRNA_seq), str(shift), gene_annotation]
                            fo.write("\t".join(string_list) + "\n")

end_time = time.time()
print("total time:", end_time - start_time)

if plot_function == 'y':
    if outfoldername!='':
        # if not os.path.exists(outfoldername):
        #     os.makedirs(outfoldername)
        pass
    if outfoldername!='':
        if plot12==1:
            print("NOTE: only ploting category 1 images!")
        elif plot12==2:
            print("NOTE: only ploting categories 1 and 2 images!")
        print("Generating t-plot images...")
        
    if filter != '':
        for gene in filter_dic:
            try:
                if 'cat' in int_[gene]:
                    for category in int_[gene]['cat']:
                        cat_number = int(category.split('_')[1])
                        y_value = []
                        x_coor = []
                        if cat_number > plot12:
                            continue
                        else:
                            max_deg_num = int_[gene]['max_deg_num']
                            median_deg_num = int_[gene]['median_deg_num']
                            for deg_loc in int_[gene]['deg']:
                                if cat_number == 1 and int_[gene]['deg'][deg_loc] != max_deg_num:
                                    x_coor.append(deg_loc)
                                    y_value.append(int_[gene]['deg'][deg_loc])
                                elif cat_number == 2 and int_[gene]['deg'][deg_loc] != median_deg_num: 
                                    x_coor.append(deg_loc)
                                    y_value.append(int_[gene]['deg'][deg_loc])
                            for vertified_sRNA in int_[gene]['cat'][category]:
                                for sRNA_loc in int_[gene]['cat'][category][vertified_sRNA]: 
                                    for shift in int_[gene]['cat'][category][vertified_sRNA][sRNA_loc]:
                                        deg_loc = sRNA_loc - 10 + 1 - shift
                                        deg_count = int_[gene]['cat'][category][vertified_sRNA][sRNA_loc][shift]
                                        red_y_value = [deg_count]
                                        red_x_coor = [deg_loc]
                                if deg_count <= min_num:
                                    continue
                                plt.figure(figsize=(8,6))
                                plt.grid(True)
                                plt.title(vertified_sRNA+'/'+gene)
                                plt.plot(x_coor,y_value,'o',color='black')
                                plt.plot(red_x_coor,red_y_value,'o',color='red')
                                plt.xlabel ("Nucleotide position (nt)")
                                plt.ylabel ("Degradome Abundance (RP10M)")
                                plt.ylim(ymin=0)
                                outfolder_sufix = outfoldername.split('/')[1]
                                if outfoldername == './':
                                    plt.savefig (outfoldername+'/'+outfolder_sufix + vertified_sRNA+'-'+gene+'.png')
                                else:
                                    plt.savefig ('/'.join(outfoldername.split('/')[:-1])+'/'+outfolder_sufix + '_' + vertified_sRNA+'-'+gene+'.png')
                                plt.close()
            except KeyError:
                continue
    else:
        for gene in int_:
            if 'cat' in int_[gene]:
                for category in int_[gene]['cat']:
                    cat_number = int(category.split('_')[1])
                    y_value = []
                    x_coor = []
                    if cat_number > plot12:
                        continue
                    else:
                        max_deg_num = int_[gene]['max_deg_num']
                        median_deg_num = int_[gene]['median_deg_num']
                        for deg_loc in int_[gene]['deg']:
                            if cat_number == 1 and int_[gene]['deg'][deg_loc] != max_deg_num:
                                x_coor.append(deg_loc)
                                y_value.append(int_[gene]['deg'][deg_loc])
                            elif cat_number == 2 and int_[gene]['deg'][deg_loc] != median_deg_num: 
                                x_coor.append(deg_loc)
                                y_value.append(int_[gene]['deg'][deg_loc])
                        for vertified_sRNA in int_[gene]['cat'][category]:
                            for sRNA_loc in int_[gene]['cat'][category][vertified_sRNA]: 
                                for shift in int_[gene]['cat'][category][vertified_sRNA][sRNA_loc]:
                                    deg_loc = sRNA_loc - 10 + 1 - shift
                                    deg_count = int_[gene]['cat'][category][vertified_sRNA][sRNA_loc][shift]
                                    red_y_value = [deg_count]
                                    red_x_coor = [deg_loc]
                            if deg_count <= min_num:
                                continue
                            plt.figure(figsize=(5.5,4), dpi=300)
                            plt.grid(True)
                            plt.title(vertified_sRNA+'/'+gene)
                            plt.plot(x_coor,y_value,'o',color='black')
                            plt.plot(red_x_coor,red_y_value,'o',color='red')
                            plt.xlabel ("Nucleotide position (nt)")
                            plt.ylabel ("Degradome Abundance (RP10M)")
                            plt.ylim(ymin=0)
                            outfolder_sufix = outfoldername.split('/')[1]
                            if outfoldername == './':
                                plt.savefig (outfoldername+'/'+outfolder_sufix + vertified_sRNA+'-'+gene+'.png')
                            else:
                                plt.savefig ('/'.join(outfoldername.split('/')[:-1])+'/'+outfolder_sufix + '_' + vertified_sRNA+'-'+gene+'.png')
                            plt.close()



# --------------------------------------------------------------> shift=-1 <---------------------------------------------------------------- #
# deg_reads
# TACGCGATATA
#           ATTTGAGAGTACGCGATATAT
#           sRNA_reads
# ==============================================================> shift=1 <================================================================ #
# deg_reads
# TACGCGATATA
#         ATATTTGAGAGTACGCGATATAT
#         sRNA_reads
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~> shift=0 <~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# deg_reads
# TACGCGATATA
#          TATTTGAGAGTACGCGATATAT
#          sRNA_reads
# sRNA_loc -10 + 1 - shift = deg_loc
