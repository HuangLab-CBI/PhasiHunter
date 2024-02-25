import time
from sympy import binomial
from collections import defaultdict
import re
import os
from concurrent.futures import ProcessPoolExecutor, wait, ALL_COMPLETED

def GetCurTime():
    """return formated current localtime 

    Returns:
        <str> -- formated current localtime
    """
    ctime = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
    return ctime

def PrintDic(dic, maxdep, prefix='\t'):
    # WARNING: not suport print list
    count = 0
    for i in dic.keys():
        count += 1
        if count > maxdep:
            break
        else:
            if isinstance(dic[i], dict):
                print(prefix + '+ ' + str(i))
                PrintDic(dic[i], maxdep, prefix + '|   ')  # 最后一个
            else:
                print(prefix + '- ' + str(i) + ' ==> ' + str(dic[i]))

def Vprint(string, enable = False, time = True):
    """verbose print

    Parameters
    ----------
    string : str
    enable : bool, optional
        by default False
    time : bool, optional
        by default True
    """    
    if time:
        if enable:
            print(string + " at " + GetCurTime())
    else:
        if enable:
            print(string)

def checkgff3(gff3file):
    gff3type = 'unknown'
    with open(gff3file, 'r') as fn:
        for i in range(10000):
            line = fn.readline()
            if re.search('NCBI_Assembly', line, flags=re.IGNORECASE):
                gff3type = 'ncbi'
                print(f'Detected gff3 annotation: NCBI')
                break
            elif re.search('phytozome', line, flags=re.I):
                gff3type = 'phytozome'
                print(f'Detected gff3 annotation: phytozome')
                break
            else:
                gff3type = 'ensemble'
    
    return gff3type

def Gff3Tobed(gff3file, outfile, gff3type):
    fo = open(outfile, 'w')
    spe = 'NA'
    with open(gff3file, 'r') as fn:
        for line in fn:
            l = line.strip().split("\t")
            if line.startswith("#"):
                try:
                    version = re.search('annot-version.*', line).group().replace('annot-version ', '')
                except AttributeError:
                    pass
            if line.startswith("#"):
                try:
                    spe = re.search('species.*', line).group().replace('species ', '')
                except AttributeError:
                    pass
                continue
            chr_, feature, start, end, longstring, strand = l[0], l[2], l[3], l[4], l[8], l[6]
            if int(start) > int(end): 
                continue
            if feature == 'region':
                continue
            if gff3type == 'ncbi':
                if re.match('NW_', chr_):
                    continue
                anno = "-".join(longstring.split(";")[0].split('-')[1:])
                if anno == '':
                    continue
                else:
                    fo.write(f'{chr_}\t{start}\t{end}\t{feature}\t{anno}\t{strand}\n')
            elif gff3type == 'phytozome':
                if re.match('chrsy|chrun', chr_, re.I):
                    continue
                if spe.lower() == 'oryza sativa':
                    anno = longstring.split(';')[0].split('=')[1].replace('.MSU' + version, '').replace('.'+feature+'.', '-')
                else:
                    anno = longstring.split(';')[0].split('=')[1].replace('.' + version, '').replace('.'+feature+'.', '-')
                fo.write(f'{chr_}\t{start}\t{end}\t{feature}\t{anno}\t{strand}\n')
            elif gff3type == 'ensemble':
                if not re.match('\d+', chr_):
                    continue
                if "ID=chromosome" in longstring:
                    continue
                if feature != 'exon':
                    try:
                        anno = longstring.split(';')[0].split(':')[1]
                    except IndexError:
                        pass
                elif feature == 'exon':
                    anno = longstring.split(';')[1].split('=')[1]
                fo.write(f'{chr_}\t{start}\t{end}\t{feature}\t{anno}\t{strand}\n')
    fo.close()

def nestedDic():
    return defaultdict(nestedDic)

def defauldict_list():
    return defaultdict(list)

def TwoDepDic():
    return defaultdict(defauldict_list)

def ThreeDepDic():
    return defaultdict(TwoDepDic)

def FourDepDic():
    return defaultdict(ThreeDepDic)

def FiveDepDic():
    return defaultdict(FourDepDic)

def LoadData(file):
    final_cluster = TwoDepDic()
    final_strand = nestedDic()
    final_phasiRNA = nestedDic()
    p_final_cluster = ThreeDepDic()
    p_final_strand = nestedDic()
    p_final_phasiRNA = nestedDic()
    hc = ['CF', 'CR']
    hg = ['GF', 'GR']
    hf = ['FF', 'FR']
    pc = ['PC']
    pg = ['PG']
    pf = ["PF"]
    with open(file, 'r') as fn:
        for line in fn:
            l = line.strip().split("\t")
            geneid, strand, pos, abun, sRNAid, seq, seqLen, pvalue, phaseRatio, phaseNumber, phaseAbun, phaseScore, recordMarker = \
                l[0], l[1], int(l[2]), l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10], l[11], l[12]
            category = recordMarker.split('_')[-1][0:2]
            if category in hc:
                final_cluster['hc'][geneid].append(pos)
                final_strand['hc'][geneid][pos] = strand
                key = geneid + '\t' + strand + '\t' + str(pos)
                value = geneid + "\t" + strand + "\t" + str(pos) + "\t" + abun + "\t" + sRNAid + "\t" + seq + "\t" + seqLen + "\t" + phaseRatio + "\t" + phaseNumber + "\t" + phaseAbun + "\t" + phaseScore + "\t" + recordMarker + "\t" + pvalue
                final_phasiRNA['hc'][geneid][key] = value
            elif category in hg:
                final_cluster['hg'][geneid].append(pos)
                final_strand['hg'][geneid][pos] = strand
                key = geneid + '\t' + strand + '\t' + str(pos)
                value = geneid + "\t" + strand + "\t" + str(pos) + "\t" + abun + "\t" + sRNAid + "\t" + seq + "\t" + seqLen + "\t" + phaseRatio + "\t" + phaseNumber + "\t" + phaseAbun + "\t" + phaseScore + "\t" + recordMarker + "\t" + pvalue
                final_phasiRNA['hg'][geneid][key] = value
            elif category in hf:
                final_cluster['hf'][geneid].append(pos)
                final_strand['hf'][geneid][pos] = strand
                key = geneid + '\t' + strand + '\t' + str(pos)
                value = geneid + "\t" + strand + "\t" + str(pos) + "\t" + abun + "\t" + sRNAid + "\t" + seq + "\t" + seqLen + "\t" + phaseRatio + "\t" + phaseNumber + "\t" + phaseAbun + "\t" + phaseScore + "\t" + recordMarker + "\t" + pvalue
                final_phasiRNA['hf'][geneid][key] = value
            elif category in pc:
                p_final_cluster['pc'][geneid][recordMarker].append(pos)
                p_final_strand['pc'][geneid][pos] = strand
                key = geneid + '\t' + strand + '\t' + str(pos)
                value = geneid + "\t" + strand + "\t" + str(pos) + "\t" + abun + "\t" + sRNAid + "\t" + seq + "\t" + seqLen + "\t" + phaseRatio + "\t" + phaseNumber + "\t" + phaseAbun + "\t" + phaseScore + "\t" + recordMarker + "\t" + pvalue
                p_final_phasiRNA['pc'][geneid][recordMarker][key] = value
            elif category in pg:
                p_final_cluster['pg'][geneid][recordMarker].append(pos)
                p_final_strand['pg'][geneid][pos] = strand
                key = geneid + '\t' + strand + '\t' + str(pos)
                value = geneid + "\t" + strand + "\t" + str(pos) + "\t" + abun + "\t" + sRNAid + "\t" + seq + "\t" + seqLen + "\t" + phaseRatio + "\t" + phaseNumber + "\t" + phaseAbun + "\t" + phaseScore + "\t" + recordMarker + "\t" + pvalue
                p_final_phasiRNA['pg'][geneid][recordMarker][key] = value
            elif category in pf:
                p_final_cluster['pf'][geneid][recordMarker].append(pos)
                p_final_strand['pf'][geneid][pos] = strand
                key = geneid + '\t' + strand + '\t' + str(pos)
                value = geneid + "\t" + strand + "\t" + str(pos) + "\t" + abun + "\t" + sRNAid + "\t" + seq + "\t" + seqLen + "\t" + phaseRatio + "\t" + phaseNumber + "\t" + phaseAbun + "\t" + phaseScore + "\t" + recordMarker + "\t" + pvalue
                p_final_phasiRNA['pf'][geneid][recordMarker][key] = value
    return (final_cluster, final_strand, final_phasiRNA, p_final_cluster, p_final_strand, p_final_phasiRNA)

def LoadallsiRNA(file):
    dic = nestedDic()
    odic = nestedDic()
    oodic = nestedDic()
    hc = ['CF', 'CR']
    hg = ['GF', 'GR']
    hf = ['FF', 'FR']
    pc = ['PC']
    pg = ['PG']
    pf = ['PF']
    with open(file, 'r') as fn:
        for line in fn:
            l = line.strip().split("\t")
            geneid, strand, pos, abun, sRNAid, seq, seqLen, pvalue, phaseRatio, phaseNumber, phaseAbun, phaseScore, recordMarker = \
                l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10], l[11], l[12]
            category = recordMarker.split('_')[-1][0:2]
            if category in hc:
                key = geneid + '\t' + strand + '\t' + pos
                value = geneid + "\t" + strand + "\t" + pos + "\t" + abun + "\t" + sRNAid + "\t" + seq + "\t" + seqLen + "\t" + phaseRatio + "\t" + phaseNumber + "\t" + phaseAbun + "\t" + phaseScore + "\t" + recordMarker + "\t" + pvalue
                dic['hc'][geneid][key] = value
                key1 = geneid + '\t' + strand + '\t' + pos + "\t" + seqLen + "\t" + recordMarker
                odic['hc'][geneid][key1] = value
            elif category in hg:
                key = geneid + '\t' + strand + '\t' + pos
                value = geneid + "\t" + strand + "\t" + pos + "\t" + abun + "\t" + sRNAid + "\t" + seq + "\t" + seqLen + "\t" + phaseRatio + "\t" + phaseNumber + "\t" + phaseAbun + "\t" + phaseScore + "\t" + recordMarker + "\t" + pvalue
                dic['hg'][geneid][key] = value
                key1 = geneid + '\t' + strand + '\t' + pos + "\t" + seqLen + "\t" + recordMarker
                odic['hg'][geneid][key1] = value
            elif category in hf:
                key = geneid + '\t' + strand + '\t' + pos
                value = geneid + "\t" + strand + "\t" + pos + "\t" + abun + "\t" + sRNAid + "\t" + seq + "\t" + seqLen + "\t" + phaseRatio + "\t" + phaseNumber + "\t" + phaseAbun + "\t" + phaseScore + "\t" + recordMarker + "\t" + pvalue
                dic['hf'][geneid][key] = value
                key1 = geneid + '\t' + strand + '\t' + pos + "\t" + seqLen + "\t" + recordMarker
                odic['hf'][geneid][key1] = value
            if category in pc:
                key = geneid + '\t' + strand + '\t' + pos
                value = geneid + "\t" + strand + "\t" + pos + "\t" + abun + "\t" + sRNAid + "\t" + seq + "\t" + seqLen + "\t" + phaseRatio + "\t" + phaseNumber + "\t" + phaseAbun + "\t" + phaseScore + "\t" + recordMarker + "\t" + pvalue
                dic['pc'][geneid][key] = value
                key1 = geneid + '\t' + strand + '\t' + pos + "\t" + seqLen + "\t" + recordMarker
                odic['pc'][geneid][key1] = value
                key2 = geneid + '\t' + strand + '\t' + pos + "\t" + seqLen + "\t" + recordMarker
                oodic['pc'][geneid][recordMarker][key2] = value
            if category in pg:
                key = geneid + '\t' + strand + '\t' + pos
                value = geneid + "\t" + strand + "\t" + pos + "\t" + abun + "\t" + sRNAid + "\t" + seq + "\t" + seqLen + "\t" + phaseRatio + "\t" + phaseNumber + "\t" + phaseAbun + "\t" + phaseScore + "\t" + recordMarker + "\t" + pvalue
                dic['pg'][geneid][key] = value
                key1 = geneid + '\t' + strand + '\t' + pos + "\t" + seqLen + "\t" + recordMarker
                odic['pg'][geneid][key1] = value
                key2 = geneid + '\t' + strand + '\t' + pos + "\t" + seqLen + "\t" + recordMarker
                oodic['pg'][geneid][recordMarker][key2] = value
            if category in pf:
                key = geneid + '\t' + strand + '\t' + pos
                value = geneid + "\t" + strand + "\t" + pos + "\t" + abun + "\t" + sRNAid + "\t" + seq + "\t" + seqLen + "\t" + phaseRatio + "\t" + phaseNumber + "\t" + phaseAbun + "\t" + phaseScore + "\t" + recordMarker + "\t" + pvalue
                dic['pf'][geneid][key] = value
                key1 = geneid + '\t' + strand + '\t' + pos + "\t" + seqLen + "\t" + recordMarker
                odic['pf'][geneid][key1] = value
                key2 = geneid + '\t' + strand + '\t' + pos + "\t" + seqLen + "\t" + recordMarker
                oodic['pf'][geneid][recordMarker][key2] = value
    return dic, odic, oodic

def SplitClusterIsland(cluster, island_number, phase_length):
    dic = defauldict_list()
    last = cluster[0]
    count = 1
    for i in cluster:
        if i - last <= island_number * phase_length:
            dic[count].append(i)
            last = i
        else:
            count += 1
            last = i
            dic[count].append(i)
    return dic

def SplitCluster(final_cluster: list, island_number=5, phase_length=21):
    """split list into subcluster

    Parameters
    ----------
    final_cluster : list
        sorted unique list
    island_number : int, optional
        island number, by default 5
    phase_length : int, optional
        phase length, by default 21

    Returns
    -------
    dict
        {
            cluster1:[...],
            cluster2:[...],
            ...
        }
    """
    if phase_length == 21:
        condition_list = [0, 2, 19]
    elif phase_length ==24:
        condition_list = [0, 2, 22]
    index = 0
    cluster = 1
    list_ = final_cluster[::-1]
    duplication = {}
    elem = list_.pop()
    index += 1
    if cluster not in duplication:
        duplication[cluster] = [elem]
    for i in final_cluster[index:]:
        if (i - elem)%phase_length in condition_list and (i-elem) <= island_number * phase_length + 2:
            elem = list_.pop()
            duplication[cluster].append(elem)
        else:
            cluster += 1
            elem = list_.pop()
            if cluster not in duplication:
                duplication[cluster] = [elem]
    return duplication

def split_dic(dic):
    """split dict according to the geneid

    Arguments:
        dic {file} -- dictionary
    """
    return_dic = {}
    for k,v in dic.items():
        geneid = k.split("\t")[0]
        if geneid not in return_dic.keys():
            return_dic[geneid] = {}
        if k not in return_dic[geneid].keys():
            return_dic[geneid][k] = v
    return return_dic

def SecondScaning(final_allsiRNA: dict, setting: tuple):
    """second scaning

    Parameters
    ----------
    final_allsiRNA : dict
        first scaning output dict
    setting : tuple
        setting = (final_cluster, final_strand, phase_number, phase_length, pvalue_cutoff, min_read_num)

    Returns
    -------
    tuple
        (out_phasiRNA, out_allsiRNA, out_phasiRNA_cluster, out_allsiRNA_cluster)
    """
    final_gene = {}
    duplication_phasiRNA = nestedDic()
    duplication_allsiRNA = nestedDic()
    record_marker = 0
    final_cluster = setting[0]
    final_strand = setting[1]
    phase_number = setting[2]
    phase_length = setting[3]
    pvalue_cutoff = setting[4]
    min_read_num = setting[5]
    window_length = phase_length*11
    out_allsiRNA = nestedDic()
    out_phasiRNA = nestedDic()
    out_phasiRNA_cluster = nestedDic()
    out_allsiRNA_cluster = nestedDic()
    out_duplication_phasiRNA_cluster = nestedDic()
    out_duplication_allsiRNA_cluster = nestedDic()

    for k in final_allsiRNA:
        max_readn = 0
        sorted_list = sorted(list(set(final_cluster[k])))
        candidate_cluster = SplitCluster(sorted_list, phase_length)
        for j in candidate_cluster:
            allsiRNA = []
            phasiRNA = []
            tag = False
            pair_pos = PairPos(list(range(len(candidate_cluster[j]))), candidate_cluster[j])
            for e in list(sorted(pair_pos.keys()))[::-1]:
                if e >= 4:
                    for length in list(sorted(pair_pos[e].keys()))[::-1]:
                        for tuple_ in pair_pos[length][e]:
                            total_n = 0
                            total_k = 0
                            pvalue = 0
                            if tuple_[0] > tuple_[1]:
                                start_pos = candidate_cluster[j][tuple_[1]]
                                tmp_pos = candidate_cluster[j][tuple_[0]]
                            else:
                                start_pos = candidate_cluster[j][tuple_[0]]
                                tmp_pos = candidate_cluster[j][tuple_[1]]
                            remainder = (tmp_pos - start_pos)%phase_length
                            if remainder == 0:
                                end_pos = tmp_pos
                            elif remainder == 2:
                                end_pos = tmp_pos + 19
                            elif remainder == 19:
                                end_pos = tmp_pos + 2
                            strand = final_strand[k][start_pos]
                            if strand == '+':
                                for i in range(start_pos, end_pos):
                                    query = k + '\t' + '+' + str(i)
                                    if query in final_allsiRNA[k]:
                                        total_n += 1
                                        allsiRNA.append(str(final_allsiRNA[k][k][query]))
                                for i in range(start_pos - 2, end_pos - 2):
                                    query = k + '\t' + '-' + str(i)
                                    if query in final_allsiRNA[k]:
                                        total_n += 1
                                        allsiRNA.append(str(final_allsiRNA[k][k][query]))
                                for i in range(start_pos, end_pos, phase_length):
                                    query = k + '\t' + '+' + '\t' + str(i)
                                    if query in final_allsiRNA[k] and final_allsiRNA[k][query].split('\t')[6] == str(phase_length):
                                        total_k += 1
                                        phasiRNA.append(final_allsiRNA[k][query])
                                for i in range(start_pos - 2, end_pos - 2, phase_length):
                                    query = k + '\t' + '-' + '\t' + str(i)
                                    if query in final_allsiRNA[k] and final_allsiRNA[k][query].split('\t')[6] == str(phase_length):
                                        total_k += 1
                                        phasiRNA.append(final_allsiRNA[k][query])
                            elif strand == '-':
                                for i in range(start_pos, end_pos):
                                    query = k + '\t' + '-' + '\t' + str(i)
                                    if query in final_allsiRNA[k]:
                                        total_n += 1
                                        allsiRNA.append(str(final_allsiRNA[k][query]))
                                for i in range(start_pos + 2, end_pos + 2):
                                    query = k + '\t' + '+' + '\t' + str(i)
                                    if query in final_allsiRNA[k]:
                                        total_n += 1
                                        allsiRNA.append(str(final_allsiRNA[k][query]))
                                for i in range(start_pos, end_pos, phase_length):
                                    query = k + '\t' + '-' + '\t' + str(i)
                                    if query in final_allsiRNA[k] and final_allsiRNA[k][query].split('\t')[6] == str(phase_length):
                                        total_k += 1
                                        phasiRNA.append(str(final_allsiRNA[k][query]))
                                for i in range(start_pos + 2, end_pos + 2, phase_length):
                                    query = k + '\t' + '+' + '\t' + str(i)
                                    if query in final_allsiRNA[k] and final_allsiRNA[k][query].split('\t')[6] == str(phase_length):
                                        total_k += 1
                                        phasiRNA.append(str(final_allsiRNA[k][query]))

                            # phase_number and pvalue cutoff
                            if total_k >= phase_number:
                                for i in range(total_k, phase_length):
                                    pvalue += float(binomial((end_pos - start_pos) * 2 - 1 - phase_length, total_n - i) * binomial(phase_length, i) / binomial((end_pos - start_pos) * 2 - 1, total_n))
                                if pvalue <= pvalue_cutoff and pvalue != 0.0:
                                    for i in range(len(phasiRNA)):
                                        value_list = phasiRNA[i].split('\t')
                                        sRNA_no = float(value_list[3])
                                        if sRNA_no > max_readn:
                                            max_readn = sRNA_no

                                    if min_read_num == 0 or max_readn >= min_read_num:  # min_read_num cutoff
                                        record_marker += 1
                                        tag = True
                                        for i in range(len(phasiRNA)):
                                            value_list = phasiRNA[i].split('\t')
                                            geneid = value_list[0]
                                            strand = value_list[1]
                                            readpos = int(value_list[2])
                                            if geneid in final_gene:
                                                if pvalue < float(final_gene[geneid].split('\t')[1]):
                                                    final_gene[geneid] = str(max_readn) + '\t' + str(pvalue)
                                                out_phasiRNA[geneid][record_marker][readpos] = "\t".join(phasiRNA[i].split("\t")[:-2]) + "\t" + str(pvalue)
                                            else:
                                                final_gene[geneid] = str(max_readn) + '\t' + str(pvalue)
                                                out_phasiRNA[geneid][record_marker][readpos] = "\t".join(phasiRNA[i].split("\t")[:-2]) + "\t" + str(pvalue)
                                        start_pos = list(sorted(out_phasiRNA[geneid][record_marker].keys()))[0]
                                        end_pos = list(sorted(out_phasiRNA[geneid][record_marker].keys()))[-1]
                                        out_phasiRNA_cluster[geneid][record_marker] = [start_pos, end_pos]

                                        for i in range(len(allsiRNA)):
                                            value_list = allsiRNA[i].split('\t')
                                            geneid = value_list[0]
                                            readpos = int(value_list[2])
                                            query = geneid + '\t' + '+' + '\t' + str(i)
                                            out_allsiRNA[geneid][record_marker][readpos] = "\t".join(allsiRNA[i].split("\t")[:-2]) + "\t" + str(pvalue)
                                        start_pos = list(sorted(out_allsiRNA[geneid][record_marker].keys()))[0]
                                        end_pos = list(sorted(out_allsiRNA[geneid][record_marker].keys()))[-1]
                                        out_allsiRNA_cluster[geneid][record_marker] = [start_pos, end_pos]
                            if tag:
                                break
                        if tag:
                            break
                if tag:
                    break


            for e in list(sorted(pair_pos.keys()))[::-1]:
                if e >= 4:
                    for length in list(sorted(pair_pos[e].keys())):
                        for tuple_ in pair_pos[length][e]:
                            total_n = 0
                            total_k = 0
                            pvalue = 0
                            if tuple_[0] > tuple_[1]:
                                start_pos = candidate_cluster[j][tuple_[1]]
                                tmp_pos = candidate_cluster[j][tuple_[0]]
                            else:
                                start_pos = candidate_cluster[j][tuple_[0]]
                                tmp_pos = candidate_cluster[j][tuple_[1]]
                            remainder = (tmp_pos - start_pos)%phase_length
                            if remainder == 0:
                                end_pos = tmp_pos
                            elif remainder == 2:
                                end_pos = tmp_pos + 19
                            elif remainder == 19:
                                end_pos = tmp_pos + 2
                            strand = final_strand[k][start_pos]
                            if strand == '+':
                                for i in range(start_pos, window_length + start_pos):
                                    query = k + '\t' + '+' + str(i)
                                    if query in final_allsiRNA[k]:
                                        total_n += 1
                                        allsiRNA.append(str(final_allsiRNA[k][k][query]))
                                for i in range(start_pos - 2, start_pos + window_length - 2):
                                    query = k + '\t' + '-' + str(i)
                                    if query in final_allsiRNA[k]:
                                        total_n += 1
                                        allsiRNA.append(str(final_allsiRNA[k][k][query]))
                                for i in range(start_pos, start_pos + window_length, phase_length):
                                    query = k + '\t' + '+' + '\t' + str(i)
                                    if query in final_allsiRNA[k] and final_allsiRNA[k][query].split('\t')[6] == str(phase_length):
                                        total_k += 1
                                        phasiRNA.append(final_allsiRNA[k][query])
                                for i in range(start_pos - 2, start_pos + window_length - 2, phase_length):
                                    query = k + '\t' + '-' + '\t' + str(i)
                                    if query in final_allsiRNA[k] and final_allsiRNA[k][query].split('\t')[6] == str(phase_length):
                                        total_k += 1
                                        phasiRNA.append(final_allsiRNA[k][query])
                            elif strand == '-':
                                for i in range(start_pos, start_pos + window_length):
                                    query = k + '\t' + '-' + '\t' + str(i)
                                    if query in final_allsiRNA[k]:
                                        total_n += 1
                                        allsiRNA.append(str(final_allsiRNA[k][query]))
                                for i in range(start_pos + 2, start_pos + window_length + 2):
                                    query = k + '\t' + '+' + '\t' + str(i)
                                    if query in final_allsiRNA[k]:
                                        total_n += 1
                                        allsiRNA.append(str(final_allsiRNA[k][query]))
                                for i in range(start_pos, start_pos + window_length, phase_length):
                                    query = k + '\t' + '-' + '\t' + str(i)
                                    if query in final_allsiRNA[k] and final_allsiRNA[k][query].split('\t')[6] == str(phase_length):
                                        total_k += 1
                                        phasiRNA.append(str(final_allsiRNA[k][query]))
                                for i in range(start_pos + 2, start_pos + window_length + 2, phase_length):
                                    query = k + '\t' + '+' + '\t' + str(i)
                                    if query in final_allsiRNA[k] and final_allsiRNA[k][query].split('\t')[6] == str(phase_length):
                                        total_k += 1
                                        phasiRNA.append(str(final_allsiRNA[k][query]))

                            # phase_number and pvalue cutoff
                            if total_k >= phase_number:
                                for i in range(total_k, phase_length):
                                    pvalue += float(binomial((window_length) * 2 - 1 - phase_length, total_n - i) * binomial(phase_length, i) / binomial((window_length) * 2 - 1, total_n))
                                if pvalue <= pvalue_cutoff and pvalue != 0.0:
                                    for i in range(len(phasiRNA)):
                                        value_list = phasiRNA[i].split('\t')
                                        sRNA_no = float(value_list[3])
                                        if sRNA_no > max_readn:
                                            max_readn = sRNA_no

                                    if min_read_num == 0 or max_readn >= min_read_num:  # min_read_num cutoff
                                        record_marker += 1
                                        tag = True
                                        for i in range(len(phasiRNA)):
                                            value_list = phasiRNA[i].split('\t')
                                            geneid = value_list[0]
                                            strand = value_list[1]
                                            readpos = int(value_list[2])
                                            if geneid in final_gene:
                                                if pvalue < float(final_gene[geneid].split('\t')[1]):
                                                    final_gene[geneid] = str(max_readn) + '\t' + str(pvalue)
                                                out_phasiRNA[geneid][record_marker][readpos] = "\t".join(phasiRNA[i].split("\t")[:-2]) + "\t" + str(pvalue)
                                            else:
                                                final_gene[geneid] = str(max_readn) + '\t' + str(pvalue)
                                                out_phasiRNA[geneid][record_marker][readpos] = "\t".join(phasiRNA[i].split("\t")[:-2]) + "\t" + str(pvalue)
                                        start_pos = list(sorted(out_phasiRNA[geneid][record_marker].keys()))[0]
                                        out_phasiRNA_cluster[geneid][record_marker] = [start_pos, start_pos + window_length]

                                        for i in range(len(allsiRNA)):
                                            value_list = allsiRNA[i].split('\t')
                                            geneid = value_list[0]
                                            readpos = int(value_list[2])
                                            query = geneid + '\t' + '+' + '\t' + str(i)
                                            out_allsiRNA[geneid][record_marker][readpos] = "\t".join(allsiRNA[i].split("\t")[:-2]) + "\t" + str(pvalue)
                                        start_pos = list(sorted(out_allsiRNA[geneid][record_marker].keys()))[0]
                                        out_allsiRNA_cluster[geneid][record_marker] = [start_pos, start_pos + window_length]
                            if tag:
                                break
                        if tag:
                            break
                if tag:
                    break

    return (out_phasiRNA, out_allsiRNA, out_phasiRNA_cluster, out_allsiRNA_cluster)

def PairPos(list_: list, list1_):
    """generating PairPos

    Parameters
    ----------
    list_ : list
        index_list, e.g. [0, 1, 2, 3, 4]

    Returns
    -------
    tuple list
        e.g. [(0,4), (0,3), (4,1), (4,2), (1,3), (1,2), (3,2)]
    """
    out_list = []
    final_out_list = TwoDepDic()
    rev_list = list_[::-1]
    for i in rev_list:
        length = i + 1
        for index in list_:
            distance = length - index
            if distance - 1 in list_ and distance >= 4:
                out_list.append((index, length-1))
    for i in out_list:
        length = i[1] - i[0] + 1
        distance = list1_[i[1]] - list1_[i[0]]
        final_out_list[length][distance].append(i)
    

    return final_out_list

def LoadAnnoData(anno):
    dic = nestedDic()
    with open(anno, 'r') as fn:
        for line in fn:
            l = line.strip().split("\t")
            id_ = l[0]
            anno = l[1]
            dic[id_] = anno

def AnnoGdna(dic_: dict, refgff3: dict, type_):
    """ add annotation for gdna based phasiRNA cluster

    Parameters
    ----------
    dic_ : dict
        out_phasiRNA_cluster[geneid][record_marker] = [start_pos, end_pos]
    refgff3 : dict
        
    Returns
    -------
        out_cluster_anno
            out_cluster_anno[chr][cluster].append(exonid)
    """
    out_cluster_anno = nestedDic()
    if type_ == 'h':
        final_out_cluster_anno = TwoDepDic()
    elif type_ == 'p':
        final_out_cluster_anno = ThreeDepDic()
    if type_ == 'h':
        for chr in dic_:
            if chr in refgff3:
                for geneid in refgff3[chr]:
                    if geneid == 'attribute':
                        continue
                    for cluster in dic_[chr]:
                        if cluster not in out_cluster_anno[chr]:
                            out_cluster_anno[chr][cluster] = []
                        if IsOverlap(dic_[chr][cluster], refgff3[chr][geneid]['attribute']):
                            for rnaid in refgff3[chr][geneid]:
                                if rnaid == 'attribute':
                                    continue
                                if IsOverlap(dic_[chr][cluster], refgff3[chr][geneid][rnaid]['attribute']):
                                    introns = []
                                    start = 0
                                    for type_ in refgff3[chr][geneid][rnaid]:
                                        if type_ == 'exon':
                                            for exonid in refgff3[chr][geneid][rnaid]['exon']:
                                                if start == 0:
                                                    start = refgff3[chr][geneid][rnaid]['exon'][exonid]['attribute'][1] + 1
                                                else:
                                                    introns.append([start, refgff3[chr][geneid][rnaid]['exon'][exonid]['attribute'][0] - 1])
                                                    start = refgff3[chr][geneid][rnaid]['exon'][exonid]['attribute'][1] + 1
                                                if IsOverlap(dic_[chr][cluster], refgff3[chr][geneid][rnaid]['exon'][exonid]['attribute']):
                                                    out_cluster_anno[chr][cluster].append(exonid)
                                            for intron_idx in range(0, len(introns)):
                                                if IsOverlap(dic_[chr][cluster], introns[intron_idx]):
                                                    out_cluster_anno[chr][cluster].append(rnaid+'_intron_'+str(intron_idx + 1))
                                        elif type_ == 'CDS':
                                            for exonid in refgff3[chr][geneid][rnaid]['CDS']:
                                                if IsOverlap(dic_[chr][cluster], refgff3[chr][geneid][rnaid]['CDS'][exonid]['attribute']):
                                                    out_cluster_anno[chr][cluster].append(exonid)

        for chr in out_cluster_anno:
            for record_marker in out_cluster_anno[chr]:
                if len(out_cluster_anno[chr][record_marker]) > 0:
                    coor = (dic_[chr][record_marker][0], dic_[chr][record_marker][1])
                    final_out_cluster_anno[chr][coor] = out_cluster_anno[chr][record_marker]
                else:
                    coor = (dic_[chr][record_marker][0], dic_[chr][record_marker][1])
                    final_out_cluster_anno[chr][coor] = ['Intergenic']
    elif type_ == 'p':
        for chromosome in list(dic_.keys()):
            chr__ = chromosome.split('_')[:-1]
            chr_ = "_".join(chr__)
            for record in dic_[chromosome]:
                if chr_ in refgff3:
                    out_cluster_anno[chromosome][record] = []
                    for geneid in refgff3[chr_]:
                        if geneid == 'attribute':
                            continue
                        if IsOverlap(dic_[chromosome][record], refgff3[chr_][geneid]['attribute']):
                            for rnaid in refgff3[chr_][geneid]:
                                if rnaid == 'attribute':
                                    continue
                                if IsOverlap(dic_[chromosome][record], refgff3[chr_][geneid][rnaid]['attribute']):
                                    introns = []
                                    start = 0
                                    for type_ in refgff3[chr_][geneid][rnaid]:
                                        if type_ == 'exon':
                                            for exonid in refgff3[chr_][geneid][rnaid]['exon']:
                                                if start == 0:
                                                    start = refgff3[chr_][geneid][rnaid]['exon'][exonid]['attribute'][1] + 1
                                                else:
                                                    introns.append([start, refgff3[chr_][geneid][rnaid]['exon'][exonid]['attribute'][0] - 1])
                                                    start = refgff3[chr_][geneid][rnaid]['exon'][exonid]['attribute'][1] + 1
                                                if IsOverlap(dic_[chromosome][record], refgff3[chr_][geneid][rnaid]['exon'][exonid]['attribute']):
                                                    out_cluster_anno[chromosome][record].append(exonid)
                                            for intron_idx in range(0, len(introns)):
                                                if IsOverlap(dic_[chromosome][record], introns[intron_idx]):
                                                    out_cluster_anno[chromosome][record].append(rnaid+'_intron_'+str(intron_idx + 1))
                                        elif type_ == 'CDS':
                                            for exonid in refgff3[chr_][geneid][rnaid]['CDS']:
                                                if IsOverlap(dic_[chromosome][record], refgff3[chr_][geneid][rnaid]['CDS'][exonid]['attribute']):
                                                    out_cluster_anno[chromosome][record].append(exonid)

        for chr_ in out_cluster_anno:
            for record in out_cluster_anno[chr_]:
                if len(out_cluster_anno[chr_][record]) > 0:
                    coor = (dic_[chr_][record][0], dic_[chr_][record][1])
                    final_out_cluster_anno[chr_][record][coor] = out_cluster_anno[chr_][record]
                else:
                    coor = (dic_[chr_][record][0], dic_[chr_][record][1])
                    final_out_cluster_anno[chr_][record][coor] = ['Intergenic']

    return final_out_cluster_anno

def IsOverlap(list1: list, list2: list):
    """ is overlap?

    Parameters
    ----------
    list1 : list
        e.g. [1,100]
    list2 : list
        e.g. [50,150]

    Returns
    -------
    bool
        bool value
    """
    a = list1[0]
    b = list1[-1]
    c = list2[0]
    d = list2[-1]
    list1 = set(list(range(a, b+1)))
    list2 = set(list(range(c, d+1)))
    test = list1 & list2
    if len(test) > 0:
        return True
    else:
        return False

def ParseGff3(file: str):
    """parse gff3

    Parameters
    ----------
    file : file
        _description_

    Returns
    -------
    dict
    """
    with open(file, 'r') as fn:
        dic = defaultdict()
        grandParentList = ['region']
        pararentList = ['gene', 'cDNA_match', 'match', 'pseudogene']
        sonList = ['lnc_RNA', 'transcript', 'mRNA', 'primary_transcript', 'snoRNA', 'tRNA', 'rRNA', 'snRNA']
        grandSonList = ['exon', 'CDS', 'miRNA', 'intron', 'three_prime_UTR', 'five_prime_UTR']
        for line in fn:
            if line.startswith("#"):
                continue
            l = line.strip().split("\t")
            feature = l[2]
            start = int(l[3])
            end = int(l[4])
            id = l[8].split(';')[0]
            # key = id + "\t" + str(start) + "\t" + str(end)
            if feature in grandParentList:
                chr_id = id.split('=')[1].split(':')[0]
                dic[chr_id] = {'attribute': [start, end]}
                curGrandParentFeature = chr_id
            elif feature in pararentList:
                geneid = id.split('-')[1]
                dic[curGrandParentFeature][geneid] = {'attribute': [start, end]}
                curParentFeature = geneid
            elif feature in sonList:
                rnaid = id.split('-')[1]
                dic[curGrandParentFeature][curParentFeature][rnaid] = {'attribute': [start, end]}
                curSonFeature = rnaid
            elif feature in grandSonList:
                exonid = id.split('-')[1:]
                if len(exonid) > 1:
                    exonid = '-'.join(exonid)
                elif len(exonid) == 1:
                    exonid = exonid[0]
                if feature == 'exon':
                    try:
                        if 'exon' not in dic[curGrandParentFeature][curParentFeature][curSonFeature]:
                            dic[curGrandParentFeature][curParentFeature][curSonFeature]['exon'] = {}
                        else:
                            dic[curGrandParentFeature][curParentFeature][curSonFeature]['exon'][exonid] = {'attribute': [start, end]}
                    except KeyError:
                        try:
                            curSonFeature = curParentFeature + 'rna'
                            dic[curGrandParentFeature][curParentFeature][curSonFeature] = {'attribute': [start, end]}
                            if 'exon' not in dic[curGrandParentFeature][curParentFeature][curSonFeature]:
                                dic[curGrandParentFeature][curParentFeature][curSonFeature]['exon'] = {}
                            else:
                                dic[curGrandParentFeature][curParentFeature][curSonFeature]['exon'][exonid] = {'attribute': [start, end]}
                        except KeyError:
                            curParentFeature = curGrandParentFeature + 'gene'
                            dic[curGrandParentFeature][curParentFeature] = {'attribute': [start, end]}
                            curSonFeature = curParentFeature + 'rna'
                            dic[curGrandParentFeature][curParentFeature][curSonFeature] = {'attribute': [start, end]}
                            if 'exon' not in dic[curGrandParentFeature][curParentFeature][curSonFeature]:
                                dic[curGrandParentFeature][curParentFeature][curSonFeature]['exon'] = {}
                            else:
                                dic[curGrandParentFeature][curParentFeature][curSonFeature]['exon'][exonid] = {'attribute': [start, end]}
                elif feature == 'CDS':
                    exonid = exonid + '-' + str(start) + '-' + str(end)
                    try:
                        if 'CDS' not in dic[curGrandParentFeature][curParentFeature][curSonFeature]:
                            dic[curGrandParentFeature][curParentFeature][curSonFeature]['CDS'] = {}
                        else:
                            dic[curGrandParentFeature][curParentFeature][curSonFeature]['CDS'][exonid] = {'attribute': [start, end]}
                    except KeyError:
                        try:
                            curSonFeature = curParentFeature + 'rna'
                            dic[curGrandParentFeature][curParentFeature][curSonFeature] = {'attribute': [start, end]}
                            if 'CDS' not in dic[curGrandParentFeature][curParentFeature][curSonFeature]:
                                dic[curGrandParentFeature][curParentFeature][curSonFeature]['CDS'] = {}
                            else:
                                dic[curGrandParentFeature][curParentFeature][curSonFeature]['CDS'][exonid] = {'attribute': [start, end]}
                        except KeyError:
                            curParentFeature = curGrandParentFeature + 'gene'
                            dic[curGrandParentFeature][curParentFeature] = {'attribute': [start, end]}
                            curSonFeature = curParentFeature + 'rna'
                            dic[curGrandParentFeature][curParentFeature][curSonFeature] = {'attribute': [start, end]}
                            if 'CDS' not in dic[curGrandParentFeature][curParentFeature][curSonFeature]:
                                dic[curGrandParentFeature][curParentFeature][curSonFeature]['CDS'] = {}
                            else:
                                dic[curGrandParentFeature][curParentFeature][curSonFeature]['CDS'][exonid] = {'attribute': [start, end]}
                elif feature == 'miRNA':
                    try:
                        if 'miRNA' not in dic[curGrandParentFeature][curParentFeature][curSonFeature]:
                            dic[curGrandParentFeature][curParentFeature][curSonFeature]['miRNA'] = {}
                        else:
                            dic[curGrandParentFeature][curParentFeature][curSonFeature]['miRNA'][exonid] = {'attribute': [start, end]}
                    except KeyError:
                        try:
                            curSonFeature = curParentFeature + 'rna'
                            dic[curGrandParentFeature][curParentFeature][curSonFeature] = {'attribute': [start, end]}
                            if 'miRNA' not in dic[curGrandParentFeature][curParentFeature][curSonFeature]:
                                dic[curGrandParentFeature][curParentFeature][curSonFeature]['miRNA'] = {}
                            else:
                                dic[curGrandParentFeature][curParentFeature][curSonFeature]['miRNA'][exonid] = {'attribute': [start, end]}
                        except KeyError:
                            curParentFeature = curGrandParentFeature + 'gene'
                            dic[curGrandParentFeature][curParentFeature] = {'attribute': [start, end]}
                            curSonFeature = curParentFeature + 'rna'
                            dic[curGrandParentFeature][curParentFeature][curSonFeature] = {'attribute': [start, end]}
                            if 'miRNA' not in dic[curGrandParentFeature][curParentFeature][curSonFeature]:
                                dic[curGrandParentFeature][curParentFeature][curSonFeature]['miRNA'] = {}
                            else:
                                dic[curGrandParentFeature][curParentFeature][curSonFeature]['miRNA'][exonid] = {'attribute': [start, end]}
                elif feature == 'intron':
                    try:
                        if 'intron' not in dic[curGrandParentFeature][curParentFeature][curSonFeature]:
                            dic[curGrandParentFeature][curParentFeature][curSonFeature]['intron'] = {}
                        else:
                            dic[curGrandParentFeature][curParentFeature][curSonFeature]['intron'][exonid] = {'attribute': [start, end]}
                    except KeyError:
                        try:
                            curSonFeature = curParentFeature + 'rna'
                            dic[curGrandParentFeature][curParentFeature][curSonFeature] = {'attribute': [start, end]}
                            if 'intron' not in dic[curGrandParentFeature][curParentFeature][curSonFeature]:
                                dic[curGrandParentFeature][curParentFeature][curSonFeature]['intron'] = {}
                            else:
                                dic[curGrandParentFeature][curParentFeature][curSonFeature]['intron'][exonid] = {'attribute': [start, end]}
                        except KeyError:
                            curParentFeature = curGrandParentFeature + 'gene'
                            dic[curGrandParentFeature][curParentFeature] = {'attribute': [start, end]}
                            curSonFeature = curParentFeature + 'rna'
                            dic[curGrandParentFeature][curParentFeature][curSonFeature] = {'attribute': [start, end]}
                            if 'intron' not in dic[curGrandParentFeature][curParentFeature][curSonFeature]:
                                dic[curGrandParentFeature][curParentFeature][curSonFeature]['intron'] = {}
                            else:
                                dic[curGrandParentFeature][curParentFeature][curSonFeature]['intron'][exonid] = {'attribute': [start, end]}
    return dic

def GetPhaseScoreCluster(pcluster_dic, phase_length):
    dic = nestedDic() 
    for i in pcluster_dic:
        for record_marker in pcluster_dic[i]:
            maxium = max(pcluster_dic[i][record_marker])
            minium = min(pcluster_dic[i][record_marker])
            dic[i][record_marker] = [minium, maxium + phase_length]
    return dic

def ConvertFormat(hout, new_allsiRNA):
    dic = FiveDepDic()
    for i in ['hc', 'hg', 'hf']:
        for j in hout[i]['allsiRNA']:
            for gene in hout[i]['allsiRNA']:
                for cluster in hout[i]['allsiRNA'][gene]:
                    for pos in hout[i]['allsiRNA'][gene][cluster]:
                        pvalue = hout[i]['allsiRNA'][gene][cluster][pos].split("\t")[-1] 
                        break
                    if gene in new_allsiRNA[i]:
                        for key in new_allsiRNA[i][gene]:
                            l = key.split("\t")
                            strand, pos, seqLen = l[1], l[2], l[3]
                            if int(pos) in hout[i]['allsiRNA'][gene][cluster]:
                                value = "\t".join(new_allsiRNA[i][gene][key].split("\t")[:-2]) + "\t" + pvalue
                                if value not in dic[i]['allsiRNA'][gene][cluster][pos]:
                                    dic[i]['allsiRNA'][gene][cluster][pos].append(value)
                        # for pos_ in dic[i]['allsiRNA'][gene][cluster]:
                        #     dic[i]['allsiRNA'][gene][cluster][pos_] = list(set(dic[i]['allsiRNA'][gene][cluster][pos_]))
    return dic

def Write(hout, pout, hanno_info, panno_info, p_phasiRNA_cluster, fo_phasiRNA, fo_allsiRNA):
    intergration = nestedDic()
    Vprint('start writing data...')
    fo_phasiRNA.write('>Hypergeometric' + "\n")
    fo_phasiRNA.write('#H cDNA based result' + "\n")
    record_marker = 0
    for gene in hout['hc']['phasiRNA']:
        for cluster in hout['hc']['phasiRNA'][gene]:
            coor = (hout['hc']['phasiRNA_cluster'][gene][cluster][0], hout['hc']['phasiRNA_cluster'][gene][cluster][1])
            record_marker += 1
            out_marker = "HC" + str(record_marker)
            for pos in hout['hc']['phasiRNA'][gene][cluster]:
                pvalue = float(hout['hc']['phasiRNA'][gene][cluster][pos].split("\t")[-1])
                intergration[gene]['H']['C'][coor] = ['-', '-', pvalue]
                fo_phasiRNA.write(hout['hc']['phasiRNA'][gene][cluster][pos] + "\t" + "-" + "\t" + out_marker + "\n")

    fo_phasiRNA.write('#H gDNA based result' + "\n")
    record_marker = 0
    for gene in hout['hg']['phasiRNA']:
        for cluster in hout['hg']['phasiRNA'][gene]:
            coor = (hout['hg']['phasiRNA_cluster'][gene][cluster][0], hout['hg']['phasiRNA_cluster'][gene][cluster][1])
            record_marker += 1
            out_marker = "HG" + str(record_marker)
            for pos in hout['hg']['phasiRNA'][gene][cluster]:
                anno = hanno_info[gene][coor]
                pvalue = float(hout['hg']['phasiRNA'][gene][cluster][pos].split("\t")[-1])
                chr_ = hout['hg']['phasiRNA'][gene][cluster][pos].split("\t")[0]                       
                for i in hanno_info[gene][coor]:
                    i = i.split("-")[0]
                    if i == 'Intergenic':
                        intergration[chr_ + '_' + i]['H']['G'][coor] = ['-', '-', pvalue]
                    else:
                        intergration[i]['H']['G'][coor] = ['-', '-', pvalue]
                annotation = ';'.join(hanno_info[gene][coor])
                fo_phasiRNA.write(hout['hg']['phasiRNA'][gene][cluster][pos] + "\t" + annotation + "\t" + out_marker + "\n")

    fo_phasiRNA.write('#H FLNC based result' + "\n")
    record_marker = 0
    for gene in hout['hf']['phasiRNA']:
        for cluster in hout['hf']['phasiRNA'][gene]:
            coor = (hout['hf']['phasiRNA_cluster'][gene][cluster][0], hout['hf']['phasiRNA_cluster'][gene][cluster][1])
            record_marker += 1
            out_marker = "HF" + str(record_marker)
            for pos in hout['hf']['phasiRNA'][gene][cluster]:
                pvalue = float(hout['hf']['phasiRNA'][gene][cluster][pos].split("\t")[-1])
                intergration[gene]['H']['F'][coor] = ['-', '-', pvalue]
                fo_phasiRNA.write(hout['hf']['phasiRNA'][gene][cluster][pos] + "\t" + "-" + "\t" + out_marker + "\n")
    
    fo_phasiRNA.write('>PhaseScore' + "\n")
    fo_phasiRNA.write('#P cDNA based result' + "\n")
    record_marker = 0
    for gene in pout['pc']['phasiRNA']:
        for record in pout['pc']['phasiRNA'][gene]:
            coor = (p_phasiRNA_cluster[0][gene][record][0], p_phasiRNA_cluster[0][gene][record][1])
            record_marker += 1 
            out_marker = "PC" + str(record_marker)
            for key in pout['pc']['phasiRNA'][gene][record]:
                phaseRatio = float(pout['pc']['phasiRNA'][gene][record][key].split("\t")[-6])
                phaseScore = float(pout['pc']['phasiRNA'][gene][record][key].split("\t")[-3])
                intergration[gene]['P']['C'][coor] = [phaseRatio, phaseScore, '-']
                tmp  = "\t".join(pout['pc']['phasiRNA'][gene][record][key].split("\t")[:-2]) + "\t" + pout['pc']['phasiRNA'][gene][record][key].split("\t")[-1]
                fo_phasiRNA.write(tmp + "\t" + "-" + "\t" + out_marker + "\n")

    fo_phasiRNA.write('#P gDNA based result' + "\n")
    record_marker = 0
    for gene in pout['pg']['phasiRNA']:
        for record in pout['pg']['phasiRNA'][gene]:
            coor = (p_phasiRNA_cluster[1][gene][record][0], p_phasiRNA_cluster[1][gene][record][1])
            record_marker += 1
            out_marker = "PG" + str(record_marker)
            for key in pout['pg']['phasiRNA'][gene][record]:
                phaseRatio = float(pout['pg']['phasiRNA'][gene][record][key].split("\t")[-6])
                phaseScore = float(pout['pg']['phasiRNA'][gene][record][key].split("\t")[-3])
                chr_ = pout['pg']['phasiRNA'][gene][record][key].split("\t")[0]
                anno = panno_info[gene][record][coor]
                for i in anno:
                    i = i.split('-')[0]
                    if i == 'Intergenic':
                        intergration[chr_]['P']['G'][coor] = [phaseRatio, phaseScore, '-']
                    else:
                        intergration[i]['P']['G'][coor] = [phaseRatio, phaseScore, '-']
                annotation = ";".join(panno_info[gene][record][coor])
                tmp  = "\t".join(pout['pg']['phasiRNA'][gene][record][key].split("\t")[1:-2]) + "\t" + pout['pg']['phasiRNA'][gene][record][key].split("\t")[-1]
                geneid = pout['pg']['phasiRNA'][gene][record][key].split("\t")[0].split('_')[:-1]
                geneid = "_".join(geneid)
                fo_phasiRNA.write(geneid + "\t" + tmp + "\t" + annotation + "\t" + out_marker + "\n")
    
    fo_phasiRNA.write('#P FLNC based result' + "\n")
    record_marker = 0
    for gene in pout['pf']['phasiRNA']:
        for record in pout['pf']['phasiRNA'][gene]:
            coor = (p_phasiRNA_cluster[2][gene][record][0], p_phasiRNA_cluster[2][gene][record][1])
            record_marker += 1
            out_marker = "PF" + str(record_marker)
            for key in pout['pf']['phasiRNA'][gene][record]:
                phaseRatio = float(pout['pf']['phasiRNA'][gene][record][key].split("\t")[-6])
                phaseScore = float(pout['pf']['phasiRNA'][gene][record][key].split("\t")[-3])
                intergration[gene]['P']['F'][coor] = [phaseRatio, phaseScore, '-']
                tmp  = "\t".join(pout['pf']['phasiRNA'][gene][record][key].split("\t")[:-2]) + "\t" + pout['pf']['phasiRNA'][gene][record][key].split("\t")[-1]
                fo_phasiRNA.write(tmp + "\t" + "-" + "\t" + out_marker + "\n")

    fo_allsiRNA.write('>Hypergeometric' + "\n")
    fo_allsiRNA.write('#H cDNA based result' + "\n")
    record_marker = 0
    for gene in hout['hc']['allsiRNA']:
        for cluster in hout['hc']['allsiRNA'][gene]:
            record_marker += 1
            out_marker = "HC" + str(record_marker)
            for pos in hout['hc']['allsiRNA'][gene][cluster]:
                for i1 in hout['hc']['allsiRNA'][gene][cluster][pos]:
                    fo_allsiRNA.write(i1 + "\t" + '-' + "\t" + out_marker + "\n")

    fo_allsiRNA.write('#H gDNA based result' + "\n")
    record_marker = 0
    for gene in hout['hg']['allsiRNA']:
        for cluster in hout['hg']['allsiRNA'][gene]:
            coor = (hout['hg']['phasiRNA_cluster'][gene][cluster][0], hout['hg']['phasiRNA_cluster'][gene][cluster][1])
            record_marker += 1
            out_marker = "HG" + str(record_marker)
            for pos in hout['hg']['allsiRNA'][gene][cluster]:
                anno = hanno_info[gene][coor]
                annotation = ';'.join(hanno_info[gene][coor])
                for i1 in hout['hg']['allsiRNA'][gene][cluster][pos]:
                    fo_allsiRNA.write(i1 + "\t" + annotation + "\t" + out_marker + "\n")
    
    fo_allsiRNA.write('#H FLNC based result' + "\n")
    record_marker = 0
    for gene in hout['hf']['allsiRNA']:
        for cluster in hout['hf']['allsiRNA'][gene]:
            record_marker += 1
            out_marker = "HF" + str(record_marker)
            for pos in hout['hf']['allsiRNA'][gene][cluster]:
                for i1 in hout['hf']['allsiRNA'][gene][cluster][pos]:
                    fo_allsiRNA.write(i1 + "\t" + '-' + "\t" + out_marker + "\n")

    fo_allsiRNA.write('>PhaseScore' + "\n")
    fo_allsiRNA.write('#P cDNA based result' + "\n")
    record_marker = 0
    for gene in pout['pc']['allsiRNA']:
        for record in pout['pc']['allsiRNA'][gene]:
            record_marker += 1
            out_marker = "PC" + str(record_marker)
            for key in pout['pc']['allsiRNA'][gene][record]:
                tmp  = "\t".join(pout['pc']['allsiRNA'][gene][record][key].split("\t")[:-2]) + "\t" + pout['pc']['allsiRNA'][gene][record][key].split("\t")[-1]
                fo_allsiRNA.write(tmp + "\t" + "-" + "\t" + out_marker + "\n")

    fo_allsiRNA.write('#P gDNA based result' + "\n")
    record_marker = 0
    for gene in pout['pg']['allsiRNA']:
        for record in pout['pg']['allsiRNA'][gene]:
            coor = (p_phasiRNA_cluster[1][gene][record][0], p_phasiRNA_cluster[1][gene][record][1])
            record_marker += 1
            out_marker = "PG" + str(record_marker)
            for key in pout['pg']['allsiRNA'][gene][record]:
                annotation = ";".join(panno_info[gene][record][coor])
                tmp  = "\t".join(pout['pg']['allsiRNA'][gene][record][key].split("\t")[1:-2]) + "\t" + pout['pg']['allsiRNA'][gene][record][key].split("\t")[-1]
                geneid = pout['pg']['allsiRNA'][gene][record][key].split("\t")[0].split('_')[:-1]
                geneid = "_".join(geneid)
                fo_allsiRNA.write(geneid + "\t" + tmp + "\t" + annotation + "\t" + out_marker + "\n")
    
    fo_allsiRNA.write('#P FLNC based result' + "\n")
    record_marker = 0
    for gene in pout['pf']['allsiRNA']:
        for record in pout['pf']['allsiRNA'][gene]:
            record_marker += 1
            out_marker = "PF" + str(record_marker)
            for key in pout['pf']['allsiRNA'][gene][record]:
                tmp  = "\t".join(pout['pf']['allsiRNA'][gene][record][key].split("\t")[:-2]) + "\t" + pout['pf']['allsiRNA'][gene][record][key].split("\t")[-1]
                fo_allsiRNA.write(tmp + "\t" + "-" + "\t" + out_marker + "\n")
    
    return intergration

def WriteIntergration(intergration, fo):
    fo.write(f'geneid\tmethod\treference_type\th_transcriptome_coordinate\th_genome_coordinate\th_fl_transcriptome_coordinate\tp_transcriptome_coordinate\tp_genome_coordinate\tp_genome_coordinate\
        \tc_pvalue\tg_pvalue\tf_pvalue\tc_phaseratio\tg_phaseratio\tf_phaseratio\tc_phasescore\tg_phasescore\tf_phasescore\n')
    for gene in intergration:
        methods = []
        refs = []
        coors = TwoDepDic()
        pvalues = TwoDepDic()
        phaseRatios = TwoDepDic()
        phaseScores = TwoDepDic()
        for method in intergration[gene]:
            methods.append(method)
            for ref in intergration[gene][method]:
                refs.append(ref)
                for coor in intergration[gene][method][ref]:
                    coor_ = f'({coor[0]}:{coor[1]})'
                    pvalue = str(intergration[gene][method][ref][coor][2])
                    phaseRatio = str(intergration[gene][method][ref][coor][0])
                    phaseScore = str(intergration[gene][method][ref][coor][1])
                    if ref == 'C': 
                        coors[method]['C'].append(coor_)
                        pvalues[method]['C'].append(pvalue)
                        phaseRatios[method]['C'].append(phaseRatio)
                        phaseScores[method]['C'].append(phaseScore)
                    elif ref == 'G':
                        coors[method]['G'].append(coor_)
                        pvalues[method]['G'].append(pvalue)
                        phaseRatios[method]['G'].append(phaseRatio)
                        phaseScores[method]['G'].append(phaseScore)
                    elif ref == 'F':
                        coors[method]['F'].append(coor_)
                        pvalues[method]['F'].append(pvalue)
                        phaseRatios[method]['F'].append(phaseRatio)
                        phaseScores[method]['F'].append(phaseScore)

        for method in ['H', 'P']:
            if len(coors[method]['C']) == 0:
                coors[method]['C'] = ['-']
            if len(coors[method]['F']) == 0:
                coors[method]['F'] = ['-']
            if len(coors[method]['G']) == 0:
                coors[method]['G'] = ['-']
            for ref in ['C', 'G', 'F']:
                if len(pvalues[method][ref]) == 0:
                    pvalues[method][ref] = ['-']
            for ref in ['C', 'G', 'F']:
                if len(phaseRatios[method][ref]) == 0:
                    phaseRatios[method][ref] = ['-']
            for ref in ['C', 'G', 'F']:
                if len(phaseScores[method][ref]) == 0:
                    phaseScores[method][ref] = ['-']

        out_tmp = f'{gene}\t{";".join(set(methods))}\t{";".join(set(refs))}\t'
        tmp_list = []
        for method in ['H', 'P']:
            for ref in ['C', 'G', 'F']:
                tmp = []
                for i in coors[method][ref]:
                    tmp.append(i)
                tmp = ";".join(tmp)
                tmp_list.append(tmp)
        out_tmp += "\t".join(tmp_list)
        out_tmp += "\t"
        tmp_list = []
        for method in ['H']:
            for ref in ['C', 'G', 'F']:
                tmp = []
                for i in pvalues[method][ref]:
                    tmp.append(i)
                tmp = ";".join(tmp)
                tmp_list.append(tmp)
        out_tmp += "\t".join(tmp_list) 
        out_tmp += "\t"
        tmp_list = []
        for method in ["P"]:
            for ref in ['C', 'G', 'F']:
                tmp = []
                for i in phaseRatios[method][ref]:
                    tmp.append(i)
                tmp = ";".join(tmp)
                tmp_list.append(tmp)
        out_tmp += "\t".join(tmp_list) 
        out_tmp += "\t"
        tmp_list = []
        for method in ["P"]:
            for ref in ['C', 'G', 'F']:
                tmp = []
                for i in phaseScores[method][ref]:
                    tmp.append(i)
                tmp = ";".join(tmp)
                tmp_list.append(tmp)
        out_tmp += "\t".join(tmp_list) 
        out_tmp += "\n"
            
        fo.write(out_tmp)   

def GeneratingIntronInforFromGff(tmp_gff_bed):
    intron_infor = []
    with open(tmp_gff_bed, 'r') as fn:
        former_feature = ''
        former_fiter_anno = ''
        former_end = -1
        former_start = -1
        for line in fn:
            l = line.strip().split("\t")
            Chr, start, end, feature, anno, strand = l[0], int(l[1]), int(l[2]), l[3], l[4], l[5]
            filter_anno = "-".join(anno.split('-')[:-1])
            filter_order = (anno.split('-')[-1]).replace('E', '')
            if feature == 'exon':
                if former_feature == feature and former_fiter_anno == filter_anno:
                    if former_end + 1 > start - 1:
                        region = (end + 1, former_start - 1)
                    else:
                        region = (former_end + 1, start - 1)
                    if region[0] > region[1]:
                        continue
                    intron = (Chr, region, 'intron', filter_anno + '--' + str(int(filter_order) - 1), strand)
                    intron_infor.append(intron)
                    former_feature = feature
                    former_fiter_anno = filter_anno
                    former_end = end
                    former_start = start
                else:
                    former_feature = feature
                    former_fiter_anno = filter_anno
                    former_end = end
                    former_start = start
            else:
                    former_feature = feature
                    former_fiter_anno = filter_anno
                    former_end = end
                    former_start = start
    return intron_infor

def GetPhasiRNAClusterBed(hout_pout, bedout, Type, phase_length):
    fo = open(bedout, 'w')
    if Type == 'h':
        for geneid in hout_pout['hg']['phasiRNA_cluster']:
            for cluster in hout_pout['hg']['phasiRNA_cluster'][geneid]:
                start = hout_pout['hg']['phasiRNA_cluster'][geneid][cluster][0]
                end = hout_pout['hg']['phasiRNA_cluster'][geneid][cluster][-1]
                fo.write(f'{geneid}\t{start}\t{end}\n')
    elif Type == 'p':
        for geneid in hout_pout['pg']['phasiRNA']:
            for record in hout_pout['pg']['phasiRNA'][geneid]:
                region = []
                for sRNA in hout_pout['pg']['phasiRNA'][geneid][record]:
                    pos = int(sRNA.split('\t')[-1])
                    region.append(pos)
                sorted_region = sorted(region)
                geneid_ = "_".join(geneid.split('_')[:-1])
                fo.write(f'{geneid_}\t{sorted_region[0]}\t{str(int(sorted_region[-1]) + phase_length)}\n')
    fo.close()

def AnnoGdnaWithBedtools(tmp_hout_bed, tmp_gff_bed):
    fo = os.popen(f'bedtools intersect -nonamecheck -a {tmp_hout_bed} -b {tmp_gff_bed} -wao | cut -f 1,2,3,7,8,9,10')
    return fo

def WringIntronInforToGff(intron_infor, tmp_gff_bed):
    with open(tmp_gff_bed, 'a+') as fo:
        for i in intron_infor:
            string = f'{str(i[0])}\t{str(i[1][0])}\t{str(i[1][1])}\t{str(i[2])}\t{str(i[3])}\t{str(i[4])}\n'
            fo.write(string)

def ParseBedAnno(hout_anno_io):
    dic = ThreeDepDic()
    for line in hout_anno_io:
        l = line.strip().split("\t")
        Chr, start, end, feature, anno, strand, overlap = l[0], int(l[1]), int(l[2]), l[3], l[4], l[5], l[6]
        if overlap == '0':
            anno = 'Intergenic'
            feature = 'Other'
        dic[Chr][(start, end)][feature].append(anno)
    return dic

def FinalOut(tmp_out_list, relation_dic, fo, passP):
    dic = defauldict_list()
    for i in tmp_out_list:
        l = i.split("\t")
        geneid, method, ref, hcc, hgc, hfc, pcc, pgc, pfc, hcp, hgp, hfp, pcr, pgr, pfr, pcs, pgs, pfs = \
        l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10], l[11], l[12], l[13], l[14], l[15], l[16], l[17]
        if geneid == 'geneid':
            fo.write(f'cluster_number\tfeature\t' + "\t".join(l))
            continue
        feature = relation_dic[geneid] 
        if feature == "Other":
            feature = 'Intergenic'
        if passP == 'y' and 'H' not in method and feature != 'Intergenic':
            continue
        max_list = []
        for j in [hcc, hgc, hfc, pcc, pgc, pfc]:
            k = j.split(';')
            max_list.append(len(k))
        max_cluster = str(max(max_list))
        dic[feature].append([max_cluster, feature] + l)
        fo.write("\t".join([max_cluster, feature] + l))
    return dic

def FormatForWriteInterGration(intergration, list): 
    dic = nestedDic()
    relation_dic = nestedDic()
    for feature in list:
        for gene in intergration[feature]:
            for method in intergration[feature][gene]:
                for ref in intergration[feature][gene][method]:
                    for coor in intergration[feature][gene][method][ref]:
                        dic[gene][method][ref][coor] = intergration[feature][gene][method][ref][coor]
            relation_dic[gene] = feature
    return dic, relation_dic

def WriteIntergration_new(intergration):
    list_ = []
    list_.append(f'geneid\tmethod\treference_type\th_transcriptome_coordinate\th_genome_coordinate\th_fl_transcriptome_coordinate\tp_transcriptome_coordinate\tp_genome_coordinate\tp_genome_coordinate\
        \tc_pvalue\tg_pvalue\tf_pvalue\tc_phaseratio\tg_phaseratio\tf_phaseratio\tc_phasescore\tg_phasescore\tf_phasescore\n')
    for gene in intergration:
        methods = []
        refs = []
        coors = TwoDepDic()
        pvalues = TwoDepDic()
        phaseRatios = TwoDepDic()
        phaseScores = TwoDepDic()
        for method in intergration[gene]:
            methods.append(method)
            for ref in intergration[gene][method]:
                refs.append(ref)
                for coor in intergration[gene][method][ref]:
                    coor_ = f'({coor[0]}:{coor[1]})'
                    pvalue = str(intergration[gene][method][ref][coor][2])
                    phaseRatio = str(intergration[gene][method][ref][coor][0])
                    phaseScore = str(intergration[gene][method][ref][coor][1])
                    if ref == 'C': 
                        coors[method]['C'].append(coor_)
                        pvalues[method]['C'].append(pvalue)
                        phaseRatios[method]['C'].append(phaseRatio)
                        phaseScores[method]['C'].append(phaseScore)
                    elif ref == 'G':
                        coors[method]['G'].append(coor_)
                        pvalues[method]['G'].append(pvalue)
                        phaseRatios[method]['G'].append(phaseRatio)
                        phaseScores[method]['G'].append(phaseScore)
                    elif ref == 'F':
                        coors[method]['F'].append(coor_)
                        pvalues[method]['F'].append(pvalue)
                        phaseRatios[method]['F'].append(phaseRatio)
                        phaseScores[method]['F'].append(phaseScore)

        for method in ['H', 'P']:
            if len(coors[method]['C']) == 0:
                coors[method]['C'] = ['-']
            if len(coors[method]['F']) == 0:
                coors[method]['F'] = ['-']
            if len(coors[method]['G']) == 0:
                coors[method]['G'] = ['-']
            for ref in ['C', 'G', 'F']:
                if len(pvalues[method][ref]) == 0:
                    pvalues[method][ref] = ['-']
            for ref in ['C', 'G', 'F']:
                if len(phaseRatios[method][ref]) == 0:
                    phaseRatios[method][ref] = ['-']
            for ref in ['C', 'G', 'F']:
                if len(phaseScores[method][ref]) == 0:
                    phaseScores[method][ref] = ['-']

        out_tmp = f'{gene}\t{";".join(set(methods))}\t{";".join(set(refs))}\t'
        tmp_list = []
        # print coor
        for method in ['H', 'P']:
            for ref in ['C', 'G', 'F']:
                tmp = []
                for i in coors[method][ref]:
                    tmp.append(i)
                tmp = ";".join(tmp)
                tmp_list.append(tmp)
        out_tmp += "\t".join(tmp_list)
        out_tmp += "\t"
        tmp_list = []
        # print pvalue
        for method in ['H']:
            for ref in ['C', 'G', 'F']:
                tmp = []
                for i in pvalues[method][ref]:
                    tmp.append(i)
                tmp = ";".join(tmp)
                tmp_list.append(tmp)
        out_tmp += "\t".join(tmp_list) 
        out_tmp += "\t"
        tmp_list = []
        # print phase ratio
        for method in ["P"]:
            for ref in ['C', 'G', 'F']:
                tmp = []
                for i in phaseRatios[method][ref]:
                    tmp.append(i)
                tmp = ";".join(tmp)
                tmp_list.append(tmp)
        out_tmp += "\t".join(tmp_list) 
        out_tmp += "\t"
        tmp_list = []
        # print phase score
        for method in ["P"]:
            for ref in ['C', 'G', 'F']:
                tmp = []
                for i in phaseScores[method][ref]:
                    tmp.append(i)
                tmp = ";".join(tmp)
                tmp_list.append(tmp)
        out_tmp += "\t".join(tmp_list) 
        out_tmp += "\n"
            
        list_.append(out_tmp)   
    return list_

def ParseTranscriptFeature(tmp_gff_bed):
    dic = nestedDic()
    with open(tmp_gff_bed, 'r') as fn:
        for line in fn:
            l = line.strip().split("\t")
            Chr, start, end, feature, anno = l[0], l[1], l[2], l[3], l[4]
            dic[anno] = feature
    return dic

def Write_new(hout, pout, hanno, panno, p_phasiRNA_cluster, fo_phasiRNA, fo_allsiRNA, transcriptAnno, flnc_anno_dic):
    intergration_FLNC = nestedDic()
    intergration = nestedDic()
    tag_dic = nestedDic()
    Vprint('start writing data...')
    fo_phasiRNA.write('>Hypergeometric' + "\n")
    fo_phasiRNA.write('#H cDNA based result' + "\n")
    record_marker = 0
    for gene in hout['hc']['phasiRNA']:
        feature = transcriptAnno[gene]
        annotation = feature + ":" + gene
        for cluster in hout['hc']['phasiRNA'][gene]:
            coor = (hout['hc']['phasiRNA_cluster'][gene][cluster][0], hout['hc']['phasiRNA_cluster'][gene][cluster][1])
            record_marker += 1
            out_marker = "HC" + str(record_marker)
            for pos in hout['hc']['phasiRNA'][gene][cluster]:
                pvalue = float(hout['hc']['phasiRNA'][gene][cluster][pos].split("\t")[-1])
                intergration[feature][gene]['H']['C'][coor] = ['-', '-', pvalue, out_marker]
                fo_phasiRNA.write(hout['hc']['phasiRNA'][gene][cluster][pos] + "\t" + annotation + "\t" + out_marker + "\n")

    fo_phasiRNA.write('#H gDNA based result' + "\n")
    record_marker = 0
    for gene in hout['hg']['phasiRNA']:
        for cluster in hout['hg']['phasiRNA'][gene]:
            annotation = []
            coor = (hout['hg']['phasiRNA_cluster'][gene][cluster][0], hout['hg']['phasiRNA_cluster'][gene][cluster][1])
            record_marker += 1
            out_marker = "HG" + str(record_marker)
            for pos in hout['hg']['phasiRNA'][gene][cluster]:
                anno = hanno[gene][coor]
                pvalue = float(hout['hg']['phasiRNA'][gene][cluster][pos].split("\t")[-1])
                chr_ = hout['hg']['phasiRNA'][gene][cluster][pos].split("\t")[0]                       
                for feature in hanno[gene][coor]:
                    for i in hanno[gene][coor][feature]:
                        i1 = i.split("-")[0]
                        if i1 == 'Intergenic':
                            intergration[feature][chr_ + '_' + str(coor[0]) + ':' + str(coor[1])]['H']['G'][coor] = ['-', '-', pvalue, out_marker]
                            tag_dic[f'{chr_}\t{out_marker[:2]}\t{str(coor[0])}\t{str(coor[1])}'] = out_marker
                        else:
                            intergration[feature][i]['H']['G'][coor] = ['-', '-', pvalue, out_marker]
                        annotation.append(feature + ':' + i)
                annotation1 = ";".join(list(set(annotation)))
                fo_phasiRNA.write(hout['hg']['phasiRNA'][gene][cluster][pos] + "\t" + annotation1 + "\t" + out_marker + "\n")

    fo_phasiRNA.write('#H FLNC based result' + "\n")
    record_marker = 0
        # try:
        #     feature = transcriptAnno[flnc_anno_dic[gene]]
        # except TypeError:
        #     gene = gene.replace('+', '%2B')
        #     feature = transcriptAnno[flnc_anno_dic[gene]]
    for gene in hout['hf']['phasiRNA']:
        feature = transcriptAnno[flnc_anno_dic[gene]]
        try:
            annotation = feature + ":" + flnc_anno_dic[gene]
        except TypeError:
            annotation = 'FLNC' + ":" + flnc_anno_dic[gene]
        for cluster in hout['hf']['phasiRNA'][gene]:
            coor = (hout['hf']['phasiRNA_cluster'][gene][cluster][0], hout['hf']['phasiRNA_cluster'][gene][cluster][1])
            record_marker += 1
            out_marker = "HF" + str(record_marker)
            for pos in hout['hf']['phasiRNA'][gene][cluster]:
                pvalue = float(hout['hf']['phasiRNA'][gene][cluster][pos].split("\t")[-1])
                try:
                    intergration[feature][flnc_anno_dic[gene]]['H']['F'][coor] = ['-', '-', pvalue, out_marker]
                except TypeError:
                    intergration['FLNC'][flnc_anno_dic[gene]]['H']['F'][coor] = ['-', '-', pvalue, out_marker]
                    intergration_FLNC[gene]['H']['F'][coor] = ['-', '-', pvalue, out_marker]
                fo_phasiRNA.write(hout['hf']['phasiRNA'][gene][cluster][pos] + "\t" + annotation + "\t" + out_marker + "\n")
    
    fo_phasiRNA.write('>PhaseScore' + "\n")
    fo_phasiRNA.write('#P cDNA based result' + "\n")
    record_marker = 0
    for gene in pout['pc']['phasiRNA']:
        feature = transcriptAnno[gene]
        annotation = feature + ":" + gene
        for record in pout['pc']['phasiRNA'][gene]:
            coor = (p_phasiRNA_cluster[0][gene][record][0], p_phasiRNA_cluster[0][gene][record][1])
            record_marker += 1 
            out_marker = "PC" + str(record_marker)
            for key in pout['pc']['phasiRNA'][gene][record]:
                phaseRatio = float(pout['pc']['phasiRNA'][gene][record][key].split("\t")[-6])
                phaseScore = float(pout['pc']['phasiRNA'][gene][record][key].split("\t")[-3])
                intergration[feature][gene]['P']['C'][coor] = [phaseRatio, phaseScore, '-', out_marker]
                tmp  = "\t".join(pout['pc']['phasiRNA'][gene][record][key].split("\t")[:-2]) + "\t" + pout['pc']['phasiRNA'][gene][record][key].split("\t")[-1]
                fo_phasiRNA.write(tmp + "\t" + annotation + "\t" + out_marker + "\n")

    fo_phasiRNA.write('#P gDNA based result' + "\n")
    record_marker = 0
    for gene in pout['pg']['phasiRNA']:
        for record in pout['pg']['phasiRNA'][gene]:
            annotation = []
            coor = (p_phasiRNA_cluster[1][gene][record][0], p_phasiRNA_cluster[1][gene][record][1])
            record_marker += 1
            out_marker = "PG" + str(record_marker)
            for key in pout['pg']['phasiRNA'][gene][record]:
                phaseRatio = float(pout['pg']['phasiRNA'][gene][record][key].split("\t")[-6])
                phaseScore = float(pout['pg']['phasiRNA'][gene][record][key].split("\t")[-3])
                chr_ = pout['pg']['phasiRNA'][gene][record][key].split("\t")[0]
                tmp  = "\t".join(pout['pg']['phasiRNA'][gene][record][key].split("\t")[1:-2]) + "\t" + pout['pg']['phasiRNA'][gene][record][key].split("\t")[-1]
                geneid = pout['pg']['phasiRNA'][gene][record][key].split("\t")[0].split('_')[:-1]
                geneid = "_".join(geneid)
                for feature in panno[geneid][coor]:
                    for i in panno[geneid][coor][feature]:
                        i1 = i.split('-')[0]
                        if i1 == 'Intergenic':
                            intergration[feature][geneid + "_" + str(coor[0]) + ":" + str(coor[1])]['P']['G'][coor] = [phaseRatio, phaseScore, '-', out_marker]
                            tag_dic[f'{geneid}\t{out_marker[:2]}\t{str(coor[0])}\t{str(coor[1])}'] = out_marker
                        else:
                            intergration[feature][i]['P']['G'][coor] = [phaseRatio, phaseScore, '-', out_marker]
                        annotation.append(feature + ':' + i)
                annotation1 = ';'.join(list(set(annotation)))
                fo_phasiRNA.write(geneid + "\t" + tmp + "\t" + annotation1 + "\t" + out_marker + "\n")
    
    fo_phasiRNA.write('#P FLNC based result' + "\n")
    record_marker = 0
    for gene in pout['pf']['phasiRNA']:
        feature = transcriptAnno[flnc_anno_dic[gene]]
        try:
            annotation = feature + ":" + flnc_anno_dic[gene]
        except TypeError:
            annotation = 'FLNC' + ":" + flnc_anno_dic[gene]
        for record in pout['pf']['phasiRNA'][gene]:
            coor = (p_phasiRNA_cluster[2][gene][record][0], p_phasiRNA_cluster[2][gene][record][1])
            record_marker += 1
            out_marker = "PF" + str(record_marker)
            for key in pout['pf']['phasiRNA'][gene][record]:
                phaseRatio = float(pout['pf']['phasiRNA'][gene][record][key].split("\t")[-6])
                phaseScore = float(pout['pf']['phasiRNA'][gene][record][key].split("\t")[-3])
                try:
                    intergration[feature][flnc_anno_dic[gene]]['P']['F'][coor] = [phaseRatio, phaseScore, '-', out_marker]
                except TypeError:
                    intergration['FLNC'][flnc_anno_dic[gene]]['P']['F'][coor] = [phaseRatio, phaseScore, '-', out_marker]
                    intergration_FLNC[gene]['P']['F'][coor] = [phaseRatio, phaseScore, '-', out_marker]
                tmp  = "\t".join(pout['pf']['phasiRNA'][gene][record][key].split("\t")[:-2]) + "\t" + pout['pf']['phasiRNA'][gene][record][key].split("\t")[-1]
                fo_phasiRNA.write(tmp + "\t" + annotation + "\t" + out_marker + "\n")

    fo_allsiRNA.write('>Hypergeometric' + "\n")
    fo_allsiRNA.write('#H cDNA based result' + "\n")
    record_marker = 0
    for gene in hout['hc']['allsiRNA']:
        annotation = transcriptAnno[gene] + ':' + gene
        for cluster in hout['hc']['allsiRNA'][gene]:
            record_marker += 1
            out_marker = "HC" + str(record_marker)
            for pos in hout['hc']['allsiRNA'][gene][cluster]:
                # for i1 in hout['hc']['allsiRNA'][gene][cluster][pos]:
                i1 = hout['hc']['allsiRNA'][gene][cluster][pos]
                fo_allsiRNA.write(i1 + "\t" + annotation + "\t" + out_marker + "\n")

    fo_allsiRNA.write('#H gDNA based result' + "\n")
    record_marker = 0
    for gene in hout['hg']['allsiRNA']:
        for cluster in hout['hg']['allsiRNA'][gene]:
            annotation = []
            coor = (hout['hg']['phasiRNA_cluster'][gene][cluster][0], hout['hg']['phasiRNA_cluster'][gene][cluster][1])
            record_marker += 1
            out_marker = "HG" + str(record_marker)
            for pos in hout['hg']['allsiRNA'][gene][cluster]:
                # for i1 in hout['hg']['allsiRNA'][gene][cluster][pos]:
                i1 = hout['hg']['allsiRNA'][gene][cluster][pos]
                for feature in hanno[gene][coor]:
                    for i in hanno[gene][coor][feature]: 
                        annotation.append(feature + ':' + i)
                annotation1 = ';'.join(list(set(annotation)))
                fo_allsiRNA.write(i1 + "\t" + annotation1 + "\t" + out_marker + "\n")
    
    fo_allsiRNA.write('#H FLNC based result' + "\n")
    record_marker = 0
    for gene in hout['hf']['allsiRNA']:
        feature = transcriptAnno[flnc_anno_dic[gene]]
        try:
            annotation = feature + ":" + flnc_anno_dic[gene]
        except TypeError:
            annotation = 'FLNC' + ":" + flnc_anno_dic[gene]
        for cluster in hout['hf']['allsiRNA'][gene]:
            record_marker += 1
            out_marker = "HF" + str(record_marker)
            for pos in hout['hf']['allsiRNA'][gene][cluster]:
                i1 = hout['hf']['allsiRNA'][gene][cluster][pos]
                fo_allsiRNA.write(i1 + "\t" + annotation + "\t" + out_marker + "\n")

    fo_allsiRNA.write('>PhaseScore' + "\n")
    fo_allsiRNA.write('#P cDNA based result' + "\n")
    record_marker = 0
    for gene in pout['pc']['allsiRNA']:
        annotation = transcriptAnno[gene] + ':' + gene
        for record in pout['pc']['allsiRNA'][gene]:
            record_marker += 1
            out_marker = "PC" + str(record_marker)
            for key in pout['pc']['allsiRNA'][gene][record]:
                tmp  = "\t".join(pout['pc']['allsiRNA'][gene][record][key].split("\t")[:-2]) + "\t" + pout['pc']['allsiRNA'][gene][record][key].split("\t")[-1]
                fo_allsiRNA.write(tmp + "\t" + annotation + "\t" + out_marker + "\n")

    fo_allsiRNA.write('#P gDNA based result' + "\n")
    record_marker = 0
    for gene in pout['pg']['allsiRNA']:
        for record in pout['pg']['allsiRNA'][gene]:
            annotation = []
            coor = (p_phasiRNA_cluster[1][gene][record][0], p_phasiRNA_cluster[1][gene][record][1])
            record_marker += 1
            out_marker = "PG" + str(record_marker)
            for key in pout['pg']['allsiRNA'][gene][record]:
                tmp  = "\t".join(pout['pg']['allsiRNA'][gene][record][key].split("\t")[1:-2]) + "\t" + pout['pg']['allsiRNA'][gene][record][key].split("\t")[-1]
                geneid = pout['pg']['allsiRNA'][gene][record][key].split("\t")[0].split('_')[:-1]
                geneid = "_".join(geneid)
                for feature in panno[geneid][coor]:
                    for i in panno[geneid][coor][feature]:
                        annotation.append(feature + ":" + i)
                annotation1 = ';'.join(list(set(annotation)))
                fo_allsiRNA.write(geneid + "\t" + tmp + "\t" + annotation1 + "\t" + out_marker + "\n")
    
    fo_allsiRNA.write('#P FLNC based result' + "\n")
    record_marker = 0
    for gene in pout['pf']['allsiRNA']:
        feature = transcriptAnno[flnc_anno_dic[gene]]
        try:
            annotation = feature + ":" + flnc_anno_dic[gene]
        except TypeError:
            annotation = 'FLNC' + ":" + flnc_anno_dic[gene]
        for record in pout['pf']['allsiRNA'][gene]:
            record_marker += 1
            out_marker = "PF" + str(record_marker)
            for key in pout['pf']['allsiRNA'][gene][record]:
                tmp  = "\t".join(pout['pf']['allsiRNA'][gene][record][key].split("\t")[:-2]) + "\t" + pout['pf']['allsiRNA'][gene][record][key].split("\t")[-1]
                fo_allsiRNA.write(tmp + "\t" + annotation + "\t" + out_marker + "\n")
    
    return (intergration, intergration_FLNC, tag_dic)

def ParallelHypergeometric(allsiRNA, phasiRNA, setting1, parallel_number, island_number):
    # worked
    # hgdna_allsiRNA = allsiRNA['hg']
    # tmp = SecondScaning(hgdna_allsiRNA, setting1)
    # hout = nestedDic()
    # hout['hg']['phasiRNA'] = tmp[0]
    # hout['hg']['allsiRNA'] = tmp[1]
    # hout['hg']['phasiRNA_cluster'] = tmp[2]
    # hout['hg']['phasiRNA_strand'] = tmp[3]
    # h = ['hc', 'hf']
    # for i in allsiRNA:
    #     if i in h:
    #         setting3 = (phasiRNA[0][i], phasiRNA[1][i], setting2[0], setting2[1], setting2[2], setting2[3])
    #         tmp = SecondScaning(allsiRNA[i], setting3)
    #         hout[i]['phasiRNA'] = tmp[0]
    #         hout[i]['allsiRNA'] = tmp[1]
    #         hout[i]['phasiRNA_cluster'] = tmp[2]
    #         hout[i]['phasiRNA_strand'] = tmp[3]
    # return hout

    hgdna_allsiRNA = allsiRNA['hg']
    hcdna_allsiRNA = allsiRNA['hc']
    hfdna_allsiRNA = allsiRNA['hf']
    hcdna_phasiRNA_cluster = phasiRNA[0]['hc']
    hcdna_phasiRNA_strand = phasiRNA[1]['hc']
    hfdna_phasiRNA_cluster = phasiRNA[0]['hf']
    hfdna_phasiRNA_strand = phasiRNA[1]['hf']
    process_pool = ProcessPoolExecutor(max_workers=parallel_number)
    futures_hg = [] 
    futures_hc = []
    futures_hf = []
    for i in hgdna_allsiRNA:
        Input = (hgdna_allsiRNA[i], setting1[0][i], setting1[1][i], setting1[2], setting1[3], setting1[4], setting1[5], i, island_number)
        future_hg = process_pool.submit(ParallelSecondScaning, Input)
        futures_hg.append(future_hg)

    for i in hcdna_allsiRNA:
        Input = (hcdna_allsiRNA[i], hcdna_phasiRNA_cluster[i], hcdna_phasiRNA_strand[i], setting1[2], setting1[3], setting1[4], setting1[5], i, island_number)
        future_hc = process_pool.submit(ParallelSecondScaning, Input)
        futures_hc.append(future_hc)

    for i in hfdna_allsiRNA:
        Input = (hfdna_allsiRNA[i], hfdna_phasiRNA_cluster[i], hfdna_phasiRNA_strand[i], setting1[2], setting1[3], setting1[4], setting1[5], i, island_number)
        future_hf = process_pool.submit(ParallelSecondScaning, Input)
        futures_hf.append(future_hf)
    # for i in hfdna_allsiRNA:
    #     Input = (hfdna_allsiRNA[i], hfdna_phasiRNA_cluster[i], hfdna_phasiRNA_strand[i], setting1[2], setting1[3], setting1[4], setting1[5], i)
    #     future_hf = process_pool.submit(ParallelHypergeometric, Input)
    #     futures_hf.append(futures_hf)

    wait(futures_hg, return_when=ALL_COMPLETED)
    wait(futures_hc, return_when=ALL_COMPLETED)
    wait(futures_hf, return_when=ALL_COMPLETED)

    process_pool.shutdown()
    hout = nestedDic()
    for i in futures_hg:
        future_hg = i.result()
        for gene in future_hg[0]:
            hout['hg']['phasiRNA'][gene] = future_hg[0][gene]
        for gene in future_hg[1]:
            hout['hg']['allsiRNA'][gene] = future_hg[1][gene]
        for gene in future_hg[2]:
            hout['hg']['phasiRNA_cluster'][gene] = future_hg[2][gene]
        for gene in future_hg[3]:
            hout['hg']['allsiRNA_cluster'][gene] = future_hg[3][gene]
    for i in futures_hc:
        future_hc = i.result()
        for gene in future_hc[0]:
            hout['hc']['phasiRNA'][gene] = future_hc[0][gene]
        for gene in future_hc[1]:
            hout['hc']['allsiRNA'][gene] = future_hc[1][gene]
        for gene in future_hc[2]:
            hout['hc']['phasiRNA_cluster'][gene] = future_hc[2][gene]
        for gene in future_hc[3]:
            hout['hc']['allsiRNA_cluster'][gene] = future_hc[3][gene]
    for i in futures_hf:
        future_hf = i.result()
        for gene in future_hf[0]:
            hout['hf']['phasiRNA'][gene] = future_hf[0][gene]
        for gene in future_hf[1]:
            hout['hf']['allsiRNA'][gene] = future_hf[1][gene]
        for gene in future_hf[2]:
            hout['hf']['phasiRNA_cluster'][gene] = future_hf[2][gene]
        for gene in future_hf[3]:
            hout['hf']['allsiRNA_cluster'][gene] = future_hf[3][gene]

    # reorder cluster id
    new_hout = nestedDic()
    for method in hout:
        for rna in hout[method]:
            count = 0
            for gene in hout[method][rna]:
                for recordid in hout[method][rna][gene]:
                    count += 1
                    new_hout[method][rna][gene][count] = hout[method][rna][gene][recordid]

    return new_hout


def Hypergeometric(hgdna_allsiRNA, allsiRNA, phasiRNA, setting1, setting2):
    tmp = SecondScaning(hgdna_allsiRNA, setting1)
    hout = nestedDic()
    hout['hg']['phasiRNA'] = tmp[0]
    hout['hg']['allsiRNA'] = tmp[1]
    hout['hg']['phasiRNA_cluster'] = tmp[2]
    hout['hg']['phasiRNA_strand'] = tmp[3]
    h = ['hc', 'hf']
    for i in allsiRNA:
        if i in h:
            setting3 = (phasiRNA[0][i], phasiRNA[1][i], setting2[0], setting2[1], setting2[2], setting2[3])
            tmp = SecondScaning(allsiRNA[i], setting3)
            hout[i]['phasiRNA'] = tmp[0]
            hout[i]['allsiRNA'] = tmp[1]
            hout[i]['phasiRNA_cluster'] = tmp[2]
            hout[i]['phasiRNA_strand'] = tmp[3]
    return hout

def ParallelSecondScaning(parameter):
    """second scaning

    Parameters
    ----------
    final_allsiRNA : dict
        first scaning output dict
    setting : tuple
        setting = (final_cluster, final_strand, phase_number, phase_length, pvalue_cutoff, min_read_num)

    Returns
    -------
    tuple
        (out_phasiRNA, out_allsiRNA, out_phasiRNA_cluster, out_allsiRNA_cluster)
    """
    final_allsiRNA = parameter[0]
    setting = parameter[1:]
    final_gene = {}
    duplication_phasiRNA = nestedDic()
    duplication_allsiRNA = nestedDic()
    record_marker = 0
    final_cluster = setting[0]
    final_strand = setting[1]
    phase_number = setting[2]
    phase_length = setting[3]
    pvalue_cutoff = setting[4]
    min_read_num = setting[5]
    geneid = setting[6]
    island_number = setting[7]
    window_length = phase_length*11
    out_allsiRNA = nestedDic()
    out_phasiRNA = nestedDic()
    out_phasiRNA_cluster = nestedDic()
    out_allsiRNA_cluster = nestedDic()
    out_duplication_phasiRNA_cluster = nestedDic()
    out_duplication_allsiRNA_cluster = nestedDic()

    max_readn = 0
    sorted_list = sorted(list(set(final_cluster)))
    candidate_cluster0 = SplitClusterIsland(sorted_list, island_number, phase_length)
    candidate_cluster = SplitIsland1(candidate_cluster0)
    # candidate_cluster = SplitCluster(sorted_list, island_number, phase_length)
    for j in candidate_cluster:
        allsiRNA = []
        phasiRNA = []
        tag = False
        tag1 = 'n'
        pair_pos = PairPos(list(range(len(candidate_cluster[j]))), candidate_cluster[j])
        for e in list(sorted(pair_pos.keys()))[::-1]:
            if e >= 4:
                for length in list(sorted(pair_pos[e].keys()))[::-1]:
                    for tuple_ in pair_pos[e][length]:
                        # allsiRNA = []
                        # phasiRNA = []
                        total_n = 0
                        total_k = 0
                        pvalue = 0
                        if tuple_[0] > tuple_[1]:
                            start_pos = candidate_cluster[j][tuple_[1]]
                            tmp_pos = candidate_cluster[j][tuple_[0]]
                        else:
                            start_pos = candidate_cluster[j][tuple_[0]]
                            tmp_pos = candidate_cluster[j][tuple_[1]]
                        remainder = (tmp_pos - start_pos)%phase_length
                        if remainder == 0:
                            end_pos = tmp_pos
                        elif remainder == 2:
                            end_pos = tmp_pos + 19
                        elif remainder == 19 and phase_length == 21:
                            end_pos = tmp_pos + 2
                        elif remainder == 22 and phase_length == 24:
                            end_pos = tmp_pos + 2
                        else:
                            end_pos = tmp_pos
                        strand = final_strand[start_pos]
                        end_pos = end_pos + phase_length
                        if strand == '+':
                            for i in range(start_pos, end_pos):
                                # query = geneid + '\t' + '+' + str(i)
                                query = geneid + '\t' + '+' + '\t' + str(i)
                                if query in final_allsiRNA:
                                    total_n += 1
                                    allsiRNA.append(str(final_allsiRNA[query]))
                            for i in range(start_pos - 2, end_pos - 2):
                                # query = geneid + '\t' + '-' + str(i)
                                query = geneid + '\t' + '-' + '\t' + str(i)
                                if query in final_allsiRNA:
                                    total_n += 1
                                    allsiRNA.append(str(final_allsiRNA[query]))
                            for i in range(start_pos, end_pos, phase_length):
                                query = geneid + '\t' + '+' + '\t' + str(i)
                                if query in final_allsiRNA and final_allsiRNA[query].split('\t')[6] == str(phase_length):
                                    total_k += 1
                                    phasiRNA.append(final_allsiRNA[query])
                            for i in range(start_pos - 2, end_pos - 2, phase_length):
                                query = geneid + '\t' + '-' + '\t' + str(i)
                                if query in final_allsiRNA and final_allsiRNA[query].split('\t')[6] == str(phase_length):
                                    total_k += 1
                                    phasiRNA.append(final_allsiRNA[query])
                        elif strand == '-':
                            for i in range(start_pos, end_pos):
                                query = geneid + '\t' + '-' + '\t' + str(i)
                                if query in final_allsiRNA:
                                    total_n += 1
                                    allsiRNA.append(str(final_allsiRNA[query]))
                            for i in range(start_pos + 2, end_pos + 2):
                                query = geneid + '\t' + '+' + '\t' + str(i)
                                if query in final_allsiRNA:
                                    total_n += 1
                                    allsiRNA.append(str(final_allsiRNA[query]))
                            for i in range(start_pos, end_pos, phase_length):
                                query = geneid + '\t' + '-' + '\t' + str(i)
                                if query in final_allsiRNA and final_allsiRNA[query].split('\t')[6] == str(phase_length):
                                    total_k += 1
                                    phasiRNA.append(str(final_allsiRNA[query]))
                            for i in range(start_pos + 2, end_pos + 2, phase_length):
                                query = geneid + '\t' + '+' + '\t' + str(i)
                                if query in final_allsiRNA and final_allsiRNA[query].split('\t')[6] == str(phase_length):
                                    total_k += 1
                                    phasiRNA.append(str(final_allsiRNA[query]))

                        # phase_number and pvalue cutoff
                        if total_k >= phase_number and (end_pos - start_pos) >= 100:
                        # if total_k >= phase_number:
                            # for i in range(total_k, phase_length):
                            #     pvalue += float(binomial((end_pos - start_pos) * 2 - 1 - phase_length, total_n - i) * binomial(phase_length, i) / binomial((end_pos - start_pos) * 2 - 1, total_n))
                            tmp_variable = int((((end_pos - start_pos) * 2) / phase_length) - 1)
                            for i in range(total_k, tmp_variable):
                                # pvalue += float(binomial((end_pos - start_pos) * 2 - 1 - phase_length, total_n - i) * binomial(phase_length, i) / binomial((end_pos - start_pos) * 2 - 1, total_n))
                                pvalue += float(binomial((end_pos - start_pos) * 2 - 1 - tmp_variable, total_n - i) * binomial(tmp_variable, i) / binomial((end_pos - start_pos) * 2 - 1, total_n))
                            if pvalue <= pvalue_cutoff and pvalue != 0.0:
                                for i in range(len(phasiRNA)):
                                    value_list = phasiRNA[i].split('\t')
                                    sRNA_no = float(value_list[3])
                                    if sRNA_no > max_readn:
                                        max_readn = sRNA_no

                                if min_read_num == 0 or max_readn >= min_read_num:  # min_read_num cutoff
                                    record_marker += 1
                                    tag = True
                                    tag1 = 'y'
                                    for i in range(len(phasiRNA)):
                                        value_list = phasiRNA[i].split('\t')
                                        geneid = value_list[0]
                                        strand = value_list[1]
                                        readpos = int(value_list[2])
                                        if geneid in final_gene:
                                            if pvalue < float(final_gene[geneid].split('\t')[1]):
                                                final_gene[geneid] = str(max_readn) + '\t' + str(pvalue)
                                            out_phasiRNA[geneid][record_marker][readpos] = "\t".join(phasiRNA[i].split("\t")[:-2]) + "\t" + str(pvalue)
                                        else:
                                            final_gene[geneid] = str(max_readn) + '\t' + str(pvalue)
                                            out_phasiRNA[geneid][record_marker][readpos] = "\t".join(phasiRNA[i].split("\t")[:-2]) + "\t" + str(pvalue)
                                    start_pos = list(sorted(out_phasiRNA[geneid][record_marker].keys()))[0]
                                    end_pos = list(sorted(out_phasiRNA[geneid][record_marker].keys()))[-1] + phase_length
                                    out_phasiRNA_cluster[geneid][record_marker] = [start_pos, end_pos]

                                    for i in range(len(allsiRNA)):
                                        value_list = allsiRNA[i].split('\t')
                                        geneid = value_list[0]
                                        readpos = int(value_list[2])
                                        query = geneid + '\t' + '+' + '\t' + str(i)
                                        out_allsiRNA[geneid][record_marker][readpos] = "\t".join(allsiRNA[i].split("\t")[:-2]) + "\t" + str(pvalue)
                                    start_pos = list(sorted(out_allsiRNA[geneid][record_marker].keys()))[0]
                                    end_pos = list(sorted(out_allsiRNA[geneid][record_marker].keys()))[-1] + phase_length
                                    out_allsiRNA_cluster[geneid][record_marker] = [start_pos, end_pos]
                        if tag:
                            break
                    if tag:
                        break
            if tag:
                break


        if tag1 == 'n':
            for e in list(sorted(pair_pos.keys()))[::-1]:
                if e >= 4:
                    for length in list(sorted(pair_pos[e].keys())):
                        for tuple_ in pair_pos[e][length]:
                            total_n = 0
                            total_k = 0
                            pvalue = 0
                            if tuple_[0] > tuple_[1]:
                                start_pos = candidate_cluster[j][tuple_[1]]
                                tmp_pos = candidate_cluster[j][tuple_[0]]
                            else:
                                start_pos = candidate_cluster[j][tuple_[0]]
                                tmp_pos = candidate_cluster[j][tuple_[1]]
                            remainder = (tmp_pos - start_pos)%phase_length
                            if remainder == 0:
                                end_pos = tmp_pos
                            elif remainder == 2:
                                end_pos = tmp_pos + 19
                            elif remainder == 19 and phase_length == 21:
                                end_pos = tmp_pos + 2
                            elif remainder == 22 and phase_length == 24:
                                end_pos = tmp_pos + 2
                            else:
                                end_pos = tmp_pos
                            strand = final_strand[start_pos]
                            if strand == '+':
                                for i in range(start_pos, window_length + start_pos):
                                    query = geneid + '\t' + '+' + '\t' + str(i)
                                    if query in final_allsiRNA:
                                        total_n += 1
                                        allsiRNA.append(str(final_allsiRNA[query]))
                                for i in range(start_pos - 2, start_pos + window_length - 2):
                                    query = geneid + '\t' + '-' + '\t' + str(i)
                                    if query in final_allsiRNA:
                                        total_n += 1
                                        allsiRNA.append(str(final_allsiRNA[query]))
                                for i in range(start_pos, start_pos + window_length - 2 , phase_length):
                                    query = geneid + '\t' + '+' + '\t' + str(i)
                                    if query in final_allsiRNA and final_allsiRNA[query].split('\t')[6] == str(phase_length):
                                        total_k += 1
                                        phasiRNA.append(final_allsiRNA[query])
                                for i in range(start_pos - 2, start_pos + window_length - 2, phase_length):
                                    query = geneid + '\t' + '-' + '\t' + str(i)
                                    if query in final_allsiRNA and final_allsiRNA[query].split('\t')[6] == str(phase_length):
                                        total_k += 1
                                        phasiRNA.append(final_allsiRNA[query])
                            elif strand == '-':
                                for i in range(start_pos, start_pos + window_length):
                                    query = geneid + '\t' + '-' + '\t' + str(i)
                                    if query in final_allsiRNA:
                                        total_n += 1
                                        allsiRNA.append(str(final_allsiRNA[query]))
                                for i in range(start_pos + 2, start_pos + window_length + 2):
                                    query = geneid + '\t' + '+' + '\t' + str(i)
                                    if query in final_allsiRNA:
                                        total_n += 1
                                        allsiRNA.append(str(final_allsiRNA[query]))
                                for i in range(start_pos, start_pos + window_length, phase_length):
                                    query = geneid + '\t' + '-' + '\t' + str(i)
                                    if query in final_allsiRNA and final_allsiRNA[query].split('\t')[6] == str(phase_length):
                                        total_k += 1
                                        phasiRNA.append(str(final_allsiRNA[query]))
                                for i in range(start_pos + 2, start_pos + window_length + 2, phase_length):
                                    query = geneid + '\t' + '+' + '\t' + str(i)
                                    if query in final_allsiRNA and final_allsiRNA[query].split('\t')[6] == str(phase_length):
                                        total_k += 1
                                        phasiRNA.append(str(final_allsiRNA[query]))

                            # phase_number and pvalue cutoff
                            if total_k >= phase_number:
                                # for i in range(total_k, phase_length):
                                    # pvalue += float(binomial((window_length) * 2 - 1 - phase_length, total_n - i) * binomial(phase_length, i) / binomial((window_length) * 2 - 1, total_n))
                                tmp_variable = int(((window_length * 2) / phase_length) - 1)
                                for i in range(total_k, tmp_variable):
                                    pvalue += float(binomial((window_length) * 2 - 1 - tmp_variable, total_n - i) * binomial(tmp_variable, i) / binomial((window_length) * 2 - 1, total_n))
                                if pvalue <= pvalue_cutoff and pvalue != 0.0:
                                    for i in range(len(phasiRNA)):
                                        value_list = phasiRNA[i].split('\t')
                                        sRNA_no = float(value_list[3])
                                        if sRNA_no > max_readn:
                                            max_readn = sRNA_no

                                    if min_read_num == 0 or max_readn >= min_read_num:  # min_read_num cutoff
                                        record_marker += 1
                                        tag = True
                                        for i in range(len(phasiRNA)):
                                            value_list = phasiRNA[i].split('\t')
                                            geneid = value_list[0]
                                            strand = value_list[1]
                                            readpos = int(value_list[2])
                                            if geneid in final_gene:
                                                if pvalue < float(final_gene[geneid].split('\t')[1]):
                                                    final_gene[geneid] = str(max_readn) + '\t' + str(pvalue)
                                                out_phasiRNA[geneid][record_marker][readpos] = "\t".join(phasiRNA[i].split("\t")[:-2]) + "\t" + str(pvalue)
                                            else:
                                                final_gene[geneid] = str(max_readn) + '\t' + str(pvalue)
                                                out_phasiRNA[geneid][record_marker][readpos] = "\t".join(phasiRNA[i].split("\t")[:-2]) + "\t" + str(pvalue)
                                        start_pos = list(sorted(out_phasiRNA[geneid][record_marker].keys()))[0]
                                        end_pos = sorted(out_phasiRNA[geneid][record_marker].keys())[-1] + phase_length
                                        out_phasiRNA_cluster[geneid][record_marker] = [start_pos, end_pos]
                                        for i in range(len(allsiRNA)):
                                            value_list = allsiRNA[i].split('\t')
                                            geneid = value_list[0]
                                            readpos = int(value_list[2])
                                            query = geneid + '\t' + '+' + '\t' + str(i)
                                            out_allsiRNA[geneid][record_marker][readpos] = "\t".join(allsiRNA[i].split("\t")[:-2]) + "\t" + str(pvalue)
                                        start_pos = list(sorted(out_allsiRNA[geneid][record_marker].keys()))[0]
                                        out_allsiRNA_cluster[geneid][record_marker] = [start_pos, start_pos + window_length]
                            if tag:
                                break
                        if tag:
                            break
                if tag:
                    break

    return (out_phasiRNA, out_allsiRNA, out_phasiRNA_cluster, out_allsiRNA_cluster)


def loadFLNCanno(flnc_anno):
    dic = nestedDic()
    with open(flnc_anno, 'r') as fn:
        for line in fn:
            l = line.strip().split("\t")
            flnc_id = l[0]
            isotype = l[1]
            gene = l[2]
            rna = l[3]
            flnc_id = flnc_id.replace('%2B', '+')
            if rna != 'novel' and rna != 'Novel':
                dic[flnc_id] = rna
            else:
                dic[flnc_id] = isotype+'_'+gene+"_"+rna
    return dic

def overlapIntergenic(tmp_fo_intergration, tmp_Intergenic, tmp_Intergenic1, tmp_Intergenic_, tmp_Intergenic2):
    # cmd = f'tail -n +2 {tmp_fo_intergration} | grep -v Intergenic >> {Non_intergenic}'
    # print(cmd)
    # os.system(cmd)
    cmd = f'tail -n +2 {tmp_fo_intergration} | grep Intergenic > {tmp_Intergenic_}'
    os.system(cmd)
    cmd = f'tail -n +2 {tmp_fo_intergration} | grep Intergenic | cut -f 3 | phasiHunter read 1 > {tmp_Intergenic}'
    os.system(cmd)
    cmd = f'paste {tmp_Intergenic} {tmp_Intergenic_} > {tmp_Intergenic1}'
    os.system(cmd)
    cmd = f'sort -k1,1 -k2,2n -k3,3n {tmp_Intergenic1} > {tmp_Intergenic2}'
    # print(cmd)
    os.system(cmd)
    cmd = f'bedtools merge -i {tmp_Intergenic2} -c 4,7,10,13,16,19,22 -o max,collapse,collapse,collapse,collapse,collapse,collapse'
    # os.system(cmd)
    overlap_Intergration = os.popen(cmd)
    return overlap_Intergration


def Integration_overlap(overlap_Intergenic_fo, fo_intergration_tmp, passP):
    fo_intergration_tmp.write(f'cluster_number\tfeature\tPHAS_gene\tmethod\treference_type\th_transcriptome_coordinate\th_genome_coordinate\th_fl_transcriptome_coordinate\tp_transcriptome_coordinate\tp_genome_coordinate\tp_fl_transcriptome_coordinate\tc_pvalue\tg_pvalue\tf_pvalue\tc_phaseratio\tg_phaseratio\tf_phaseratio\tc_phasescore\tg_phasescore\tf_phasescore\n')
    for i in overlap_Intergenic_fo:
        l = i.strip().split("\t")
        geneid_, start, end, cluster_number, method, h_genome_coordinate, p_genome_coordinate, g_pvalue, g_phaseratio, g_phasescore = \
        l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8], l[9]
        if passP == 'y' and 'H' not in method:
            continue
        geneid_new = f'{geneid_}_{start}:{end}'
        feature = 'Intergenic'
        reference_type = 'G'
        h_transcriptome_coordinate = '-'
        h_fl_transcriptome_coordinate = '-'
        p_transcriptome_coordinate = '-'
        p_fl_transcriptome_coordinate = '-'
        c_pvalue = '-'
        f_pvalue = '-'
        c_phaseratio = '-'
        f_phaseratio = '-'
        c_phasescore = '-'
        f_phasescore = '-'
        realign = [cluster_number,feature,geneid_new,method,reference_type,h_transcriptome_coordinate,h_genome_coordinate,h_fl_transcriptome_coordinate,p_transcriptome_coordinate,p_genome_coordinate,p_fl_transcriptome_coordinate,c_pvalue,g_pvalue,f_pvalue,c_phaseratio,g_phaseratio,f_phaseratio,c_phasescore,g_phasescore,f_phasescore]
        fo_intergration_tmp.write("\t".join(realign) + "\n")

def catCombine(Non_intergenic, intergrationfile):
    cmd = f'cat {Non_intergenic} >> {intergrationfile}'
    # print(cmd)
    os.system(cmd)

def WriteIntergration_new_dup(intergration, tmp_gff_bed, relation_dic, intergration_FLNC, flnc_anno_dic, gff3type):
    trans = nestedDic()
    with open(tmp_gff_bed, 'r') as fn:
        for line in fn:
            line1 = line.strip()
            l = line.strip().split("\t")
            chr_,start,end,feature,anno,strand = l[0],int(l[1]),int(l[2]),l[3],l[4],l[5]
            if feature == 'exon':
                if gff3type == 'ensemble':
                    gene = re.sub('-E\d+$', '', anno)
                elif gff3type == 'ncbi':
                    gene = anno.split('-')[0]
                elif gff3type == 'phytozome':
                    gene = anno.split('-')[0]
                if gene in intergration:
                    recoder = anno.split('-')[1]
                    if strand == '-':
                        if recoder == '1' or recoder == '01':
                            tmp = 0
                            tmp = end - start + 1 + tmp
                            for i in range(1, end-start+2):
                                trans[gene][i] = end + 1 - i
                        else:
                            for i in range(tmp+1, end - start + 2 + tmp):
                                trans[gene][i] = end + tmp + 1 - i
                            tmp = end - start + 1 + tmp
                    elif strand == '+':
                        if recoder == '1' or recoder == '01':
                            tmp = 0
                            tmp = end - start + 1 + tmp
                            for i in range(1, end-start+2):
                                trans[gene][i] = start - 1 + i  
                        else:
                            for i in range(tmp+1, end -start + 2 + tmp):
                                trans[gene][i] = start - 1 + i - tmp 
                            tmp = end - start + 1 + tmp

    # consider only P method?
    list_ = []
    # list_.append(f'feature\tPHAS_Loci\tGenome_start\tGenome_end\ttranscript_start\ttranscript_end\tmethod_ref\tpvalue\tphase_score\tphase_ratio\ttag')
    for gene in list(intergration.keys()):
        if relation_dic[gene] == 'Other' or relation_dic[gene] == 'intron' or relation_dic[gene] == 'FLNC':
            continue
        for ref in ['C', 'G', 'F']:
            for method in ['H', 'P']:
                for coor in list(intergration[gene][method][ref].keys()):
                    genome_start, genome_end, transcript_start, transcript_end, pvalue, phase_score, phase_ratio, tag = \
                    '-', '-', '-', '-', '-', '-', '-', '-'
                    c_genome_start, c_genome_end, c_transcript_start, c_transcript_end, c_pvalue, c_phase_score, c_phase_ratio, c_tag = \
                    '-', '-', '-', '-', '-', '-', '-', '-'
                    c1_genome_start, c1_genome_end, c1_transcript_start, c1_transcript_end, c1_pvalue, c1_phase_score, c1_phase_ratio, c1_tag = \
                    '-', '-', '-', '-', '-', '-', '-', '-'
                    method_ref = []
                    method_ref.append(method+ref)
                    if intergration[gene][method][ref][coor] == '-':
                        continue
                    if ref == 'C' and method == 'H':
                        c_transcript_start = coor[0]
                        c_transcript_end = coor[1]
                        c_pvalue = intergration[gene][method][ref][coor][2]
                        c_genome_start = trans[gene][coor[0]]
                        c_genome_end = trans[gene][coor[1]]
                        c_tag = intergration[gene][method][ref][coor][3]
                    if ref == 'C' and method == 'P':
                        c_transcript_start = coor[0]
                        c_transcript_end = coor[1]
                        c_phase_score = intergration[gene][method][ref][coor][1]
                        c_phase_ratio = intergration[gene][method][ref][coor][0]
                        c_genome_start = trans[gene][coor[0]]
                        c_genome_end = trans[gene][coor[1]]
                        c_tag = intergration[gene][method][ref][coor][3]
                    if ref == 'F' and method == 'H':
                        c_transcript_start = coor[0]
                        c_transcript_end = coor[1]
                        c_pvalue = intergration[gene][method][ref][coor][2]
                        c_tag = intergration[gene][method][ref][coor][3]
                    if ref == 'F' and method == 'P':
                        c_transcript_start = coor[0]
                        c_transcript_end = coor[1]
                        c_phase_score = intergration[gene][method][ref][coor][1]
                        c_phase_ratio = intergration[gene][method][ref][coor][0]
                        c_tag = intergration[gene][method][ref][coor][3]
                    if ref == 'G' and method == 'H':
                        c_genome_start = coor[0]
                        c_genome_end = coor[1]
                        c_pvalue = intergration[gene][method][ref][coor][2]
                        c_tag = intergration[gene][method][ref][coor][3]
                    if ref == 'G' and method == 'P':
                        c_genome_start = coor[0]
                        c_genome_end = coor[1]
                        c_phase_score = intergration[gene][method][ref][coor][1]
                        c_phase_ratio = intergration[gene][method][ref][coor][0]
                        c_tag = intergration[gene][method][ref][coor][3]

                    for ref1 in ['F', 'G', 'C']:
                        for method1 in ['P', 'H']:
                            if method == method1 and ref == ref1:
                                continue
                            for coor1 in list(intergration[gene][method1][ref1].keys()):
                                if intergration[gene][method1][ref1][coor1] == '-':
                                    continue
                                if Overlap(coor, coor1, method+ref, method1+ref1, gene, trans):
                                    method_ref.append(method1+ref1)
                                    if ref == 'C' and method == 'H':
                                        c_transcript_start = coor[0]
                                        c_transcript_end = coor[1]
                                        c_pvalue = intergration[gene][method][ref][coor][2]
                                        c_genome_start = trans[gene][coor[0]]
                                        c_genome_end = trans[gene][coor[1]]
                                        c_tag = intergration[gene][method][ref][coor][3]
                                    if ref == 'C' and method == 'P':
                                        c_transcript_start = coor[0]
                                        c_transcript_end = coor[1]
                                        c_phase_score = intergration[gene][method][ref][coor][1]
                                        c_phase_ratio = intergration[gene][method][ref][coor][0]
                                        c_genome_start = trans[gene][coor[0]]
                                        c_genome_end = trans[gene][coor[1]]
                                        c_tag = intergration[gene][method][ref][coor][3]
                                    if ref == 'F' and method == 'H':
                                        c_transcript_start = coor[0]
                                        c_transcript_end = coor[1]
                                        c_pvalue = intergration[gene][method][ref][coor][2]
                                        c_tag = intergration[gene][method][ref][coor][3]
                                    if ref == 'F' and method == 'P':
                                        c_transcript_start = coor[0]
                                        c_transcript_end = coor[1]
                                        c_phase_score = intergration[gene][method][ref][coor][1]
                                        c_phase_ratio = intergration[gene][method][ref][coor][0]
                                        c_tag = intergration[gene][method][ref][coor][3]
                                    if ref == 'G' and method == 'H':
                                        c_genome_start = coor[0]
                                        c_genome_end = coor[1]
                                        c_pvalue = intergration[gene][method][ref][coor][2]
                                        c_tag = intergration[gene][method][ref][coor][3]
                                    if ref == 'G' and method == 'P':
                                        c_genome_start = coor[0]
                                        c_genome_end = coor[1]
                                        c_phase_score = intergration[gene][method][ref][coor][1]
                                        c_phase_ratio = intergration[gene][method][ref][coor][0]
                                        c_tag = intergration[gene][method][ref][coor][3]
                                    
                                    if ref1 == 'C' and method1 == 'H':
                                        c1_transcript_start = coor1[0]
                                        c1_transcript_end = coor1[1]
                                        c1_pvalue = intergration[gene][method1][ref1][coor1][2]
                                        c1_tag = intergration[gene][method1][ref1][coor1][3]
                                    if ref1 == 'C' and method1 == 'P':
                                        c1_transcript_start = coor1[0]
                                        c1_transcript_end = coor1[1]
                                        c1_phase_score = intergration[gene][method1][ref1][coor1][1]
                                        c1_phase_ratio = intergration[gene][method1][ref1][coor1][0]
                                        c1_tag = intergration[gene][method1][ref1][coor1][3]
                                    if ref1 == 'F' and method1 == 'H':
                                        c1_transcript_start = coor1[0]
                                        c1_transcript_end = coor1[1]
                                        c1_pvalue = intergration[gene][method1][ref1][coor1][2]
                                        c1_tag = intergration[gene][method1][ref1][coor1][3]
                                    if ref1 == 'F' and method1 == 'P':
                                        c1_transcript_start = coor1[0]
                                        c1_transcript_end = coor1[1]
                                        c1_phase_score = intergration[gene][method1][ref1][coor1][1]
                                        c1_phase_ratio = intergration[gene][method1][ref1][coor1][0]
                                        c1_tag = intergration[gene][method1][ref1][coor1][3]
                                    if ref1 == 'G' and method1 == 'H':
                                        c1_genome_start = coor1[0]
                                        c1_genome_end = coor1[1]
                                        c1_pvalue = intergration[gene][method1][ref1][coor1][2]
                                        c1_tag = intergration[gene][method1][ref1][coor1][3]
                                    if ref1 == 'G' and method1 == 'P':
                                        c1_genome_start = coor1[0]
                                        c1_genome_end = coor1[1]
                                        c1_phase_score = intergration[gene][method1][ref1][coor1][1]
                                        c1_phase_ratio = intergration[gene][method1][ref1][coor1][0]
                                        c1_tag = intergration[gene][method1][ref1][coor1][3]

                                    intergration[gene][method1][ref1][coor1] = '-'
                                else:
                                    # method_ref.append(method+ref)
                                    # method_ref.append(method1+ref1)

                                    if ref == 'C' and method == 'H':
                                        c_transcript_start = coor[0]
                                        c_transcript_end = coor[1]
                                        c_pvalue = intergration[gene][method][ref][coor][2]
                                        c_genome_start = trans[gene][coor[0]]
                                        c_genome_end = trans[gene][coor[1]]
                                        c_tag = intergration[gene][method][ref][coor][3]
                                    if ref == 'C' and method == 'P':
                                        c_transcript_start = coor[0]
                                        c_transcript_end = coor[1]
                                        c_phase_score = intergration[gene][method][ref][coor][1]
                                        c_phase_ratio = intergration[gene][method][ref][coor][0]
                                        c_genome_start = trans[gene][coor[0]]
                                        c_genome_end = trans[gene][coor[1]]
                                        c_tag = intergration[gene][method][ref][coor][3]
                                    if ref == 'F' and method == 'H':
                                        c_transcript_start = coor[0]
                                        c_transcript_end = coor[1]
                                        c_pvalue = intergration[gene][method][ref][coor][2]
                                        c_tag = intergration[gene][method][ref][coor][3]
                                    if ref == 'F' and method == 'P':
                                        c_transcript_start = coor[0]
                                        c_transcript_end = coor[1]
                                        c_phase_score = intergration[gene][method][ref][coor][1]
                                        c_phase_ratio = intergration[gene][method][ref][coor][0]
                                        c_tag = intergration[gene][method][ref][coor][3]
                                    if ref == 'G' and method == 'H':
                                        c_genome_start = coor[0]
                                        c_genome_end = coor[1]
                                        c_pvalue = intergration[gene][method][ref][coor][2]
                                        c_tag = intergration[gene][method][ref][coor][3]
                                    if ref == 'G' and method == 'P':
                                        c_genome_start = coor[0]
                                        c_genome_end = coor[1]
                                        c_phase_score = intergration[gene][method][ref][coor][1]
                                        c_phase_ratio = intergration[gene][method][ref][coor][0]
                                        c_tag = intergration[gene][method][ref][coor][3]

                    if c_genome_start != '-':
                        genome_start = c_genome_start
                    else:
                        genome_start = c1_genome_start
                    if c_genome_end != '-':
                        genome_end = c_genome_end
                    else:
                        genome_end = c1_genome_end
                    if c_transcript_start != '-':
                        transcript_start = c_transcript_start
                    else:
                        transcript_start = c1_transcript_start

                    if c_transcript_end != '-':
                        transcript_end = c_transcript_end
                    else:
                        transcript_end = c1_transcript_end

                    if c_pvalue != '-':
                        pvalue = c_pvalue
                    else:
                        pvalue = c1_pvalue
                    if c_phase_score != '-':
                        phase_score = c_phase_score
                    else:
                        phase_score = c1_phase_score
                    if c_phase_ratio != '-':
                        phase_ratio = c_phase_ratio
                    else:
                        phase_ratio = c1_phase_ratio
                    
                    if c_tag != '-':
                        tag = c_tag
                    else:
                        tag = c1_tag
                    
                    list_.append(f"{relation_dic[gene]}\t{gene}\t{genome_start}\t{genome_end}\t{transcript_start}\t{transcript_end}\t{','.join(list(set(method_ref)))}\t{pvalue}\t{phase_score}\t{phase_ratio}\t{tag}")
                    intergration[gene][method][ref][coor] = '-'

    for gene in list(intergration_FLNC.keys()):
        for ref in ['F']:
            for method in ['H', 'P']:
                for coor in list(intergration_FLNC[gene][method][ref].keys()):
                    genome_start, genome_end, transcript_start, transcript_end, pvalue, phase_score, phase_ratio, tag = \
                    '-', '-', '-', '-', '-', '-', '-', '-'
                    c_genome_start, c_genome_end, c_transcript_start, c_transcript_end, c_pvalue, c_phase_score, c_phase_ratio, c_tag = \
                    '-', '-', '-', '-', '-', '-', '-', '-'
                    c1_genome_start, c1_genome_end, c1_transcript_start, c1_transcript_end, c1_pvalue, c1_phase_score, c1_phase_ratio, c1_tag = \
                    '-', '-', '-', '-', '-', '-', '-', '-'
                    method_ref = []
                    method_ref.append(method+ref)
                    if intergration_FLNC[gene][method][ref][coor] == '-':
                        continue
                    if ref == 'F' and method == 'H':
                        c_transcript_start = coor[0]
                        c_transcript_end = coor[1]
                        c_pvalue = intergration_FLNC[gene][method][ref][coor][2]
                        c_tag = intergration_FLNC[gene][method][ref][coor][3]
                    if ref == 'F' and method == 'P':
                        c_transcript_start = coor[0]
                        c_transcript_end = coor[1]
                        c_phase_score = intergration_FLNC[gene][method][ref][coor][1]
                        c_phase_ratio = intergration_FLNC[gene][method][ref][coor][0]
                        c_tag = intergration_FLNC[gene][method][ref][coor][3]

                    for ref1 in ['F']:
                        for method1 in ['P', 'H']:
                            if method == method1 and ref == ref1:
                                continue
                            for coor1 in list(intergration_FLNC[gene][method1][ref1].keys()):
                                if intergration_FLNC[gene][method1][ref1][coor1] == '-':
                                    continue
                                if Overlap(coor, coor1, method+ref, method1+ref1, gene, trans):
                                    method_ref.append(method1+ref1)
                                    if ref == 'F' and method == 'H':
                                        c_transcript_start = coor[0]
                                        c_transcript_end = coor[1]
                                        c_pvalue = intergration_FLNC[gene][method][ref][coor][2]
                                        c_tag = intergration_FLNC[gene][method][ref][coor][3]
                                    if ref == 'F' and method == 'P':
                                        c_transcript_start = coor[0]
                                        c_transcript_end = coor[1]
                                        c_phase_score = intergration_FLNC[gene][method][ref][coor][1]
                                        c_phase_ratio = intergration_FLNC[gene][method][ref][coor][0]
                                        c_tag = intergration_FLNC[gene][method][ref][coor][3]
                                    
                                    if ref1 == 'F' and method1 == 'H':
                                        c1_transcript_start = coor1[0]
                                        c1_transcript_end = coor1[1]
                                        c1_pvalue = intergration_FLNC[gene][method1][ref1][coor1][2]
                                        c1_tag = intergration_FLNC[gene][method1][ref1][coor1][3]
                                    if ref1 == 'F' and method1 == 'P':
                                        c1_transcript_start = coor1[0]
                                        c1_transcript_end = coor1[1]
                                        c1_phase_score = intergration_FLNC[gene][method1][ref1][coor1][1]
                                        c1_phase_ratio = intergration_FLNC[gene][method1][ref1][coor1][0]
                                        c1_tag = intergration_FLNC[gene][method1][ref1][coor1][3]

                                    intergration_FLNC[gene][method1][ref1][coor1] = '-'
                                else:
                                    # method_ref.append(method+ref)
                                    # method_ref.append(method1+ref1)

                                    if ref == 'F' and method == 'H':
                                        c_transcript_start = coor[0]
                                        c_transcript_end = coor[1]
                                        c_pvalue = intergration_FLNC[gene][method][ref][coor][2]
                                        c_tag = intergration_FLNC[gene][method][ref][coor][3]
                                    if ref == 'F' and method == 'P':
                                        c_transcript_start = coor[0]
                                        c_transcript_end = coor[1]
                                        c_phase_score = intergration_FLNC[gene][method][ref][coor][1]
                                        c_phase_ratio = intergration_FLNC[gene][method][ref][coor][0]
                                        c_tag = intergration_FLNC[gene][method][ref][coor][3]

                    if c_transcript_start != '-':
                        transcript_start = c_transcript_start
                    else:
                        transcript_start = c1_transcript_start

                    if c_transcript_end != '-':
                        transcript_end = c_transcript_end
                    else:
                        transcript_end = c1_transcript_end

                    if c_pvalue != '-':
                        pvalue = c_pvalue
                    else:
                        pvalue = c1_pvalue
                    if c_phase_score != '-':
                        phase_score = c_phase_score
                    else:
                        phase_score = c1_phase_score
                    if c_phase_ratio != '-':
                        phase_ratio = c_phase_ratio
                    else:
                        phase_ratio = c1_phase_ratio
                    
                    if c_tag != '-':
                        tag = c_tag
                    else:
                        tag = c1_tag

                    list_.append(f"{flnc_anno_dic[gene]}\t{gene}\t{genome_start}\t{genome_end}\t{transcript_start}\t{transcript_end}\t{','.join(list(set(method_ref)))}\t{pvalue}\t{phase_score}\t{phase_ratio}\t{tag}")
                    intergration_FLNC[gene][method][ref][coor] = '-'
    return list_



def Overlap(coor, coor1, mr1, mr2, gene, trans):
    if 'C' in mr1 and 'C' in mr2:
        list1 = set(list(range(coor[0], coor[1])))
        list2 = set(list(range(coor1[0], coor1[1])))
        set_ = list1&list2
        if len(set_) >= 1:
            return True
        else:
            return False
    if 'G' in mr1 and 'G' in mr2:
        list1 = set(list(range(coor[0], coor[1])))
        list2 = set(list(range(coor1[0], coor1[1])))
        set_ = list1&list2
        if len(set_) >= 1:
            return True
        else:
            return False
    if 'F' in mr1 and 'F' in mr2:
        list1 = set(list(range(coor[0], coor[1])))
        list2 = set(list(range(coor1[0], coor1[1])))
        set_ = list1&list2
        if len(set_) >= 1:
            return True
        else:
            return False
    if 'C' in mr1 and 'HG' in mr2:
        list1 = [trans[gene][coor[0]], trans[gene][coor[1]]]
        try:
            list1 = sorted(list1)
        except TypeError:
            return False
        list2 = [coor1[0], coor1[1]]
        if IsOverlap(list1, list2):
            return True
        else:
            return False

    if 'C' in mr1 and 'PG' in mr2:
        list1 = [trans[gene][coor[0]], trans[gene][coor[1]]]
        try:
            list1 = sorted(list1)
        except TypeError:
            return False
        list2 = [coor1[0], coor1[1]]
        if IsOverlap(list1, list2):
            return True
        else:
            return False
    
    if 'C' in mr1 and 'F' in mr2:
        return True
    
    if 'G' in mr1 and 'F' in mr2:
        return True

def Intergenic_PHAS_Loci(intergrationfile, tag_dic):
    trans = {}
    intergration = nestedDic()
    with open(intergrationfile, 'r') as fn:
        for line in fn:
            line1 = line.strip()
            if line1.startswith('cluster_number'):
                continue
            l = line.strip().split("\t")
            cluster_number,feature,geneid,method,reference_type,h_transcriptome_coordinate,h_genome_coordinate,h_fl_transcriptome_coordinate,p_transcriptome_coordinate,p_genome_coordinate,p_fl_transcriptome_coordinate,c_pvalue,g_pvalue,f_pvalue,c_phaseratio,g_phaseratio,f_phaseratio,c_phasescore,g_phasescore,f_phasescore = \
                l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8],l[9],l[10],l[11],l[12],l[13],l[14],l[15],l[16],l[17],l[18],l[19]

            if 'H' in method:
                for h_coor in h_genome_coordinate.split(','):
                    if h_coor != '-':
                        idx = h_genome_coordinate.split(',').index(h_coor)
                        pvalue = g_pvalue.split(',')[idx]
                        coor0 = int(h_coor.replace('(', '').replace(')', '').split(':')[0])
                        coor1 = int(h_coor.replace('(', '').replace(')', '').split(':')[1])
                        query = geneid.split('_')[0] + '_' + geneid.split('_')[1] + '\t' + 'HG' + '\t' + str(coor0) + '\t' + str(coor1)
                        if query in tag_dic:
                            tmp_tag = tag_dic[query]
                        else:
                            tmp_tag = '-'
                        intergration[geneid]['H']['G'][(coor0, coor1)] = ['-', '-', pvalue, tmp_tag]
            if 'P' in method:
                for p_coor in p_genome_coordinate.split(','):
                    if p_coor != '-':
                        idx = p_genome_coordinate.split(',').index(p_coor)
                        phase_score = g_phasescore.split(',')[idx]
                        phase_ratio = g_phaseratio.split(',')[idx]
                        coor0 = int(p_coor.replace('(', '').replace(')', '').split(':')[0])
                        coor1 = int(p_coor.replace('(', '').replace(')', '').split(':')[1])
                        query = geneid.split('_')[0] + '_' + geneid.split('_')[1] + '\t' + 'PG' + '\t' + str(coor0) + '\t' + str(coor1)
                        if query in tag_dic:
                            tmp_tag = tag_dic[query]
                        else:
                            tmp_tag = '-'
                        intergration[geneid]['P']['G'][(coor0, coor1)] = [phase_ratio, phase_score, '-', tmp_tag]

    list_ = []
    order = 0
    for gene in list(intergration.keys()):
        for ref in ['G']:
            for method in ['H', 'P']:
                for coor in list(intergration[gene][method][ref].keys()):
                    genome_start, genome_end, transcript_start, transcript_end, pvalue, phase_score, phase_ratio, tag = \
                    '-', '-', '-', '-', '-', '-', '-', '-'
                    c_genome_start, c_genome_end, c_transcript_start, c_transcript_end, c_pvalue, c_phase_score, c_phase_ratio, c_tag = \
                    '-', '-', '-', '-', '-', '-', '-', '-'
                    c1_genome_start, c1_genome_end, c1_transcript_start, c1_transcript_end, c1_pvalue, c1_phase_score, c1_phase_ratio, c1_tag = \
                    '-', '-', '-', '-', '-', '-', '-', '-'
                    method_ref = []
                    method_ref.append(method+ref)
                    if intergration[gene][method][ref][coor] == '-':
                        continue
                    if ref == 'G' and method == 'H':
                        c_genome_start = coor[0]
                        c_genome_end = coor[1]
                        c_pvalue = intergration[gene][method][ref][coor][2]
                        c_tag = intergration[gene][method][ref][coor][3]
                    if ref == 'G' and method == 'P':
                        c_genome_start = coor[0]
                        c_genome_end = coor[1]
                        c_phase_score = intergration[gene][method][ref][coor][1]
                        c_phase_ratio = intergration[gene][method][ref][coor][0]
                        c_tag = intergration[gene][method][ref][coor][3]

                    for ref1 in ['G']:
                        for method1 in ['P', 'H']:
                            if method == method1 and ref == ref1:
                                continue
                            for coor1 in list(intergration[gene][method1][ref1].keys()):
                                if intergration[gene][method1][ref1][coor1] == '-':
                                    continue
                                if Overlap(coor, coor1, method+ref, method1+ref1, gene, trans):
                                    method_ref.append(method1+ref1)
                                    if ref == 'G' and method == 'H':
                                        c_genome_start = coor[0]
                                        c_genome_end = coor[1]
                                        c_pvalue = intergration[gene][method][ref][coor][2]
                                        c_tag = intergration[gene][method][ref][coor][3]
                                    if ref == 'G' and method == 'P':
                                        c_genome_start = coor[0]
                                        c_genome_end = coor[1]
                                        c_phase_score = intergration[gene][method][ref][coor][1]
                                        c_phase_ratio = intergration[gene][method][ref][coor][0]
                                        c_tag = intergration[gene][method][ref][coor][3]
                                    
                                    if ref1 == 'G' and method1 == 'H':
                                        c1_genome_start = coor1[0]
                                        c1_genome_end = coor1[1]
                                        c1_pvalue = intergration[gene][method1][ref1][coor1][2]
                                        c1_tag = intergration[gene][method1][ref1][coor1][3]
                                    if ref1 == 'G' and method1 == 'P':
                                        c1_genome_start = coor1[0]
                                        c1_genome_end = coor1[1]
                                        c1_phase_score = intergration[gene][method1][ref1][coor1][1]
                                        c1_phase_ratio = intergration[gene][method1][ref1][coor1][0]
                                        c1_tag = intergration[gene][method1][ref1][coor1][3]

                                    intergration[gene][method1][ref1][coor1] = '-'
                                else:
                                    # method_ref.append(method+ref)
                                    # method_ref.append(method1+ref1)

                                    if ref == 'G' and method == 'H':
                                        c_genome_start = coor[0]
                                        c_genome_end = coor[1]
                                        c_pvalue = intergration[gene][method][ref][coor][2]
                                        c_tag = intergration[gene][method][ref][coor][3]
                                    if ref == 'G' and method == 'P':
                                        c_genome_start = coor[0]
                                        c_genome_end = coor[1]
                                        c_phase_score = intergration[gene][method][ref][coor][1]
                                        c_phase_ratio = intergration[gene][method][ref][coor][0]
                                        c_tag = intergration[gene][method][ref][coor][3]

                    if c_genome_start != '-':
                        genome_start = c_genome_start
                    else:
                        genome_start = c1_genome_start
                    if c_genome_end != '-':
                        genome_end = c_genome_end
                    else:
                        genome_end = c1_genome_end

                    if c_pvalue != '-':
                        pvalue = c_pvalue
                    else:
                        pvalue = c1_pvalue
                    if c_phase_score != '-':
                        phase_score = c_phase_score
                    else:
                        phase_score = c1_phase_score
                    if c_phase_ratio != '-':
                        phase_ratio = c_phase_ratio
                    else:
                        phase_ratio = c1_phase_ratio
                    
                    if c_tag != '-':
                        tag = c_tag
                    else:
                        tag = c1_tag

                    order += 1
                    new_gene = gene.split('_')[0]+'_'+gene.split('_')[1] + '#' + str(order)
                    list_.append(f"Intergenic\t{new_gene}\t{genome_start}\t{genome_end}\t{transcript_start}\t{transcript_end}\t{','.join(list(set(method_ref)))}\t{pvalue}\t{phase_score}\t{phase_ratio}\t{tag}")
                    intergration[gene][method][ref][coor] = '-'
    return list_

def PHAS_Loci_out_write(PHAS_Loci_out, PHAS_Loci, PHAS_Loci1, passP):
    recorder = 0
    with open(PHAS_Loci_out, 'w') as fo:
        fo.write(f'feature\tPHAS_Loci\tGenome_start\tGenome_end\ttranscript_start\ttranscript_end\tmethod_ref\tpvalue\tphase_score\tphase_ratio\ttag\trecorder\n')
        for i in PHAS_Loci:
            l = i.split('\t')
            method_ref = l[6]
            if passP == 'y' and 'H' not in method_ref:
                continue
            if i.startswith('feature'):
                fo.write(f'{i}\trecorder\n')
                continue
            recorder += 1
            fo.write(f'{i}\t{recorder}\n')
        for i in PHAS_Loci1:
            l = i.split('\t')
            method_ref = l[6]
            if passP == 'y' and 'H' not in method_ref:
                continue
            recorder += 1
            fo.write(f'{i}\t{recorder}\n')

def SplitIsland1(candidate_cluster0, phase_length=21):
    candidate_cluster = nestedDic()
    parent_order = 1
    for i in candidate_cluster0:
        dic = [[]]
        order = 0
        object = (dic, order, candidate_cluster0[i])
        out_list = literation_func(object)
        for j in out_list:
            if len(j) != 0:
                candidate_cluster[parent_order] = j
                parent_order += 1

    return candidate_cluster

def literation_func(lists, phase_length=21): 
    if phase_length == 21:
        candidate_list = [0, 2, 19]
    elif phase_length == 24:
        candidate_list = [0, 2, 19]
    dic = lists[0]
    order = lists[1]
    list = lists[2]
    left_list = []
    if len(list) != 0:
        former  = list[0]
    else:
        return dic


    for i in list:
        test = (i - former) % phase_length
        if test in candidate_list:
            dic[order].append(i)
            former = i
        else:
            left_list.append(i)
    order += 1
    dic.append([])
    
    litetation_object = (dic, order, left_list)
    return literation_func(litetation_object, phase_length)

def GetBaseName(dir: str):
    """return filename

    Parameters
    ----------
    dir : directory
        filename directory
    Returns
    -------
    str
        filename
    """    
    return os.path.basename(dir)
