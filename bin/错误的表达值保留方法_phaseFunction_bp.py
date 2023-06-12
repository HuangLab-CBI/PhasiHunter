# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~> import modules <~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
import time
import math
import os
from Bio import SeqIO
from sympy import binomial
from scipy.special import comb
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, wait, ALL_COMPLETED
# --------------------------------------------------------------> function defination <---------------------------------------------------------------- #
def GetCurTime():
    """return formated current localtime 

    Returns:
        <str> -- formated current localtime
    """
    ctime = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
    return ctime

def defauldict_list():
    return defaultdict(list)

def TwoDepDic():
    return defaultdict(defauldict_list)

def ThreeDepDic():
    return defaultdict(TwoDepDic)

def FourDepDic():
    return defaultdict(ThreeDepDic)

def FiveDepDic():
    return defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list)))))

def nestedDic():
    """generating nested dictory

    Returns
    -------
    dict
        nested dictory
    """    
    return defaultdict(nestedDic)

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

def ParseRef(file: str, type: str, method = ''):
    """for parse reference fasta file

    Arguments:
        file {str} -- reference sequence fasta
        type {str} -- cdna | gdna

    Raises:
        AttributeError: _description_

    Returns:
        dic -- for cdna: {transcriptID: [start_position, end_position, sequence, annotation]}
                for gdna: {chromesomeID: [start_position, end_position, annotation]}
    """
    dic = {}
    if type == 'cdna' or type == 'flnc':
        for query in SeqIO.parse(file,"fasta"):
            name = query.id
            anno = " ".join((query.description).split(" ")[1:])
            seq = str(query.seq)
            seqLength = len(seq)
            if name not in dic:
                if method == 'PhaseScore':
                    dic[name] = [0, seqLength - 1, seq, anno]
                else:
                    dic[name] = [1, seqLength, seq, anno]
            else:
                raise AttributeError('exits same transcript in ref sequence fasta file')
        return dic
    elif type == 'gdna':
        for query in SeqIO.parse(file,"fasta"):
            name = query.id
            anno = " ".join((query.description).split(" ")[1:])
            seq = str(query.seq)
            seqLength = len(seq)
            if name not in dic:
                if method == 'PhaseScore':
                    dic[name] = [0, seqLength - 1, anno]
                else:
                    dic[name] = [1, seqLength, anno]
            else:
                raise AttributeError('exits same transcript in ref sequence fasta file')
        return dic
    else:
        raise Exception('input file type erro, must cnda or gdna')

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
        grandSonList = ['exon', 'CDS', 'miRNA', 'intron']
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
                if feature == 'CDS':
                    exonid = exonid + '-' + str(start) + '-' + str(end)
                try:
                    dic[curGrandParentFeature][curParentFeature][curSonFeature][exonid] = {'attribute': [start, end]}
                except KeyError:
                    try:
                        curSonFeature = curParentFeature + 'rna'
                        dic[curGrandParentFeature][curParentFeature][curSonFeature] = {'attribute': [start, end]}
                        dic[curGrandParentFeature][curParentFeature][curSonFeature][exonid] = {'attribute': [start, end]}
                    except KeyError:
                        curParentFeature = curGrandParentFeature + 'gene'
                        dic[curGrandParentFeature][curParentFeature] = {'attribute': [start, end]}
                        curSonFeature = curParentFeature + 'rna'
                        dic[curGrandParentFeature][curParentFeature][curSonFeature] = {'attribute': [start, end]}
                        dic[curGrandParentFeature][curParentFeature][curSonFeature][exonid] = 0
    return dic

def ParseHypergeometricMap(map, phase_length):
    """parse map file

    Parameters
    ----------
    map : file
        map file which generated by bowtie

    Returns
    -------
    dict
        key - geneid+'\t'+strand+'\t'+str(sRNA_pos)+\n
        value - str(sRNA_num)+'\t'+sRNA+'\t'+line1[4]+'\t'+str(sRNA_len)
    """
    genewithhits = {}
    n = 0
    for line in open(map, 'r'):
        n+=1
        line1=line.split('\t')
        sRNA=line1[0]
        try:
            sRNA_num=float(sRNA.split('@')[1])
        except IndexError:
            sRNA_num=float(sRNA.split('_x')[1])
        strand=line1[1]
        geneid=line1[2]
        sRNA_pos=int(line1[3])
        sRNA_len = int(len(line1[4]))
        key = geneid+'\t'+strand+'\t'+str(sRNA_pos)
        if (key) in genewithhits.keys():
            no=float((genewithhits[key]).split('\t')[0])
            no+=sRNA_num
            if sRNA_len == phase_length:
                genewithhits[key]=str(sRNA_num)+'\t'+sRNA+'\t'+line1[4]+'\t'+str(sRNA_len)
            elif genewithhits[key].split('\t')[3] == str(phase_length):
                #WARNING: bug 来源地,因为这里只是保留了phase length，但是没有保留和phase length 相匹配的序列， 因此sRNA序列被覆盖了，所以结果表现出 输出的序列和打印的序列长度不一致的情况。
                abun = genewithhits[key].split('\t')[0]
                seq = genewithhits[key].split('\t')[2]
                genewithhits[key]=str(abun)+'\t'+sRNA+'\t'+ seq +'\t'+str(phase_length)
            else:
                genewithhits[key]=str(no)+'\t'+sRNA+'\t'+line1[4]+'\t'+str(sRNA_len)
        else:
            genewithhits[key]=str(sRNA_num)+'\t'+sRNA+'\t'+line1[4]+'\t'+str(sRNA_len)
    return genewithhits

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
    elif phase_length == 24:
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

def PairPos(list_: list):
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
    loop_list = list_
    for i in loop_list:
        former = list_.pop(0)
        if len(list_) >= 2:
            out_list.append((former, list_[-1]))
            out_list.append((former, list_[-2]))
        else:
            break
        laster = list_.pop(-1)
        out_list.append((laster, list_[0]))
        try:
            out_list.append((laster, list_[1]))
        except IndexError:
            break
        if len(list_) == 0:
            break
    return out_list

def IsOverlap(list1: list, list2: list):
    """ 判断两个list是否有overlap

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

def AnnoGdna(dic_: dict, refgff3: dict):
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
    for chr in dic_:
        if chr in refgff3:
            for geneid in refgff3[chr]:
                if geneid == 'attribute':
                    continue
                refgff3[chr][geneid]['attribute']
                for cluster in dic_[chr]:
                    if cluster not in out_cluster_anno[chr]:
                        out_cluster_anno[chr][cluster] = []
                    if IsOverlap(dic_[chr][cluster], refgff3[chr][geneid]['attribute']):
                        for rnaid in refgff3[chr][geneid]:
                            if rnaid == 'attribute':
                                continue
                            if IsOverlap(dic_[chr][cluster], refgff3[chr][geneid][rnaid]['attribute']):
                                for exonid in refgff3[chr][geneid][rnaid]:
                                    if exonid == 'attribute':
                                        continue
                                    if IsOverlap(dic_[chr][cluster], refgff3[chr][geneid][rnaid][exonid]['attribute']):
                                        out_cluster_anno[chr][cluster].append(exonid)
    return out_cluster_anno

def FirstScaning(mapdata: dict, setting: tuple):
    """first scaning

    Parameters
    ----------
    mapdata : dict
        _description_
    setting : tuple
        setting = (windowLength, phase_length, phase_number, pvalue_cutoff, min_read_num, type)

    Returns
    -------
    tuple
        (final_cluster, final_strand, final_allsiRNA)
    """
    windowLength = setting[0]
    phase_length = setting[1]
    phase_number = setting[2]
    pvalue_cutoff = setting[3]
    min_read_num = setting[4]
    type_ = setting[5]
    record_marker = 0
    # final_cluster = defaultdict(list)
    # final_strand = nestedDic()
    final_phasiRNA = nestedDic()
    final_allsiRNA = nestedDic()
    for k in mapdata:
        max_readn = 0
        pvalue = 0.0
        line = mapdata[k].split('\t')
        key = k.split('\t')
        sRNA_no = float(line[0])
        sRNA = line[1]
        strand = key[1]
        geneid = key[0]
        sRNA_pos = int(key[2])
        total_n = 0
        total_k = 0
        phasiRNA = []
        allsiRNA = []
        if strand == '+':
            for i in range(sRNA_pos, sRNA_pos + windowLength):
                query = geneid + '\t' + '+' + '\t' + str(i)
                if query in mapdata:
                    total_n += 1
                    allsiRNA.append(query + '\t' + mapdata[query])
            for i in range(sRNA_pos - 2, sRNA_pos + windowLength - 2):
                query = geneid + '\t' + '-' + '\t' + str(i)
                if query in mapdata:
                    total_n += 1
                    allsiRNA.append(query + '\t' + mapdata[query])
            for i in range(sRNA_pos, sRNA_pos + windowLength, phase_length):
                query = geneid + '\t' + '+' + '\t' + str(i)
                if query in mapdata and mapdata[query].split('\t')[-1] == str(phase_length):
                    total_k += 1
                    phasiRNA.append(query + '\t' + mapdata[query])
            for i in range(sRNA_pos - 2, sRNA_pos + windowLength - 2, phase_length):
                query = geneid + '\t' + '-' + '\t' + str(i)
                if query in mapdata and mapdata[query].split('\t')[-1] == str(phase_length):
                    total_k += 1
                    phasiRNA.append(query + '\t' + mapdata[query])
        elif strand == '-':
            for i in range(sRNA_pos, sRNA_pos + windowLength):
                query = geneid + '\t' + '-' + '\t' + str(i)
                if query in mapdata:
                    total_n += 1
                    allsiRNA.append(query + '\t' + mapdata[query])
            for i in range(sRNA_pos + 2, sRNA_pos + windowLength + 2):
                query = geneid + '\t' + '+' + '\t' + str(i)
                if query in mapdata:
                    total_n += 1
                    allsiRNA.append(query + '\t' + mapdata[query])
            for i in range(sRNA_pos, sRNA_pos + windowLength, phase_length):
                query = geneid + '\t' + '-' + '\t' + str(i)
                if query in mapdata and mapdata[query].split('\t')[-1] == str(phase_length):
                    total_k += 1
                    phasiRNA.append(query + '\t' + mapdata[query])
            for i in range(sRNA_pos + 2, sRNA_pos + windowLength + 2, phase_length):
                query = geneid + '\t' + '+' + '\t' + str(i)
                if query in mapdata and mapdata[query].split('\t')[-1] == str(phase_length):
                    total_k += 1
                    phasiRNA.append(query + '\t' + mapdata[query])

        # phase_number and pvalue cutoff
        if total_k >= phase_number:
            for i in range(total_k, phase_length):
                pvalue += float(binomial(windowLength * 2 - 1 - phase_length, total_n - i) * binomial(phase_length, i) / binomial(windowLength * 2 - 1, total_n))
            if pvalue <= pvalue_cutoff and pvalue != 0.0:
                for i in range(len(phasiRNA)):
                    value_list = phasiRNA[i].split('\t')
                    sRNA_no = float(value_list[3])
                    if sRNA_no > max_readn:
                        max_readn = sRNA_no

                if min_read_num == 0 or max_readn >= min_read_num:  # min_read_num cutoff
                    record_marker += 1
                    if type_ == 'cdna':
                        out_record_marker = "CF"+str(record_marker) 
                    elif type_ == 'gdna':
                        out_record_marker = "GF"+str(record_marker)
                    elif type_ == 'flnc':
                        out_record_marker = "FF"+str(record_marker)
                    for i in range(len(phasiRNA)):
                        value_list = phasiRNA[i].split('\t')
                        geneid = value_list[0]
                        strand = value_list[1]
                        readpos = int(value_list[2])
                        query = geneid + '\t' + strand + '\t' + str(readpos)
                        final_phasiRNA[geneid][query][record_marker] = phasiRNA[i] + '\t' + str('%.3e' % pvalue) + '\t' + geneid+"_"+str(out_record_marker)
                        # final_cluster[geneid].append(readpos)
                        # final_strand[geneid][readpos] = strand
                    for i in range(len(allsiRNA)):
                        value_list = allsiRNA[i].split('\t')
                        geneid = value_list[0]
                        readpos = int(value_list[2])
                        strand = value_list[1]
                        query = geneid + '\t' + strand + '\t' + str(readpos)
                        final_allsiRNA[geneid][query][record_marker] = allsiRNA[i] + "\t" + str('%.3e' % pvalue) + '\t' + geneid+"_"+str(out_record_marker)
                        # if query not in final_allsiRNA[geneid]:
                        #     final_allsiRNA[geneid][query] = allsiRNA[i] + "\t" + str(pvalue)
                        # elif query in final_allsiRNA[geneid] and pvalue < float(final_allsiRNA[geneid][query].split("\t")[-1]):
                        #     final_allsiRNA[geneid][query] = allsiRNA[i] + "\t" + str(pvalue)
    return (final_phasiRNA, final_allsiRNA)

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
    record_marker = 1
    final_cluster = setting[0]
    final_strand = setting[1]
    phase_number = setting[2]
    phase_length = setting[3]
    pvalue_cutoff = setting[4]
    min_read_num = setting[5]
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
            pair_pos = PairPos(list(range(len(candidate_cluster[j]))))
            for tuple_ in pair_pos:
                total_n = 0
                total_k = 0
                pvalue = 0
                if tuple_[0] > tuple_[1]:
                    start_pos = candidate_cluster[j][tuple_[1]]
                    tmp_pos = candidate_cluster[j][tuple_[0]]
                else:
                    start_pos = candidate_cluster[j][tuple_[0]]
                    tmp_pos = candidate_cluster[j][tuple_[1]]
                if abs(tuple_[1] - tuple_[0]) < phase_number:
                    break
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
                        if query in final_allsiRNA[k] and final_allsiRNA[k][query].split('\t')[-2] == str(phase_length):
                            total_k += 1
                            phasiRNA.append(final_allsiRNA[k][query])
                    for i in range(start_pos - 2, end_pos - 2, phase_length):
                        query = k + '\t' + '-' + '\t' + str(i)
                        if query in final_allsiRNA[k] and final_allsiRNA[k][query].split('\t')[-2] == str(phase_length):
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
                        if query in final_allsiRNA[k] and final_allsiRNA[k][query].split('\t')[-2] == str(phase_length):
                            total_k += 1
                            phasiRNA.append(str(final_allsiRNA[k][query]))
                    for i in range(start_pos + 2, end_pos + 2, phase_length):
                        query = k + '\t' + '+' + '\t' + str(i)
                        if query in final_allsiRNA[k] and final_allsiRNA[k][query].split('\t')[-2] == str(phase_length):
                            total_k += 1
                            phasiRNA.append(str(final_allsiRNA[k][query]))

                # phase_number and pvalue cutoff
                if total_k >= phase_number:
                    for i in range(total_k, phase_length):
                        pvalue += float(binomial((end_pos - start_pos) * 2 - 1 - phase_length, total_n - i) * binomial(phase_length, i) / binomial((end_pos - start_pos) * 2 - 1, total_n))
                    if pvalue <= pvalue_cutoff and pvalue != 0.0:
                        for i in range(len(phasiRNA)):
                            value_list = phasiRNA[i].split('\t')
                            sRNA_no = float(value_list[-5])
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
                                    out_phasiRNA[geneid][record_marker][readpos] = "\t".join(phasiRNA[i].split("\t")[:-1]) + "\t" + str(pvalue)
                                else:
                                    final_gene[geneid] = str(max_readn) + '\t' + str(pvalue)
                                    out_phasiRNA[geneid][record_marker][readpos] = "\t".join(phasiRNA[i].split("\t")[:-1]) + "\t" + str(pvalue)
                            start_pos = list(sorted(out_phasiRNA[geneid][record_marker].keys()))[0]
                            end_pos = list(sorted(out_phasiRNA[geneid][record_marker].keys()))[-1]
                            out_phasiRNA_cluster[geneid][record_marker] = [start_pos, end_pos]

                            for i in range(len(allsiRNA)):
                                value_list = allsiRNA[i].split('\t')
                                geneid = value_list[0]
                                readpos = int(value_list[2])
                                query = geneid + '\t' + '+' + '\t' + str(i)
                                out_allsiRNA[geneid][record_marker][readpos] = "\t".join(allsiRNA[i].split("\t")[:-1]) + "\t" + str(pvalue)
                            start_pos = list(sorted(out_allsiRNA[geneid][record_marker].keys()))[0]
                            end_pos = list(sorted(out_allsiRNA[geneid][record_marker].keys()))[-1]
                            out_allsiRNA_cluster[geneid][record_marker] = [start_pos, end_pos]
                if tag:
                    break
            if tag == False:
                for i in range(len(phasiRNA)):
                    value_list = phasiRNA[i].split('\t')
                    geneid = value_list[0]
                    strand = value_list[1]
                    readpos = int(value_list[2])
                    duplication_phasiRNA[geneid][record_marker][readpos] = phasiRNA[i]
                    start_pos = list(sorted(duplication_phasiRNA[geneid][record_marker].keys()))[0]
                    end_pos = list(sorted(duplication_phasiRNA[geneid][record_marker].keys()))[-1]
                    out_duplication_phasiRNA_cluster[geneid][record_marker] = [start_pos, end_pos]
                for i in range(len(allsiRNA)):
                    value_list = allsiRNA[i].split('\t')
                    geneid = value_list[0]
                    readpos = int(value_list[2])
                    duplication_allsiRNA[geneid][record_marker][readpos] = allsiRNA[i]
                    start_pos = list(sorted(duplication_allsiRNA[geneid][record_marker].keys()))[0]
                    end_pos = list(sorted(duplication_allsiRNA[geneid][record_marker].keys()))[-1]
                    out_duplication_allsiRNA_cluster[geneid][record_marker] = [start_pos, end_pos]

    return (out_phasiRNA, out_allsiRNA, out_phasiRNA_cluster, out_allsiRNA_cluster, duplication_phasiRNA, duplication_allsiRNA, out_duplication_phasiRNA_cluster, out_duplication_allsiRNA_cluster)

def TwoStrandMappingDetect(list_: list, dic: dict):
    """detecting two strand mapping

    Arguments:
        list_ {list} -- condicatePhaseCluster, 
        dic {dict} -- {'-': [1196, 1301, 1175, 1343, 52], '+': [1364, 1280, 1406, 1385]}
    
    Return:
        bool
    """

    judge = []
    for i in list_:
        if int(i) in dic['-']:
            judge.append('-')
        if int(i) in dic['+']:
            judge.append('+')

    judge = set(judge)
    if len(judge) == 2:
        return True
    else:
        return False

def AlignBins(start: int, end: int, phaseLength: int):
    """align bin to transcript

    Arguments:
        start {int} -- transcript start position,
        end {int} -- transcript end position,
        phaseLength {int} -- 21 | 24

    Returns:
        dic -- {transcriptPosition: bin, ...}
    """
    dic = {}
    bin = 1
    for i in list(range(start, end)):
        if bin <= phaseLength:
            if i not in dic:
                dic[i] = bin
                bin += 1
        else:
            bin = 1
            if i not in dic:
                dic[i] = bin
                bin += 1
    return dic

def ParseBins(dic: dict, bins: dict, phase_length):
    """Parse bins

    Arguments:
        dic {dict} -- {'attribute': [(id, strand, pos, seq, abun), ...]}
        bins {dict} -- {pos: bin, ...}

    Returns:
        dict -- binDic = {'binAbun': {bin: abun, ...}, 'binPos': {bin: pos,pos, ...}, 'bin21id': {bin: id,id, ...}, 'bin21Pos': {bin: pos,pos, ...}, 'bin21Abun': {bin: abun, ...}, 'hypergeometric': {length: [pos1, pos1, pos2, ...], length: [pos1, pos2, pos3, ...], ...}}
    """

    # %bin_abun = {bin: abun, ..}
    # %bin_pos = {bin: pos,pos,.., ..}
    # %bin_query_21nt = {bin:{pos: sRNA}, ..}
    # %bin_pos_21nt = {bin: pos, ..}
    # %hypergeometric = {bin: {length: [pos1, pos1, pos2, ...], lenght: [pos1, pos2, ...]}, ...}, 最多只能存两个相同的元素， 一个是正链的， 另一个是负链的
    # %bin_strand = {'+': [pos, pos, pos, ...], '-': [pos, pos, pos], ...}
    # %bin_strand_21nt = {'+': [pos, pos, pos, ...], '-': [pos, pos, pos], ...}
    # %bin21out = {pos + "\t" + strand: [sRNAid: seq]}
    # %binout = {pos + "\t" + strand: [sRNAid: seq]}

    binDic = {'binAbun': {}, 'binPos': {}, 'bin21id': {}, 'bin21Pos': {}, 'bin21Abun': {}, 'hypergeometric': {}, 'binStrand': TwoDepDic(), 'bin21Strand': {}, 'bin21out': {}, 'binout': nestedDic(), 'binlenpos': TwoDepDic()}
    for k in dic:
        id = k[0]
        pos = k[2]
        abun = k[4]
        seq = k[3]
        seq_len = len(seq)
        strand = k[1]
        length = len(k[3])
        if strand == '-':
            key = str(pos - 2) + "\t" + strand
            # if key not in binDic['binout']:
            #     binDic['binout'][key] = [id, seq, str(seq_len), str(abun)]
            binDic['binout'][seq_len][key] = [id, seq, str(seq_len), str(abun)]
        elif strand == '+':
            key = str(pos) + "\t" + strand
            # if key not in binDic['binout']:
            #     binDic['binout'][key] = [id, seq, str(seq_len), str(abun)]
            binDic['binout'][seq_len][key] = [id, seq, str(seq_len), str(abun)]
        if length == phase_length:
            if strand == '-':
                key = str(pos - 2) + "\t" + strand
                if key not in binDic['bin21out']:
                    binDic['bin21out'][key] = [id, seq, str(seq_len), str(abun)]
            elif strand == '+':
                key = str(pos) + "\t" + strand
                if key not in binDic['bin21out']:
                    binDic['bin21out'][key] = [id, seq, str(seq_len), str(abun)]
        # if strand not in binDic['binStrand']:
        #     binDic['binStrand'][strand] = [pos]
        # else:
        #     binDic['binStrand'][strand].append(pos)
        binDic['binStrand'][seq_len][strand].append(pos)
        if length == phase_length:
            if strand not in binDic['bin21Strand']:
                binDic['bin21Strand'][strand] = [pos]
            else:
                binDic['bin21Strand'][strand].append(pos)
        if length not in binDic['hypergeometric']:
            binDic['hypergeometric'][length] = [pos]
        else:
            binDic['hypergeometric'][length].append(pos)
        try:
            if bins[pos] not in binDic['binAbun']:
                binDic['binAbun'][bins[pos]] = abun
            else:
                binDic['binAbun'][bins[pos]] += abun
        except KeyError:
            continue
        try:
            if bins[pos] not in binDic['binPos']:
                binDic['binPos'][bins[pos]] = str(pos)
            else:
                binDic['binPos'][bins[pos]] = str(binDic['binPos'][bins[pos]]) + ',' + str(pos)
        except KeyError:
            continue
        binDic['binlenpos'][seq_len][bins[pos]].append(str(pos))
        try:
            if length == phase_length:
                if bins[pos] not in binDic['bin21id']:
                    binDic['bin21id'][bins[pos]] = {}
                    if pos not in binDic['bin21id'][bins[pos]]:
                        binDic['bin21id'][bins[pos]][pos] = id
                    else:
                        binDic['bin21id'][bins[pos]][pos] = binDic['bin21id'][bins[pos]][pos] + ',' + id
                else:
                    if pos not in binDic['bin21id'][bins[pos]]:
                        binDic['bin21id'][bins[pos]][pos] = id
                    else:
                        binDic['bin21id'][bins[pos]][pos] = binDic['bin21id'][bins[pos]][pos] + ',' + id
                if bins[pos] not in binDic['bin21Abun']:
                    binDic['bin21Abun'][bins[pos]] = abun
                else:
                    binDic['bin21Abun'][bins[pos]] += abun
                if bins[pos] not in binDic['bin21Pos']:
                    binDic['bin21Pos'][bins[pos]] = str(pos)
                else:
                    binDic['bin21Pos'][bins[pos]] = str(binDic['bin21Pos'][bins[pos]]) + ',' + str(pos)
        except KeyError:
            continue
    return binDic

def FindMaxAbun(dict_: dict):
    """find max value from a dict

    Arguments:
        dict_ {dict} -- {1:22, 2:33, 3:44, 5:99, ...}

    Returns:
        list -- [(k, maxValue)]
    """
    tmp = [(k,v) for k,v in dict_.items() if v == max(dict_.values())]
    return tmp

def ContainDuplication(list_: list):
    """determine whether duplicate values are included

    Arguments:
        list_ {list} -- _description_

    Returns:
        bool -- True | False
    """
    setLength = len(set(list_))
    listLength = len(list_)
    if setLength < listLength:
        return True
    else:
        return False

def DistanceFilter(condidateCluster: list, phaseNumber: int, maxAbunBinId, island_length):
    """returns the subsets of the cluster that satisfies the distance condition

    Arguments:
        condidateCluster {list} -- the candidate cluster for a transcript

        phaseNumber {int} -- phaseNumber cutoff

        maxAbunBinId {_type_} -- the id set of the candidate transcript cluster

    Returns:
        tuple -- ({subClusterOrder: [id1, id2, id3, ...], ...}, {UnfilterSubClusterOrder: [pos1, pos2, pos3, pos4, ...], ...})
    """
    list_ = sorted(list(set(sorted([int(i) for i in condidateCluster]))))
    firstValue = list_[0]
    tmpDic = {}
    numKey = 1
    for i in list_:
        if i == firstValue:
            formerI = i
            if numKey not in tmpDic:
                tmpDic[numKey] = [firstValue]
            continue
        else:
            distance = i - formerI
            if distance <= island_length:
                tmpDic[numKey].append(i)
                formerI = i
            else:
                numKey += 1
                if numKey not in tmpDic:
                    tmpDic[numKey] = [i]
                    formerI = i
    dic = {k:v for k,v in tmpDic.items() if len(v) >= (phaseNumber/2)}
    idDic = {}
    for i in dic:
        if i not in idDic:
            idDic[i] = []
        for j in dic[i]:
            ids = maxAbunBinId[j].split(',')
            for l in ids:
                idDic[i].append(l)
    outDic = {k:v for k,v in idDic.items() if len(v) >= phaseNumber}
    return outDic, dic

def PhaseScoreAnalysis(tuple):
    """function for multiple process

    Arguments:
        tuple {tuple} -- (mapDic, phase_length, k, phase_number, out_RNA, phaseScore_cutoff, phaseRatio_cutoff, island_length, type_)

    Returns:
        int -- 1
    """
    out_phasiRNA = TwoDepDic()
    out_allsiRNA = TwoDepDic()
    record_marker = 0
    mapDic = tuple[0]
    phase_length = tuple[1]
    k = tuple[2]
    phase_number = tuple[3]
    phaseScore_cutoff = tuple[4]
    phaseRatio_cutoff = tuple[5]
    island_length = tuple[6]
    type_ = tuple[7]
    start = mapDic['length'][0]
    end = mapDic['length'][1]
    bins = AlignBins(start, end + 1, phase_length)
    binDic = ParseBins(mapDic['attribute'], bins, phase_length)
    Vprint(f'{k} phaseScore analysis start...')
    maxAbunBin = FindMaxAbun(binDic['binAbun'])[0][0]
    try:
        condicatePhaseCluster = binDic['bin21Pos'][maxAbunBin].split(',')
        # if not TwoStrandMappingDetect(condicatePhaseCluster, binDic['bin21Strand']):
        #     Vprint(f'*WARNING, {k} has been skiped phaseScore analysis due to absence of two strand mapping')
        #     Vprint(f'{k} phaseScore analysis finished...')
        #     return 1

        # phaseNumber filter
        PhaseNumberFilter(condicatePhaseCluster, phase_number, k)
    except KeyError:
        Vprint(f'*WARNING, {k} has been skiped phaseScore analysis due to the absence of 21nt seq in maxAbun bin')
        Vprint(f'{k} phaseScore analysis finished...')
        return ({}, {})
    # distance filter
    maxAbunBinId = binDic['bin21id'][maxAbunBin]
    filteredCluster = DistanceFilter(condicatePhaseCluster, phase_number, maxAbunBinId, island_length)[0]
    unFilteredClusterPos = DistanceFilter(condicatePhaseCluster, phase_number, maxAbunBinId, island_length)[1]
    for i in filteredCluster:
        start = unFilteredClusterPos[i][0]
        end = unFilteredClusterPos[i][-1]
        phaseClusterBins = AlignBins(start, end + 1, phase_length)
        modifiedBinDic = ParseBins(mapDic['attribute'], phaseClusterBins, phase_length)
        maxAbunBin = FindMaxAbun(modifiedBinDic['binAbun'])[0][0]
        highestAbun = modifiedBinDic['binAbun'][maxAbunBin]
        # ! get the total abundance
        totalAbun = 0
        for j in modifiedBinDic['binAbun']:
            totalAbun += modifiedBinDic['binAbun'][j]
        phaseRatio = round(highestAbun/totalAbun, 3)
        try:
            phaseNumber = len(modifiedBinDic['bin21Pos'][maxAbunBin].split(','))
        except KeyError:
            continue
        if phaseNumber < phase_number:
            continue
        # get phaseAbun
        try:
            phaseAbun = modifiedBinDic['bin21Abun'][maxAbunBin]
        except KeyError:
            Vprint('skip this subCluster due to absence of 21nt seq in mostAbun bin')
            continue
        # conculate phaseScore
        phaseScore = round(phaseRatio*phaseNumber*math.log(phaseAbun), 3)
        # ! phaseScore and phaseRatio filter
        if phaseScore >= phaseScore_cutoff and phaseRatio >= phaseRatio_cutoff:
            # print(f'phaseScore\t{k}\t{start}\t{end}\t{phaseRatio}\t{phaseNumber}\t{phaseAbun}\t{phaseScore}')
            record_marker += 1
            PhaseScoresRNA(binDic, 'phase_length', start, end, record_marker, modifiedBinDic, maxAbunBin, phaseRatio, phaseNumber, phaseAbun, phaseScore, out_phasiRNA, k, type_)
            PhaseScoresRNA(binDic, 'all', start, end, record_marker, modifiedBinDic, maxAbunBin, phaseRatio, phaseNumber, phaseAbun, phaseScore, out_allsiRNA, k, type_)
    Vprint(f'{k} phaseScore analysis finished...')
    return (out_phasiRNA, out_allsiRNA)

def ParsePhaseScoreMap(file):
    """parse phase score map file

    Parameters
    ----------
    file : file

    Returns
    -------
    dict
    """
    dic = {}
    out_dic = {}
    with open(file, 'r') as fn:
        for line in fn:
            l = line.strip().split("\t")
            id = l[0]
            strand = l[1]
            gene = l[2]
            if strand == '+':
                pos = int(l[3])
            elif strand == '-':
                pos = int(l[3]) + 2
            seq = l[4]
            try:
                abun = float(id.split("_x")[-1])
            except ValueError:
                abun = float(id.split("@")[-1])
            key = str(pos) + "\t" + strand + "\t" + seq
            if gene not in dic:
                dic[gene] = {}
            if key not in dic[gene]:
                dic[gene][key] = str(abun) + "\t" + id
            else:
                abun = float(dic[gene][key].split("\t")[0]) + abun
                dic[gene][key] = str(abun) + "\t" + id
        for i in dic:
            if i not in out_dic:
                out_dic[i] = {}
            for j in dic[i]:
                tmp = j.split("\t")
                pos_ = int(tmp[0])
                strand_ = tmp[1]
                seq_ = tmp[2]
                tmp1 = dic[i][j].split("\t")
                abun_ = float(tmp1[0])
                id_ = tmp1[1]
                if 'attribute' not in out_dic[i]:
                    out_dic[i]['attribute'] = [(id_, strand_, pos_, seq_, abun_)]
                else:
                    out_dic[i]['attribute'].append((id_, strand_, pos_, seq_, abun_))
    return out_dic

def HypergeometricAnalysis(functionParament):
    """hypergeometric analysis

    Parameters
    ----------
    functionParament : tuple
        (mapdata, setting1, direct)
    mapdata : dict
        _description_
    setting1 : tuple
        setting1 = (windowLength, phase_length, phase_number, pvalue_cutoff, min_read_num, record_marker, type_)
    direct : str
        forward | reverse
        

    Returns
    -------
    tuple
        (out_phasiRNA, out_allsiRNA, out_phasiRNA_cluster, out_allsiRNA_cluster)
    """
    mapdata = functionParament[0]
    setting1 = functionParament[1]
    direct = functionParament[2]
    if direct.upper() == 'FORWARD':
        tmp_tuple_= FirstScaning(mapdata, setting1)
    elif direct.upper() == 'REVERSE':
        tmp_tuple_ = ReverseFirstScaning(mapdata, setting1)
    # final_cluster = tmp_tuple_[0]
    # final_strand = tmp_tuple_[1]
    # final_allsiRNA = tmp_tuple_[2]
    # phase_number = setting1[2]
    # phase_length = setting1[1]
    # pvalue_cutoff = setting1[3]
    # min_read_num = setting1[4]
    # setting2 = (final_cluster, final_strand, phase_number, phase_length, pvalue_cutoff, min_read_num)
    # tmp_tuple_ = SecondScaning(final_allsiRNA, setting2)

    return tmp_tuple_

def WritingData(outputfilename: str, dic_: dict, dic1_: dict, type_: str, gdna_anno='', duplication_ = False):
    """Writing data

    Parameters
    ----------
    outputfilename : str
        phasiRNA outputfile name or allsiRNA ouputfile name
    dic_ : dict
        output_phasiRNA or output_allsiRNA
    dic1_ : dict
        output_phasiRNA_cluster or output_allsiRNA_cluster
    type : dict
        cdna | gdna | Fl-cdna
    """

    o_ = open(outputfilename, 'a+')
    if outputfilename != '':
        o_.write('#' + ' ' + type_.upper() + ' based results' + "\n")
        for i in dic_:
            for j in dic_[i]:
                if gdna_anno != '':
                    if len(gdna_anno[i][j]) == 0:
                        tmp = ['>', 'cluster:', str(dic1_[i][j][0]), str(dic1_[i][j][-1]), "Intergenic"]
                        o_.write(' '.join(tmp) + "\n")
                    else:
                        if duplication_:
                            tmp = ['>', 'cluster:', str(dic1_[i][j][0]), str(dic1_[i][j][0] + 231), "\t".join(gdna_anno[i][j])]
                        else:
                            tmp = ['>', 'cluster:', str(dic1_[i][j][0]), str(dic1_[i][j][-1]), "\t".join(gdna_anno[i][j])]
                        o_.write(' '.join(tmp) + "\n")
                    for k in dic_[i][j]:
                        tmp = dic_[i][j][k]
                        out = tmp + "\t" + "\n"
                        o_.write(out)
                else:
                    if duplication_:
                        tmp = ['>', 'cluster:', str(dic1_[i][j][0]), str(dic1_[i][j][0] + 231), i]
                    else:
                        tmp = ['>', 'cluster:', str(dic1_[i][j][0]), str(dic1_[i][j][-1]), i]
                    o_.write(' '.join(tmp) + "\n")
                    for k in dic_[i][j]:
                        tmp = dic_[i][j][k]
                        out = tmp + "\t" + "\n"
                        o_.write(out)
    o_.close()

def WritingDuplication(output: str, dic_: dict, type_: str, gdna_anno = ''):
# discard function
    """duplication data writting

    Parameters
    ----------
    output : str
        outputfilename
    dic_ : dict
        out_duplication_phasiRNA, out_duplication_allsiRNA
    """
    o_ = open(output, 'a+')
    if output != '':
        o_.write('#' + ' ' + 'unduplication cluster' + ' based results ' + type_.upper() +  "\n")
        for i in dic_:
            for j in dic_[i]:
                if gdna_anno != '':
                    tmp = ['>' + 'cluster:', "\t".join(gdna_anno[i][j])]
                    o_.write(' '.join(tmp) + "\n")
                    for k in dic_[i][j]:
                        tmp = dic_[i][j][k]
                        out = tmp + "\t" + "\n"
                        o_.write(out)
                else:
                    tmp = ['>' + 'cluster:']
                    o_.write(' '.join(tmp) + "\n")
                    for k in dic_[i][j]:
                        tmp = dic_[i][j][k]
                        out = tmp + "\t" + "\n"
                        o_.write(out)
    o_.close()

def PrePhasiScoreAnalysis(mapDic: dict, ref: dict):
    """running before phaseScoreAnalysis

    Parameters
    ----------
    mapDic : dict
        _description_
    ref : dict
        refcDNA | refgDNA | refFlnc
    """
    for k in list(mapDic.keys()):
        if 'length' not in mapDic:
            start = ref[k][0]
            end = ref[k][1]
            mapDic[k]['length'] = [start, end]

def ReverseFirstScaning(mapdata: dict, setting: tuple):
    """first scaning

    Parameters
    ----------
    mapdata : dict
        _description_
    setting : tuple
        setting = (windowLength, phase_length, phase_number, pvalue_cutoff, min_read_num, type_)

    Returns
    -------
    tuple
        (final_cluster, final_strand, final_allsiRNA)
    """
    windowLength = setting[0]
    phase_length = setting[1]
    phase_number = setting[2]
    pvalue_cutoff = setting[3]
    min_read_num = setting[4]
    record_marker = 0
    type_ = setting[5]
    final_phasiRNA = nestedDic()
    final_allsiRNA = nestedDic()
    for k in mapdata:
        max_readn = 0
        pvalue = 0.0
        line = mapdata[k].split('\t')
        key = k.split('\t')
        sRNA_no = float(line[0])
        sRNA = line[1]
        strand = key[1]
        geneid = key[0]
        sRNA_pos = int(key[2])
        total_n = 0
        total_k = 0
        phasiRNA = []
        allsiRNA = []
        if strand == '+':
            for i in range(sRNA_pos - windowLength - 1, sRNA_pos + 1)[::-1]:
                query = geneid + '\t' + '+' + '\t' + str(i)
                if query in mapdata:
                    total_n += 1
                    allsiRNA.append(query + '\t' + mapdata[query])
            for i in range(sRNA_pos - windowLength - 2 - 1, sRNA_pos - 2 + 1)[::-1]:
                query = geneid + '\t' + '-' + '\t' + str(i)
                if query in mapdata:
                    total_n += 1
                    allsiRNA.append(query + '\t' + mapdata[query])
            for i in range(sRNA_pos - windowLength - 1, sRNA_pos + 1, phase_length)[::-1]:
                query = geneid + '\t' + '+' + '\t' + str(i)
                if query in mapdata and mapdata[query].split('\t')[-1] == str(phase_length):
                    total_k += 1
                    phasiRNA.append(query + '\t' + mapdata[query])
            for i in range(sRNA_pos - 2 - windowLength - 1, sRNA_pos - 2 + 1, phase_length)[::-1]:
                query = geneid + '\t' + '-' + '\t' + str(i)
                if query in mapdata and mapdata[query].split('\t')[-1] == str(phase_length):
                    total_k += 1
                    phasiRNA.append(query + '\t' + mapdata[query])
        elif strand == '-':
            for i in range(sRNA_pos - windowLength - 1, sRNA_pos + 1)[::-1]:
                query = geneid + '\t' + '-' + '\t' + str(i)
                if query in mapdata:
                    total_n += 1
                    allsiRNA.append(query + '\t' + mapdata[query])
            for i in range(sRNA_pos - windowLength + 2 - 1, sRNA_pos + 2 + 1)[::-1]:
                query = geneid + '\t' + '+' + '\t' + str(i)
                if query in mapdata:
                    total_n += 1
                    allsiRNA.append(query + '\t' + mapdata[query])
            for i in range(sRNA_pos - windowLength - 1, sRNA_pos + 1, phase_length)[::-1]:
                query = geneid + '\t' + '-' + '\t' + str(i)
                if query in mapdata and mapdata[query].split('\t')[-1] == str(phase_length):
                    total_k += 1
                    phasiRNA.append(query + '\t' + mapdata[query])
            for i in range(sRNA_pos - windowLength + 2 - 1, sRNA_pos + 2 + 1, phase_length)[::-1]:
                query = geneid + '\t' + '+' + '\t' + str(i)
                if query in mapdata and mapdata[query].split('\t')[-1] == str(phase_length):
                    total_k += 1
                    phasiRNA.append(query + '\t' + mapdata[query])

        # phase_number and pvalue cutoff
        if total_k >= phase_number:
            for i in range(total_k, phase_length):
                pvalue += float(binomial(windowLength * 2 - 1 - phase_length, total_n - i) * binomial(phase_length, i) / binomial(windowLength * 2 - 1, total_n))
            if pvalue <= pvalue_cutoff and pvalue != 0.0:
                for i in range(len(phasiRNA)):
                    value_list = phasiRNA[i].split('\t')
                    sRNA_no = float(value_list[3])
                    if sRNA_no > max_readn:
                        max_readn = sRNA_no

                if min_read_num == 0 or max_readn >= min_read_num:  # min_read_num cutoff
                    record_marker += 1
                    if type_ == 'cdna':
                        out_record_marker = 'CR'+str(record_marker)
                    elif type_ == 'gdna':
                        out_record_marker = 'GR'+str(record_marker)
                    elif type_ == 'flnc':
                        out_record_marker = "FR"+str(record_marker)
                    for i in range(len(phasiRNA)):
                        value_list = phasiRNA[i].split('\t')
                        geneid = value_list[0]
                        strand = value_list[1]
                        readpos = int(value_list[2])
                        query = geneid + '\t' + strand + '\t' + str(readpos)
                        final_phasiRNA[geneid][query][record_marker] = phasiRNA[i] + '\t' + str('%.3e' % pvalue) + '\t' + geneid+"_"+str(out_record_marker)
                        # final_cluster[geneid].append(readpos)
                        # final_strand[geneid][readpos] = strand
                    for i in range(len(allsiRNA)):
                        value_list = allsiRNA[i].split('\t')
                        geneid = value_list[0]
                        readpos = int(value_list[2])
                        strand = value_list[1]
                        query = geneid + '\t' + strand + '\t' + str(readpos)
                        final_allsiRNA[geneid][query][record_marker] = allsiRNA[i] + "\t" + str('%.3e' % pvalue) + '\t' + geneid+"_"+str(out_record_marker)
    return (final_phasiRNA, final_allsiRNA)

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

def PhaseNumberFilter(condicatePhaseCluster, phase_number, k):
    """phase number filter

    Parameters
    ----------
    condicatePhaseCluster : 
    phase_number : int
    k : geneid or chrid

    Returns
    -------
    int
    """
    if len(condicatePhaseCluster) < phase_number:
        Vprint(f'*WARNING, {k} has been skiped phaseScore analysis due to phaseNumber less than cutoff')
        Vprint(f'{k} phaseScore analysis finished...')
        return {}

def GeneratingSRNACluster(refgdna, gdna_map, phase_length, extended_maplen, outfile):
    """generating RNA cluster based on gdna data for phaseScore method, will generating a fa file

    Parameters
    ----------
    refgdna : file
    gdna_map : map file
    phase_length : int
    extended_maplen : int
    outfile : outfile
    """
    fo = open(outfile, 'w')
    dic = nestedDic()
    dic1 = nestedDic()
    for query in SeqIO.parse(refgdna, 'fasta'):
        id = query.id
        seq = str(query.seq)
        len_ = len(seq)
        dic[id] = len_
        dic1[id] = seq
    # Parse map simply
    map_dic = nestedDic()
    with open(gdna_map, 'r') as fn:
        map_count = 0
        for line in fn:
            line1 = line.strip()
            l = line.strip().split("\t")
            id = l[0]
            strand = l[1]
            chr_ = l[2]
            pos = int(l[3])
            seq = l[4]
            map_count += 1
            if strand == '-':
                pos += 2
            key = id + "\t" + strand + "\t" + str(map_count)
            map_dic[chr_][key] = pos
    for i in dic:
        loci = []
        for j in sorted(map_dic[i].items(), key = lambda kv:(kv[1], kv[0])):
            v = j[1]
            loci.append(v)
        DefineCluster(i, loci, phase_length, extended_maplen, dic, dic1, fo)
    fo.close()

def DefineCluster(chr, loci, phase_length, extended_maplen, chr_len, chr_seq, fo):
    """cluster define

    Parameters
    ----------
    chr : str
    loci : 
    phase_length : int
    extended_maplen : 
    chr_len : 
    chr_seq : 
    fo : 
    """
    loci_length = len(loci)
    cluster_num = 0
    index = []
    for i in range(0, loci_length):
        if i == 0:
            index.append(i)
        else:
            distance = loci[i] - loci[i-1]
            if distance <= 100:
                index.append(i)
            else:
                num = len(index)
                if num <= 4:
                    index = []
                    index.append(i)
                    continue
                start_read_pos = loci[index[0]]
                last_read_pos = loci[index[num - 1]]
                if last_read_pos - start_read_pos <= 3 * phase_length:
                    continue
                if start_read_pos > extended_maplen:
                    cluster_beg = start_read_pos - extended_maplen
                else:
                    cluster_beg = 1
                if last_read_pos + extended_maplen + 30 < chr_len[chr]:
                    cluster_end = last_read_pos + extended_maplen + 30
                else:
                    cluster_end = chr_len[chr]

                cluster_num += 1
                out_str = Extract_cluster_seq(chr, cluster_beg, cluster_end, cluster_num, chr_seq)
                fo.write(out_str+"\n")
                index = []
                index.append(i)

def Extract_cluster_seq(chr, beg, end, cluster_num, seq_dic):
    """generating cluster from defined cluster

    Parameters
    ----------
    chr : str
    beg : int
    end : int
    cluster_num : cluster id
    seq_dic : dict

    Returns
    -------
    str
        sequence
    """
    seq = seq_dic[chr][beg-1:end]
    out_str = '>'+chr+'_'+str(cluster_num)+"\t"+str(beg)+"\t"+str(end)+"\n"+seq
    return out_str

def BuildIndex(tmpfile):
    """build index

    Parameters
    ----------
    tmpfile : file

    Returns
    -------
    str
        index prefix
    """
    tuple_ = os.path.split(os.path.realpath(tmpfile))
    wd = tuple_[0]
    basename = tuple_[1]
    cmd = '[ -e ' + wd + '/.phasiHuter_bowtieIndex ] || mkdir ' + wd + '/.phasiHuter_bowtieIndex'
    os.system(cmd)
    cmd = f'[ -e {wd}/.phasiHuter_bowtieIndex/{basename}.1.ebwt ] || bowtie-build {tmpfile} {wd}/.phasiHuter_bowtieIndex/{basename} > /dev/null 2>&1'
    os.system(cmd)
    return wd + '/.phasiHuter_bowtieIndex/' + basename

def LoadPhaseScoreData(type: str, map, ref, phase_length, extended_maplen, tmp_file, max_hits, real_gdna_map, fa):    
    """loading phase Score data, ParallelPhaseScore first step

    Parameters
    ----------
    type : str
    map : 
    ref : file
    phase_length : int
    extended_maplen : int
    tmp_file : file
    max_hits : int
    real_gdna_map : file
    fa : file

    Returns
    -------
    dict
    """
    if type == 'cdna':
        dna_mapDic = ParsePhaseScoreMap(map)
        refDNA = ParseRef(ref, 'cdna', 'PhaseScore')
        PrePhasiScoreAnalysis(dna_mapDic, refDNA)
    elif type == 'flnc':
        dna_mapDic = ParsePhaseScoreMap(map)
        refDNA = ParseRef(ref, 'flnc', 'PhaseScore')
        PrePhasiScoreAnalysis(dna_mapDic, refDNA)
    elif type == 'gdna':
        GeneratingSRNACluster(ref, map, phase_length, extended_maplen, tmp_file)
        index = BuildIndex(tmp_file)
        cmd = f'bowtie -f -a -v 0 -m {max_hits} {index} {fa} {real_gdna_map} > /dev/null 2>&1'
        os.system(cmd)
        dna_mapDic = ParsePhaseScoreMap(real_gdna_map)
        refDNA = ParseRef(tmp_file, 'gdna', 'PhaseScore')
        PrePhasiScoreAnalysis(dna_mapDic, refDNA)
    return dna_mapDic

def ParallelHypergeometric(mapfile, phase_length, parallel_number, setting1, o, all):
    """parallel hypergeometric analysis

    Parameters
    ----------
    mapfile : 
    phase_length : int
    parallel_number : int
    setting1 : tupel
    type : 
    o : io
    all : io
    """
    cdna_mapdata = ParseHypergeometricMap(mapfile, phase_length)
    cdna_mapdata_dic = split_dic(cdna_mapdata)
    process_poll = ProcessPoolExecutor(max_workers=parallel_number)
    
    futures_F = []
    futures_R = [] 

    for i in cdna_mapdata_dic:
        forwardFunctionInput = (cdna_mapdata_dic[i], setting1, 'forward')
        future_F = process_poll.submit(HypergeometricAnalysis, forwardFunctionInput)
        futures_F.append(future_F)

    for i in cdna_mapdata_dic:
        reverseFunctionInput = (cdna_mapdata_dic[i], setting1, 'reverse')
        future_R = process_poll.submit(HypergeometricAnalysis, reverseFunctionInput)
        futures_R.append(future_R)

    wait(futures_F, return_when=ALL_COMPLETED)
    wait(futures_R, return_when=ALL_COMPLETED)

    process_poll.shutdown()
    for i in futures_F:
        future_F = i.result()
        HypergeometricOutput(future_F[0], o)
        HypergeometricOutput(future_F[1], all)

    for i in futures_R:
        future_R = i.result()
        HypergeometricOutput(future_R[0], o)
        HypergeometricOutput(future_R[1], all)

def HypergeometricOutput(final_sRNAs, fo):
    """writing hypergeometric output

    Parameters
    ----------
    final_sRNAs : dict
    fo : io
    """
    for gene in final_sRNAs:
        for i in final_sRNAs[gene]:
            for j in final_sRNAs[gene][i]:
                l = final_sRNAs[gene][i][j].split("\t")
                record = l[-1]
                outstring = "\t".join(l[:-1]) + "\t" + "-" + "\t" + "-" + "\t" + "-" + "\t" + "-" + "\t" + record + "\n"
                fo.write(outstring)

def ParallelPhaseScore(parallel_number, type_, mapfile, ref, phase_length, extended_maplen, tmpfile, max_hits, real_gdna_map, phase_number, o, fa, all, phaseScore_cutoff, phaseRatio_cutoff, island):
    """parallel run phaseScore
    分为两步, 第一步LoadPhaseScoreData, 第二步PhaseScoreAnalysis

    Parameters
    ----------
    parallel_number : int
    type_ : gdna | cdna
    mapfile : 
    ref : 
    phase_length : int
    extended_maplen : int
    tmpfile : file
    max_hits : int
    real_gdna_map : file
    phase_number : int
    o : io
    fa : file
    all : io
    """
    futures = []
    mapDic = LoadPhaseScoreData(type_, mapfile, ref, phase_length, extended_maplen, tmpfile, max_hits, real_gdna_map, fa)
    process_poll = ProcessPoolExecutor(max_workers = parallel_number)
    for k in mapDic:
        tuple = (mapDic[k], phase_length, k, phase_number, phaseScore_cutoff, phaseRatio_cutoff, island*phase_length, type_)
        future = process_poll.submit(PhaseScoreAnalysis, tuple)
        futures.append(future)

    wait(futures, return_when=ALL_COMPLETED)
    process_poll.shutdown()

    if type_ == 'gdna':
        dic = ParseTmpfile(tmpfile)
    elif type_ == 'cdna' or type_ == 'flnc':
        dic = {}

    for i in futures:
        future = i.result()
        phasiRNA = future[0]
        allsiRNA = future[1]
        phaseScoreOut(phasiRNA, o, type_, dic)
        phaseScoreOut(allsiRNA, all, type_, dic)


def phaseScoreOut(final_sRNAs, fo, type_, dic):
    """writing phaseScore output

    Parameters
    ----------
    final_sRNAs : dict
    fo : io
    """
    if type_ == 'cdna' or type_ == 'flnc':
        for gene in final_sRNAs:
            for record_marker in final_sRNAs[gene]:
                for i in final_sRNAs[gene][record_marker]:
                    l = i.strip().split("\t")
                    geneid, sRNAid, seq, seq_len, abun, pos, strand, marker, phaseRatio, phaseNumber, phaseAbun, phaseScore = l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10], l[11]
                    realign = [geneid, strand, pos, abun, sRNAid, seq, seq_len, '-', phaseRatio, phaseNumber, phaseAbun, phaseScore, marker]
                    fo.write("\t".join(realign)+"\n")
    elif type_ == 'gdna':
        for gene in final_sRNAs:
            for record_marker in final_sRNAs[gene]:
                for i in final_sRNAs[gene][record_marker]:
                    l = i.strip().split("\t")
                    geneid, sRNAid, seq, seq_len, abun, pos, strand, marker, phaseRatio, phaseNumber, phaseAbun, phaseScore = l[0], l[1], l[2], l[3], l[4], int(l[5]), l[6], l[7], l[8], l[9], l[10], l[11]
                    pos = pos + dic[geneid]
                    realign = [geneid, strand, str(pos-1), abun, sRNAid, seq, seq_len, '-', phaseRatio, phaseNumber, phaseAbun, phaseScore, marker]
                    fo.write("\t".join(realign)+"\n")


def PhaseScoresRNA(binDic, type_, start, end, record_marker, modifiedBinDic, maxAbunBin, phaseRatio, phaseNumber, phaseAbun, phaseScore, out_RNA, k, type__):
    """generating phaseScore small RNA detail information

    Parameters
    ----------
    binDic : dict
    type_ : phase_length | all
    start : int
    end : int
    record_marker : int
    modifiedBinDic : dict
    maxAbunBin : int
    phaseRatio : float
    phaseNumber : float
    phaseAbun : float
    phaseScore : float
    out_RNA : 
    k : geneid
    type__ : gdna | cdna
    """
    if type_ == 'phase_length':
        negative_strand = True
        positive_strand = True
        try:
            binDic['bin21Strand']['-']
        except KeyError:
            negative_strand = False
        try:
            binDic['bin21Strand']['+']
        except KeyError:
            positive_strand = False

        if positive_strand and negative_strand:
            if start in binDic['bin21Strand']['-']:
                start -= 2
            if end in binDic['bin21Strand']['-'] and end in binDic['bin21Strand']['+']:
                pass
            elif end in binDic['bin21Strand']['-']:
                end -= 2
        elif negative_strand and not positive_strand:
            if start in binDic['bin21Strand']['-']:
                start -= 2
            elif end in binDic['bin21Strand']['-']:
                end -= 2
        if type__ == 'gdna':
            out_record_marker = str(k) + "_" + "PG" + str(record_marker)
        elif type__ == 'cdna':
            out_record_marker = str(k) + "_" + "PC" + str(record_marker)
        elif type__ == 'flnc':
            out_record_marker = str(k) + "_" + "PF" + str(record_marker)
        # ! write phasiRNA to file
        loop_cluster = modifiedBinDic['bin21Pos'][maxAbunBin].split(',')
        loop_cluster = sorted([ int(ii) for ii in loop_cluster ])
        tmp_variable = [str(phaseRatio), str(phaseNumber), str(phaseAbun), str(phaseScore)]
        out_cluster = []
        if positive_strand and negative_strand:
            for i1 in [0]*len(loop_cluster):
                if loop_cluster.count(loop_cluster[i1]) == 2 and len(loop_cluster) >= 2:
                    loop_cluster[i1] = loop_cluster[i1] - 2
                    out_cluster.append(str(loop_cluster[i1]) + "\t" + '-')
                    loop_cluster.pop(i1)
                elif loop_cluster[i1] in binDic['bin21Strand']['-'] and loop_cluster[i1] in binDic['bin21Strand']['+']:
                    out_cluster.append(str(loop_cluster[i1]) + '\t' + '+')
                    loop_cluster.pop(i1)
                elif loop_cluster[i1] in binDic['bin21Strand']['-']:
                    out_cluster.append(str(loop_cluster[i1] - 2) + '\t' + '-')
                    loop_cluster.pop(i1)
                elif loop_cluster[i1] in binDic['bin21Strand']['+']:
                    out_cluster.append(str(loop_cluster[i1]) + '\t' + '+')
                    loop_cluster.pop(i1)
            for i1 in out_cluster:
                outstring = k + "\t" + "\t".join(binDic['bin21out'][i1]) + "\t" + str(i1) + "\t" + out_record_marker + "\t" +  "\t".join(tmp_variable) + "\n"
                out_RNA[k][record_marker].append(outstring)
        elif negative_strand and not positive_strand:
            for i1 in [0]*len(loop_cluster):
                if loop_cluster.count(loop_cluster[i1]) == 2 and len(loop_cluster) >= 2:
                    loop_cluster[i1] = loop_cluster[i1] - 2
                    out_cluster.append(str(loop_cluster[i1]) + "\t" + '-')
                    loop_cluster.pop(i1)
                elif loop_cluster[i1] in binDic['bin21Strand']['-']:
                    out_cluster.append(str(loop_cluster[i1] - 2) + '\t' + '-')
                    loop_cluster.pop(i1)
            for i1 in out_cluster:
                outstring = k + "\t" + "\t".join(binDic['bin21out'][i1]) + "\t" + str(i1) + "\t" + out_record_marker + "\t" +  "\t".join(tmp_variable) +"\n"
                out_RNA[k][record_marker].append(outstring)
        elif positive_strand and not negative_strand:
            for i1 in [0]*len(loop_cluster):
                if loop_cluster.count(loop_cluster[i1]) == 2 and len(loop_cluster) >= 2:
                    loop_cluster[i1] = loop_cluster[i1] - 2
                    out_cluster.append(str(loop_cluster[i1]) + "\t" + '-')
                    loop_cluster.pop(i1)
                elif loop_cluster[i1] in binDic['bin21Strand']['+']:
                    out_cluster.append(str(loop_cluster[i1]) + '\t' + '+')
                    loop_cluster.pop(i1)
            for i1 in out_cluster:
                outstring = k + "\t" + "\t".join(binDic['bin21out'][i1]) + "\t" + str(i1) + "\t" + out_record_marker + "\t" +  "\t".join(tmp_variable) +"\n"
                out_RNA[k][record_marker].append(outstring)
    elif type_ == 'all':
        if type__ == 'gdna':
            out_record_marker = str(k) + "_" + "PG" + str(record_marker)
        elif type__ == 'cdna':
            out_record_marker = str(k) + "_" + "PC" + str(record_marker)
        elif type__ == 'flnc':
            out_record_marker = str(k) + "_" + "PF" + str(record_marker)
        for length in binDic['binStrand']:
            negative_strand = True
            positive_strand = True
            try:
                binDic['binStrand'][length]['-']
            except KeyError:
                negative_strand = False
            try:
                binDic['binStrand'][length]['+']
            except KeyError:
                positive_strand = False

            if positive_strand and negative_strand:
                if start in binDic['binStrand'][length]['-']:
                    start -= 2
                if end in binDic['binStrand'][length]['-'] and end in binDic['binStrand'][length]['+']:
                    pass
                elif end in binDic['binStrand'][length]['-']:
                    end -= 2
            elif negative_strand and not positive_strand:
                if start in binDic['binStrand'][length]['-']:
                    start -= 2
                elif end in binDic['binStrand'][length]['-']:
                    end -= 2
            # ! write phasiRNA to file
            for binid in modifiedBinDic['binlenpos'][length]:
                loop_cluster = modifiedBinDic['binlenpos'][length][binid]
                loop_cluster = sorted([ int(ii) for ii in loop_cluster ])
                tmp_variable = [str(phaseRatio), str(phaseNumber), str(phaseAbun), str(phaseScore)]
                out_cluster = []
                if positive_strand and negative_strand:
                    for i1 in [0]*len(loop_cluster):
                        if loop_cluster.count(loop_cluster[i1]) == 2 and len(loop_cluster) >= 2:
                            loop_cluster[i1] = loop_cluster[i1] - 2
                            out_cluster.append(str(loop_cluster[i1]) + "\t" + '-')
                            loop_cluster.pop(i1)
                        elif loop_cluster[i1] in binDic['binStrand'][length]['-'] and loop_cluster[i1] in binDic['binStrand'][length]['+']:
                            out_cluster.append(str(loop_cluster[i1]) + '\t' + '+')
                            loop_cluster.pop(i1)
                        elif loop_cluster[i1] in binDic['binStrand'][length]['-']:
                            out_cluster.append(str(loop_cluster[i1] - 2) + '\t' + '-')
                            loop_cluster.pop(i1)
                        elif loop_cluster[i1] in binDic['binStrand'][length]['+']:
                            out_cluster.append(str(loop_cluster[i1]) + '\t' + '+')
                            loop_cluster.pop(i1)
                    for i1 in out_cluster:
                        outstring = k + "\t" + "\t".join(binDic['binout'][length][i1]) + "\t" + str(i1) + "\t" + out_record_marker + "\t" +  "\t".join(tmp_variable) + "\n"
                        out_RNA[k][record_marker].append(outstring)
                elif negative_strand and not positive_strand:
                    for i1 in [0]*len(loop_cluster):
                        if loop_cluster.count(loop_cluster[i1]) == 2 and len(loop_cluster) >= 2:
                            loop_cluster[i1] = loop_cluster[i1] - 2
                            out_cluster.append(str(loop_cluster[i1]) + "\t" + '-')
                            loop_cluster.pop(i1)
                        elif loop_cluster[i1] in binDic['binStrand'][length]['-']:
                            out_cluster.append(str(loop_cluster[i1] - 2) + '\t' + '-')
                            loop_cluster.pop(i1)
                    for i1 in out_cluster:
                        outstring = k + "\t" + "\t".join(binDic['binout'][length][i1]) + "\t" + str(i1) + "\t" + out_record_marker + "\t" +  "\t".join(tmp_variable) +"\n"
                        out_RNA[k][record_marker].append(outstring)
                elif positive_strand and not negative_strand:
                    for i1 in [0]*len(loop_cluster):
                        if loop_cluster.count(loop_cluster[i1]) == 2 and len(loop_cluster) >= 2:
                            loop_cluster[i1] = loop_cluster[i1] - 2
                            out_cluster.append(str(loop_cluster[i1]) + "\t" + '-')
                            loop_cluster.pop(i1)
                        elif loop_cluster[i1] in binDic['binStrand'][length]['+']:
                            out_cluster.append(str(loop_cluster[i1]) + '\t' + '+')
                            loop_cluster.pop(i1)
                    for i1 in out_cluster:
                        outstring = k + "\t" + "\t".join(binDic['binout'][length][i1]) + "\t" + str(i1) + "\t" + out_record_marker + "\t" +  "\t".join(tmp_variable) +"\n"
                        out_RNA[k][record_marker].append(outstring)

def ParseTmpfile(tmpfile):
    """parse phasescore gdna tmpfile sequence file

    Parameters
    ----------
    tmpfile : 
    """
    dic = nestedDic()
    for query in SeqIO.parse(tmpfile, 'fasta'):
        name = query.description.split("\t")
        id_ = name[0]
        start = int(name[1])
        dic[id_] = start
    return dic