import sys  
import re
from Bio import SeqIO
from collections import defaultdict

def nestedDic():
    """generating nested dictory

    Returns
    -------
    dict
        nested dictory
    """    
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
    return defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list)))))

def nestedDic():
    """generating nested dictory

    Returns
    -------
    dict
        nested dictory
    """    
    return defaultdict(nestedDic)

def readStdin():
    """read stdin and return a list

    Returns
    -------
    list
        \n striped list
    """
    fo = sys.stdin.readlines()
    fo = [i.strip() for i in fo]
    return fo

def reSplit(pattern, string):
    """return resplit list

    Parameters
    ----------
    pattern : 
        reg pattern
    string : 

    Returns
    -------
    list
        splited list
    """
    l = re.split(pattern, string)
    return l 

def convertFormat(fo):
    # chr_start:end to chr start end
    for i in fo:
        l = i.split('_')
        Chr = "_".join(l[:-1])
        coor = l[-1]
        coor1 = reSplit('\s+', coor)[0]
        coor2 = coor1.replace(':', "\t")
        print(f'{Chr}\t{coor2}')

def convertForIGVVisulization(fo):
    # chr start end to chr:start-end
    for i in fo:
        l = i.split("\t")
        Chr = l[0]
        start = l[1]
        end = l[2]
        print(f'{Chr}:{start}-{end}')


def convertToMultiplePhaseTankInputFormat(fo):
    # a.fa\nb.fa to a.fa,b.fa
    out = ''
    for i in fo:
        out = out + i + ','
    out = out[:-1] 
    print(out)



def main():
    import sys
    method = 1
    help = '''
    usage:
    ls * | phasiHunter read [n]
        option:
            # necessary options:
    
            # options with default value
            [0]: format combien category result
            [1]: chr_start:end to chr start end
            [2]: chr start end to chr:start-end
            [3]: a.fa\\nb.fa to a.fa,b.fa
            [4]: parse bedtools -a -b -wao | cut -f 1,2,3,7,8,9 stdout
            [5]: extract feature
            [6]: unitas output format, location:XR_004857687.1\\tstart:2096\\tend:2244 to XR_004857687\\t2096\\t2244
            [7]: XR_002261863.2\\t(1929:2055) to XR_002261863.2\\t1929\\t2055
            [8]: 添加chr\\tstart\\tend 在每一行的前面，1功能的增强版本
            [9]: NC_050096.1     34607   35318   exon    XM_020544715.3-1 **TO** XM_020544715.3  1       711; 此功能被排除在pipline了
            [10]: XM_008683071.3  (1422:1508)     (1332:1521) **TO** XM_008683071.3  1422    1508
            [11]: [transcript_id=XM_020544715.3]\\t[location=join(34607..35318,36037..36174,36259..36504,36600..36713,36822..37004,37416..37633,38021..39618,39701..40208)]\\t[gbkey=mRNA] to XM_020544715.3\\t1\\t711
            [12]: filter bed line with start > end
            [13]: filter XRR_1234--1 to XRR_1234
            [14]: extractForphasiRNA
            [15]: extractForPHASLoci
            [16]: phasiRNAnalyzer filter
            [17]: tissueAnalysis
            [18]: tissueFormat
    
            # optional options
    
            # other
            -v:       --  print version information
            -h:       --  print help information
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

    method = int(sys.argv[1])

    fo = readStdin()
    rDic = {
        0: catFormat,
        1 : convertFormat,
        2 : convertForIGVVisulization,
        3 : convertToMultiplePhaseTankInputFormat,
        4 : FormatBedtoolsAnnotate,
        5 : extractFeature,
        6 : unitasOutputFormat,
        7 : phasiHunterHypergeometricCDNAFormat,
        8 : overlapIntergenic,
        9 : generatingExonRegion,
        10: cdnaFormat,
        11: exonInfo,
        12: filterStartLargerThanEnd,
        13: caculategDNAClusterNumber,
        14: extractForphasiRNA,
        15: extractForPHASLoci,
        16: phasiRNAnalyzer_filter,
        17: tissueAnalysis,
        18: tissueFormat,
            }

    if method in rDic:
        rDic[method](fo)

def catFormat(fo):
    tag = 'y'
    for i in sorted(fo[:1]):
        l = i.split("\t")
        Category,Small_RNA,Target_gene,sRNA_loc,Deg_loc,Deg_count,sRNA_seq,Shift,Gene_annotation = l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8]
        if Category == 'Category' and tag == 'y':
            print(i)
            tag = 'n'
            continue
        if Category == 'Category' and tag == 'n':
            continue
        print(i)
    for i in sorted(fo[1:]):
        l = i.split("\t")
        Category,Small_RNA,Target_gene,sRNA_loc,Deg_loc,Deg_count,sRNA_seq,Shift,Gene_annotation = l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8]
        if Category == 'Category':
            continue
        print(i)

def phasiHunterHypergeometricCDNAFormat(fo):
    for i in fo:
        l = i.split("\t")
        Chr = l[0]
        l1 = l[1].split(":")
        start = l1[0].replace('(', '')
        end = l1[1].replace(')', '')
        print(f'{Chr}\t{start}\t{end}')

def FormatBedtoolsAnnotate(fo):
    dic = OneDepDic()
    for i in fo:
        l = i.split('\t')
        candidate = f'{l[0]}\t{l[1]}\t{l[2]}'
        anno = f'{l[3]}\t{l[4]}\t{l[5]}'
        dic[candidate].append(anno)

    for i in dic:
        print(i, end='')
        for j in dic[i]:
            print(f'\t{j}', end='')
        print()

def extractFeature(fo):
    dic = OneDepDic()
    for i in fo:
        l = i.split('\t')
        length = len(l)
        group = int(length/3)
        for i in range(0, group):
            if i == 0: 
                candidate = f'{l[0 + i*3]}\t{l[ 1 + i*3]}\t{l[2 + i*3]}'
            else:
                anno = f'{l[0 + i*3]}\t{l[ 1 + i*3]}\t{l[2 + i*3]}'
                dic[candidate].append(anno)

    for i in dic:
        print(i, end='')
        for j in dic[i]:
            j1 = j.split('\t')
            print(f'\t{j1[0]}', end='')
        print()

def unitasOutputFormat(fo):
    for i in fo:
        ele = i.split("\t")
        Chr = ele[0].split(':')[1]
        start = ele[1].split(":")[1]
        end = ele[2].split(":")[1]
        print(f'{Chr}\t{start}\t{end}')

def overlapIntergenic(fo):
    # chr_start:end to chr start end
    for i in fo:
        l = i.split('_')
        Chr = "_".join(l[:-1])
        coor = l[-1]
        coor1 = reSplit('\s+', coor)[0]
        coor2 = coor1.replace(':', "\t")
        print(f'{Chr}\t{coor2}\t{i}')

def generatingExonRegion(fo):
    former_rna = ''
    for i in fo:
        l = i.strip().split("\t")
        chr_ = l[0]
        start = int(l[1])
        end = int(l[2])
        try:
            j = l[4].split('-')
            rna = j[0]
        except IndexError:
            rna = j
        if rna != former_rna: # 进入第一个RNA
            o_start = start
            print(f'{rna}\t1\t{end - start}')
            former_rna = rna
            former_end = end - start
        else:
            if end - start + former_end > former_end + 1:
                print(f'{rna}\t{former_end + 1}\t{end - start + former_end}')
            former_end = end - start + former_end


def cdnaFormat(fo):
    for i in fo:
        l = i.strip().split("\t")
        rna = l[0]
        hcoor = l[1].split(";")
        pcoor = l[2].split(";")
        if hcoor[0] != '-':
            for coor in hcoor:
                coora = coor.replace('(', '')
                coorb = coora.replace(')', '')
                coorc = coorb.split(':')
                coor1 = coorc[0]
                coor2 = coorc[1]
                print(f'{rna}\t{coor1}\t{coor2}')
        if pcoor[0] != '-':
            for coor in pcoor:
                coora = coor.replace('(', '')
                coorb = coora.replace(')', '')
                coorc = coorb.split(':')
                coor1 = coorc[0]
                coor2 = coorc[1]
                print(f'{rna}\t{coor1}\t{coor2}')

def exonInfo(fo):
    for i in fo:
        l = i.strip().split("\t")
        rna_, loc_, feature_ = l[0], l[1], l[2]
        rna = rna_.replace('[', '').replace(']', '').split('=')[1]
        loc_1 = re.search("\d+[^)\]]*", loc_).group()
        locs = loc_1.split(',')
        feature = feature_.replace('[', '').replace(']', '').split('=')[1]
        former_start = 0
        former_end = 0
        count = 1
        for i in locs:
            region = i.split('..')
            try:
                start = int(region[0])
            except ValueError:
                start = int(region[0].replace('>', '').replace('<', ''))

            try:
                end = int(region[1])
            except ValueError:
                end = int(region[1].replace('>', '').replace('<', ''))
            except IndexError:
                end = int(region[0])

            # print(f'end:{end}')
            # print(f'start:{start}')
            if count == 1:
                former_start = start
                print(f'{rna}\t{start - former_start + 1}\t{end - start}')
                count += 1
                former_end = end - start
            else:
                print(f'{rna}\t{former_end+1}\t{end - start + former_end}')
                former_end = end - start + former_end
                former_start = start

def filterStartLargerThanEnd(fo):
    for i in fo:
        l = i.strip().split('\t')
        chr_ = l[0]
        start = int(l[1])
        end = int(l[2])
        if end > start:
            print(f'{chr_}\t{start}\t{end}')

def caculategDNAClusterNumber(fo):
    for i in fo:
        l = i.strip().split('\t')
        rna = l[0].split('--')[0]
        print(f'{rna}')



def extractForphasiRNA(fo):
    for i in fo:
        l = i.strip().split('\t')
        cluster_number, feature, geneid, method, ref, htc, hgc, hfc, ptc, pgc, pfc, htp, hgp, hfp, ptr, pgr, pfr, pts, pgs, pfs = \
        l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10], l[11], l[12], l[13], l[14], l[15], l[16], l[17], l[18], l[19]
        if cluster_number == 'cluster_number':
            continue
        else:
            if 'C' in ref:
                if 'H' in method:
                    print(f'{geneid} {feature} HC')
                elif 'P' in method:
                    print(f'{geneid} {feature} PC')
            elif 'F' in ref:
                if 'H' in method:
                    print(f'{geneid} {feature} HF')
                elif 'P' in method:
                    print(f'{geneid} {feature} PF')
            elif 'G' in ref:
                if 'H' in method:
                    print(f'{geneid} {feature} HG')
                elif 'P' in method:
                    print(f'{geneid} {feature} PG')

def extractForPHASLoci(fo):
    for i in fo:
        l = i.strip().split('\t')
        cluster_number, feature, geneid, method, ref, htc, hgc, hfc, ptc, pgc, pfc, htp, hgp, hfp, ptr, pgr, pfr, pts, pgs, pfs = \
        l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10], l[11], l[12], l[13], l[14], l[15], l[16], l[17], l[18], l[19]
        if cluster_number == 'cluster_number':
            continue
        else:
            if 'C' in ref:
                if 'H' in method:
                    for coor in htc.split(','):
                        if coor != '-':
                            tmp = coor.replace('(', '').replace(')', '').split(':')
                            start = tmp[0]
                            end = tmp[1]
                            if feature == 'Intergenic':
                                print(f'{geneid} Other:{feature} HC {start} {end}')
                            else:
                                print(f'{geneid} {feature}:{geneid} HC {start} {end}')
                elif 'P' in method:
                    for coor in ptc.split(','):
                        if coor != '-':
                            tmp = coor.replace('(', '').replace(')', '').split(':')
                            start = tmp[0]
                            end = tmp[1]
                            if feature == 'Intergenic':
                                print(f'{geneid} Other:{feature} PC {start} {end}')
                            else:
                                print(f'{geneid} {feature}:{geneid} PC {start} {end}')
            elif 'F' in ref:
                if 'H' in method:
                    for coor in hfc.split(','):
                        if coor != '-':
                            tmp = coor.replace('(', '').replace(')', '').split(':')
                            start = tmp[0]
                            end = tmp[1]
                            
                            if feature == 'Intergenic':
                                print(f'{geneid} Other:{feature} HF {start} {end}')
                            else:
                                print(f'{geneid} {feature}:{geneid} HF {start} {end}')
                elif 'P' in method:
                    for coor in pfc.split(','):
                        if coor != '-':
                            tmp = coor.replace('(', '').replace(')', '').split(':')
                            start = tmp[0]
                            end = tmp[1]
                            
                            if feature == 'Intergenic':
                                print(f'{geneid} Other:{feature} PF {start} {end}')
                            else:
                                print(f'{geneid} {feature}:{geneid} PF {start} {end}')
            elif 'G' in ref:
                if 'H' in method:
                    for coor in hgc.split(','):
                        if coor != '-':
                            tmp = coor.replace('(', '').replace(')', '').split(':')
                            start = tmp[0]
                            end = tmp[1]
                            
                            if feature == 'Intergenic':
                                print(f'{geneid} Other:{feature} HG {start} {end}')
                            else:
                                print(f'{geneid} {feature}:{geneid} HG {start} {end}')
                elif 'P' in method:
                    for coor in pgc.split(','):
                        if coor != '-':
                            tmp = coor.replace('(', '').replace(')', '').split(':')
                            start = tmp[0]
                            end = tmp[1]
                            
                            if feature == 'Intergenic':
                                print(f'{geneid} Other:{feature} PG {start} {end}')
                            else:
                                print(f'{geneid} {feature}:{geneid} PG {start} {end}')


def phasiRNAnalyzer_filter(fo):
    dic = OneDepDic()
    for line in fo:
        line1 = line.strip()
        l = line.strip().split("\t")
        id, strand, start, abun, sran, seq, anno, pvalue, record = l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8]
        length = len(seq)
        if length == 21:
            dic[record].append(line1)
        else:
            continue
    for record in dic:
        for line1 in dic[record]:
            if len(dic[record]) >= 4:
                print(line1)

def tissueAnalysis(fo):
    dic = OneDepDic()
    total_abun = 0
    total_count = 0
    for line in fo:
        line1 = line.strip()
        l = line.strip().split("\t")
        abun = float(l[1])
        total_abun += abun
        total_count += 1
    print(f'Total phasiRNA number: {total_count}')
    print(f'Total phasiRNA abun: {total_abun}')


def tissueFormat(fo):
    dic = OneDepDic()
    for line in fo:
        line1 = line.strip()
        l = line.strip().split("\t")
        lib = l[2]
        id = l[0]
        seq = l[1]
        abun = id.split('__')[3]
        print(f'{lib}\t{seq}\t{abun}')

if __name__ == "__main__":
    main()