# %%
import os
import sys
from sympy import binomial
from concurrent.futures import ProcessPoolExecutor, wait, ALL_COMPLETED
from collections import defaultdict
from phaseFunction import (
    GetBaseName,
    GetCurTime,
    Vprint,
    ParseRef,
    ParseGff3,
    ParseHypergeometricMap,
    nestedDic,
    AnnoGdna,
    PhaseScoreAnalysis,
    ParsePhaseScoreMap,
    HypergeometricAnalysis,
    WritingData,
    PrePhasiScoreAnalysis,
    split_dic,
    GeneratingSRNACluster,
    BuildIndex,
    LoadPhaseScoreData,
    ParallelHypergeometric,
    ParallelPhaseScore,
)

if __name__ == '__main__':
    # ----------> default variables <------------
    ctime = GetCurTime()
    programName = GetBaseName(__file__)
    Vprint(f'start run {programName}', enable=True)
    island_number = 5
    phase_length = 21
    method='b'
    phase_number = 4
    phaseScore_cutoff = 15
    phaseRatio_cutoff = 0.4
    max_hits = 10
    delete_phasiHuter_bowtieIndex = 'y'
    parallel_number = 1
    extended_maplen = 80
    pvalue_cutoff = 0.0005
    min_read_num = 0
    cdna_mapfile = ''
    gdna_mapfile = ''
    flnc_mapfile = ''
    cdna = ''
    gdna = ''
    flnc = ''
    tmpoutputfilename = 'phase_o.txt'
    tmpalloutputfilename = 'phase_a.txt'
    # enable_cdna = True
    # PHASgene_name = []
    # hypergeometric_result = nestedDic()
    # duplication = nestedDic()

    help = '''
    phase usage:
        option:
            # necessary options:
            -cm: file  --  map file based on reference transcriptome sequence
            -c:  file  --  reference transcritome sequence, fasta file
            -a:  out   --  allsiRNA cluster output file
            -o:  out   --  phasiRNA cluster output file
            -fa: file  --  sRNA file

            # options with default value
            -me: str   --  phasiRNA prediction method, h(hypergeometric test) | p(phase score) | b (both), default=b
            -il: int   --  phasiRNA cluster island, default=5
            -pl: int   --  phase length, 21 | 24, default=21
            -pn: int   --  phase number, default=4
            -mh: int   --  max hits when mapping to ref sequence, default=10
            -j:  int   --  parallel number, default=1
            -pv: float --  pvalue cutoff, default=0.001, only function with h method applied
            -ps: float --  phase score cutoff, default=15, only function with p method applied
            -pr: float --  phase ratio cutoff, default=0.4, only function with p method applied

            # optional options
            -f:  file  --  full length transcriptome sequence, fasta file
            -g:  file  --  reference genome sequence, fasta file
            -gm: file  --  map file based on reference genome sequence
            -fm: file  --  map file based on full length transcriptome sequence

            # other
            -cl: str  --  delete .phasiHuter_bowtieIndex, y|n, default=y
            -v:       --  print version information
            -h:       --  print help information

    '''
    version = '''
    version v1.0
    '''
    for i in range(1, len(sys.argv)):
        if sys.argv[i] == '-cm':
            cdna_mapfile = sys.argv[i+1]
        elif sys.argv[i] == '-il':
            island_number = int(sys.argv[i+1])
        elif sys.argv[i] == '-pl':
            phase_length = int(sys.argv[i+1])
        elif sys.argv[i] == '-pn':
            phase_number = int(sys.argv[i+1])
        elif sys.argv[i] == '-me':
            method = sys.argv[i+1]
        elif sys.argv[i] == '-ps':
            phaseScore_cutoff = float(sys.argv[i+1])
        elif sys.argv[i] == '-pr':
            phaseRatio_cutoff = float(sys.argv[i+1])
        elif sys.argv[i] == '-mh':
            max_hits = int(sys.argv[i+1])
        elif sys.argv[i] == '-j':
            parallel_number = int(sys.argv[i+1])
        elif sys.argv[i] == '-pv':
            pvalue_cutoff = float(sys.argv[i+1])
        elif sys.argv[i] == '-mn':
            min_read_num = int(sys.argv[i+1])
        elif sys.argv[i] == '-c':
            cdna = sys.argv[i+1]
        elif sys.argv[i] == '-cl':
            delete_phasiHuter_bowtieIndex = sys.argv[i+1]
        elif sys.argv[i] == '-g':
            gdna = sys.argv[i+1]
        elif sys.argv[i] == '-f':
            flnc = sys.argv[i+1]
        elif sys.argv[i] == '-gm':
            gdna_mapfile = sys.argv[i+1]
        elif sys.argv[i] == '-fm':
            flnc_mapfile = sys.argv[i+1]
        elif sys.argv[i] == '-fa':
            fa = sys.argv[i+1]
        elif sys.argv[i] == '-o':
            tmpoutputfilename = sys.argv[i+1]
        elif sys.argv[i] == '-a':
            tmpalloutputfilename = sys.argv[i+1]
        elif sys.argv[i] == '-v':
            print(version)
            sys.exit()
        elif sys.argv[i] == '-h':
            print(help)
            sys.exit()

    windowLength = phase_length * 11
    leastPhaseClusterLength = phase_length * island_number
    setting1 = (windowLength, phase_length, phase_number, pvalue_cutoff, min_read_num, 'cdna')
    setting2 = (windowLength, phase_length, phase_number, pvalue_cutoff, min_read_num, 'gdna')
    setting3 = (windowLength, phase_length, phase_number, pvalue_cutoff, min_read_num, 'flnc')
    TMP_WD = os.getcwd()
    tmp_file = TMP_WD + '/' + tmpalloutputfilename + '_tmpfile.fa'
    real_gdna_map = TMP_WD + '/' + tmpalloutputfilename + '_real_gdna_map.map'
    o_ = open(tmpoutputfilename, 'w+')
    all_ = open(tmpalloutputfilename, 'w+')
    # // Hypergeometric 
    if method == 'b' or method == 'h':
        if cdna_mapfile != '':
            ParallelHypergeometric(cdna_mapfile, phase_length, parallel_number, setting1, o_, all_)
        if gdna_mapfile != '':
            ParallelHypergeometric(gdna_mapfile, phase_length, parallel_number, setting2, o_, all_)
        if flnc_mapfile != '':
            ParallelHypergeometric(flnc_mapfile, phase_length, parallel_number, setting3, o_, all_)
    if method == 'b' or method == 'p':
        # // phaseScore
        # ! loadding data
        # ! phaseScore analysis
        if cdna_mapfile != '' and cdna != '':
            ParallelPhaseScore(parallel_number, 'cdna', cdna_mapfile, cdna, phase_length, extended_maplen, tmp_file, max_hits, real_gdna_map, phase_number, o_, fa, all_, phaseScore_cutoff, phaseRatio_cutoff, island_number)
        if gdna_mapfile != '' and gdna != '':
            ParallelPhaseScore(parallel_number, 'gdna', gdna_mapfile, gdna, phase_length, extended_maplen, tmp_file, max_hits, real_gdna_map, phase_number, o_, fa, all_, phaseScore_cutoff, phaseRatio_cutoff, island_number)
        if flnc_mapfile != '' and flnc != '':
            ParallelPhaseScore(parallel_number, 'flnc', flnc_mapfile, flnc, phase_length, extended_maplen, tmp_file, max_hits, real_gdna_map, phase_number, o_, fa, all_, phaseScore_cutoff, phaseRatio_cutoff, island_number)
    
    if gdna_mapfile != '' and gdna != '':
        if method == 'p' or method == 'b':
            os.system('rm ' +  tmp_file)
            os.system('rm ' + real_gdna_map)

    o_.close()
    all_.close()

    Vprint('Analysis finished', enable=True)