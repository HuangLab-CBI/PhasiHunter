# %%
from integrationFunction import *
import sys

def main1():
    Vprint('start intergrate data ...', enable=True)
    # ----------> variable declaration <------------
    phase_number = 4
    phase_length = 21
    pvalue_cutoff = 0.001
    min_read_num = 0
    parallel_number = 1
    island_number = 5
    flnc_anno = ''
    gdna_enable = 'n'
    help = '''
    phase usage:
        option:
            # necessary options:
            -io: file  --  phase module -o output file
            -ia: file  --  phase module -a output file
            -an: file  --  reference genome gff3 file
            -o:  out  --  integration phasiRNA cluster
            -a:  out  --  integration all siRNA cluster
            -s:  out  --  integration summary
            -g:  str  --  y | n, whether exist gdna based PHAS Loci

            # options with default value
            -j:  int   --  parallel number, default=1
            -pn: int   --  phase number, default=4
            -pl: int   --  phase length, 21 | 24, default=21
            -pv: float --  pvalue cutoff, default=0.001
            -mn: int   --  min read num, default=0  
            -il: int   --  phasiRNA cluster island, default=5

            # optional options
            -fn: file  --  full length transcript annotation file

            # other
            -v:       --  print version information
            -h:       --  print help information

    '''
    version = 'v1.0'

    for i in range(1, len(sys.argv)):
        if sys.argv[i] == '-io':
            phasiRNAfile = sys.argv[i+1]
        elif sys.argv[i] == '-ia':
            allsiRNAfile = sys.argv[i+1]
        elif sys.argv[i] == '-fn':
            flnc_anno = sys.argv[i+1]
        elif sys.argv[i] == '-il':
            island_number = int(sys.argv[i+1])
        elif sys.argv[i] == '-an':
            gff3 = sys.argv[i+1]
        elif sys.argv[i] == '-g':
            gdna_enable = sys.argv[i+1]
        elif sys.argv[i] == '-o':
            outfile = sys.argv[i+1]
        elif sys.argv[i] == '-a':
            allfile = sys.argv[i+1]
        elif sys.argv[i] == '-s':
            intergrationfile = sys.argv[i+1]
        elif sys.argv[i] == '-pn':
            phase_number = int(sys.argv[i+1])
        elif sys.argv[i] == '-pl':
            phase_length = int(sys.argv[i+1])
        elif sys.argv[i] == '-pv':
            pvalue_cutoff = float(sys.argv[i+1])
        elif sys.argv[i] == '-mn':
            min_read_num = float(sys.argv[i+1])
        elif sys.argv[i] == '-v':
            print(version)
            sys.exit()
        elif sys.argv[i] == '-h':
            print(help)
            sys.exit()

    TMP_WD = os.getcwd()
    ctime = GetCurTime()
    tmp_folder = f'{ctime}_tmpfolder'
    cmd = f'mkdir {tmp_folder}'
    os.system(cmd)
    tmp_gff_bed = f'{TMP_WD}/{tmp_folder}/{intergrationfile}_{ctime}_tmp_gff.bed'
    tmp_hout_bed = f'{TMP_WD}/{tmp_folder}/{intergrationfile}_{ctime}_tmp_hout.bed'
    tmp_pout_bed = f'{TMP_WD}/{tmp_folder}/{intergrationfile}_{ctime}_tmp_pout.bed'
    tmp_Intergenic = f'{TMP_WD}/{tmp_folder}/{intergrationfile}_{ctime}_tmp.Intergenic'
    tmp_Intergenic2 = f'{TMP_WD}/{tmp_folder}/{intergrationfile}_{ctime}_tmp.Intergenic2'
    tmp_Intergenic1 = f'{TMP_WD}/{tmp_folder}/{intergrationfile}_{ctime}_tmp.Intergenic1'
    tmp_Intergenic_ = f'{TMP_WD}/{tmp_folder}/{intergrationfile}_{ctime}_tmp.Intergenic_'
    Non_intergenic = f'{TMP_WD}/{tmp_folder}/{intergrationfile}_{ctime}_tmp.Non_intergenic'
    fo_phasiRNA = open(outfile, 'w')
    fo_allsiRNA = open(allfile, 'w')
    tmp_fo_intergration = f'{TMP_WD}/{tmp_folder}/{intergrationfile}_{ctime}_integration'
    fo_intergration_tmp = open(tmp_fo_intergration, 'w')
    fo_intergration = open(intergrationfile, 'w')
    # ----------> loading data <------------
    Gff3Tobed(gff3, tmp_gff_bed)
    TMP_WD = os.getcwd()
    ctime = GetCurTime()
    tmp_phasiRNAfile = f'{TMP_WD}/{tmp_folder}/{ctime}{intergrationfile}_phasiRNA.txt'
    tmp_allsiRNAfile = f'{TMP_WD}/{tmp_folder}/{ctime}{intergrationfile}_allsiRNA.txt'
    os.system(f'sort -k1,1 -k13,13 -k3,3 {phasiRNAfile} > {tmp_phasiRNAfile}')
    os.system(f'sort -k1,1 -k13,13 -k3,3 {allsiRNAfile} > {tmp_allsiRNAfile}')
    phasiRNA = LoadData(tmp_phasiRNAfile)
    tmp = LoadallsiRNA(tmp_allsiRNAfile)
    allsiRNA = tmp[0]
    p_allsiRNA = tmp[1]
    pp_allsiRNA = tmp[2]
    new_allsiRNA = tmp[1]
    hgdna_allsiRNA = allsiRNA['hg']
    hgdna_phasiRNA_cluster = phasiRNA[0]['hg']
    hgdna_phasiRNA_strand = phasiRNA[1]['hg']
    pout = nestedDic()
    pout['pg']['phasiRNA'] = phasiRNA[5]['pg']
    pout['pg']['allsiRNA'] = pp_allsiRNA['pg']
    pout['pc']['phasiRNA'] = phasiRNA[5]['pc']
    pout['pf']['phasiRNA'] = phasiRNA[5]['pf']
    pout['pc']['allsiRNA'] = pp_allsiRNA['pc']
    pout['pf']['allsiRNA'] = pp_allsiRNA['pf']
    pg_phasiRNA_cluster = GetPhaseScoreCluster(phasiRNA[3]['pg'])
    pc_phasiRNA_cluster = GetPhaseScoreCluster(phasiRNA[3]['pc'])
    pf_phasiRNA_cluster = GetPhaseScoreCluster(phasiRNA[3]['pf'])
    p_phasiRNA_cluster = (pc_phasiRNA_cluster, pg_phasiRNA_cluster, pf_phasiRNA_cluster)
    setting1 = (hgdna_phasiRNA_cluster, hgdna_phasiRNA_strand, phase_number, phase_length, pvalue_cutoff, min_read_num)
    hout = ParallelHypergeometric(allsiRNA, phasiRNA, setting1, parallel_number, island_number)

    def ConvertFormat_New(allsiRNAfile, hout):
        with open(allsiRNAfile, 'r') as fn:
            for line in fn:
                l = line.strip().split("\t")

    ConvertFormat_New(allsiRNAfile, hout)
    intron_infor = GeneratingIntronInforFromGff(tmp_gff_bed)
    WringIntronInforToGff(intron_infor, tmp_gff_bed)
    transcriptAnno = ParseTranscriptFeature(tmp_gff_bed)
    GetPhasiRNAClusterBed(hout, tmp_hout_bed, 'h')
    GetPhasiRNAClusterBed(pout, tmp_pout_bed, 'p')
    if flnc_anno != '':
        flnc_anno_dic = loadFLNCanno(flnc_anno)
    else:
        flnc_anno_dic = ''
    hout_anno_io = AnnoGdnaWithBedtools(tmp_hout_bed, tmp_gff_bed)
    pout_anno_io = AnnoGdnaWithBedtools(tmp_pout_bed, tmp_gff_bed)
    hanno = ParseBedAnno(hout_anno_io)
    panno = ParseBedAnno(pout_anno_io)
    intergration = Write_new(hout, pout, hanno, panno, p_phasiRNA_cluster, fo_phasiRNA, fo_allsiRNA, transcriptAnno, flnc_anno_dic)
    fo_phasiRNA.close()
    feature_list = ['lncRNA', 'Other', 'mRNA', 'rRNA', 'transcript', 'intron', 'FLNC', 'lnc_RNA']
    filter_intergration = FormatForWriteInterGration(intergration, feature_list)
    relation_dic = filter_intergration[1]
    tmp_out_list = WriteIntergration_new(filter_intergration[0])
    if gdna_enable == 'n':
        duplication_dic = FinalOut(tmp_out_list, relation_dic, fo_intergration)
    elif gdna_enable == 'y':
        duplication_dic = FinalOut(tmp_out_list, relation_dic, fo_intergration_tmp)
        fo_intergration_tmp.close()
        tmp_variable = overlapIntergenic(tmp_fo_intergration, tmp_Intergenic, tmp_Intergenic1, tmp_Intergenic_, tmp_Intergenic2)
        fo_Non_intergenic = open(Non_intergenic, 'w')
        with open(tmp_fo_intergration, 'r') as fn:
            for line in fn:
                line1 = line.strip()
                l = line.strip().split("\t")
                feature = l[1]
                if feature != 'Intergenic':
                    fo_Non_intergenic.write(line1 + '\n')
        fo_Non_intergenic.close()
        overlap_Intergenic_fo = tmp_variable
        Integration_overlap(overlap_Intergenic_fo, fo_intergration)
        fo_intergration.close()
        catCombine(Non_intergenic, intergrationfile)
    Vprint('intergration finished!', enable=True)
    fo_allsiRNA.close()
    os.system(f'rm -r {TMP_WD}/{tmp_folder}')

if __name__ == "__main__":
    main1()