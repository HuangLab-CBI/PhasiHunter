import sys
#!/usr/bin/env python
# this program aims to identify the putative microRNA/sRNA initiator triggering the biogenesis of phased siRNAs
# By Ji Huang, Crop Bioinformatics Group, NJAU, 2018. (C)Copyright reserved.
# Only for internal test within CBG. No free distribution without the writing permission from CBG.
def main():
    help = '''
    phase usage:
        option:
            # necessary options:
            -i:  file  --  integration -o output
            -j:  file  --  the target predicted by psRNAtarget server or TarHunter.pl
            -o:  out  --  output file name

            # options with default value
            -pl: int   --  phase length, 21 | 24, default=21
            -pd: int   --  phasiRNA cluster extended island, default=5, extended_length=island*phase_length
            -ps: int   --  the position of cleavage at 10-11(0) or 9-12(1), default=1

            # optional options

            # other
            -v:       --  print version information
            -h:       --  print help information

    '''
    version = '''
    version v1.0
    '''

    for i in range(1, len(sys.argv)):
        if sys.argv[i] == '-i':
            inputfile = sys.argv[i+1]
        elif sys.argv[i] == '-j':
            psRNAtarget = sys.argv[i+1]
        elif sys.argv[i] == '-pd':
            distance = int(sys.argv[i+1])
        elif sys.argv[i] == '-pl':
            phaselength = int(sys.argv[i+1])
        elif sys.argv[i] == '-ps':
            shift = sys.argv[i+1]
        elif sys.argv[i] == '-o':
            outputfilename = sys.argv[i+1]
        elif sys.argv[i] == '-v':
            print(version)
            sys.exit()
        elif sys.argv[i] == '-h':
            print(help)
            sys.exit()

    if distance==0:
        if phaselength == 21:
            distance=105
        else:
            distance=120
    
    read_seq={}
    result=['miRNA\tGeneID\tCleavage_site\tTarget_pos_end\tPosmax/min\tDistance\tInitiator_seq\tInitiator_len\tRecordmarker\n'] 
    print("Loading psRNAtarget analysis data..")
    m_all={}
    for line in open (psRNAtarget, 'r'):
        line1=line.split('\t')
        if line.startswith('#') or line1[0] == 'miRNA_Acc.':
            continue
        else:
            Target_gene = line1[1]
            if Target_gene in m_all.keys():
                m_all[Target_gene]+=[line]
            else:
                m_all[Target_gene]=[line]

    record_all={}

    print("Loading phased siRNA file...")
    for line in open(inputfile, 'r'):
        if line.startswith('>') or line.startswith('#'):
            continue
        record=str(line.split('\t')[-1].strip())
        if record in record_all.keys():
            record_all[record]+=[line]
        else:
            record_all[record]=[line] #add record according to the record_markers
    print("Phase siRNA loaded.")
    print("Analysing...")
    minf={}
    maxf={}

    for key in record_all.keys(): #loop with record_markers
        # each record cluster, seperately
        minf['-']=99999999990 # what does minf means
        minf['+']=99999999999
        maxf['-']=0
        maxf['+']=0
        for i in range(len(record_all[key])):
            line=record_all[key][i].split('\t')
            geneid=line[0]
            strand=line[1]
            pos=int(line[2])
            if strand=='-': # for get the minimum pos at - strand
                if minf['-']>=pos:
                    minf['-']=pos
            elif strand=='+': # for get the minimum pos at - strand
                if minf['+']>=pos:
                    minf['+']=pos
            recordmarker=line[-1].strip()
            #print recordmarker
        if minf['-']<=minf['+']: 
            posmin=minf['-']+2
        else:
            posmin=minf['+']
        if geneid in m_all.keys(): 
            for j in range(len(m_all[geneid])): 
                line2=m_all[geneid][j].split('\t')
                miRNAid=line2[0]
                miRNAseq=line2[8].replace('U','T')
                miRNA_pos=int(line2[4])
                target_pos_start=int(line2[6])
                target_pos_end=int(line2[7])-1
                cleavage_site=target_pos_end-9
                if shift==0: 
                    #cleavage_site = target_pos_end - 9+ miRNA_pos-1
                    if posmin-cleavage_site>=0 and posmin-cleavage_site<=distance: 
                        if (posmin-cleavage_site)%int(phaselength)==0: 
                            result.append(miRNAid+'\t'+geneid+'\t'+str(cleavage_site)+'\t'+str(target_pos_end)+'\t'+str(posmin)+'\t'+str(posmin-cleavage_site)+'\t'+miRNAseq+'\t'+str(len(miRNAseq))+'\t'+str(recordmarker)+'\n')
                else:
                    cleavage_sites=[cleavage_site-1, cleavage_site, cleavage_site+1] 
                    for k in cleavage_sites:
                        if posmin-k>=0 and posmin-k<=distance:
                            
                            if (posmin-k)%int(phaselength)==0:
                                result.append(miRNAid+'\t'+geneid+'\t'+str(cleavage_site)+'\t'+str(target_pos_end)+'\t'+str(posmin)+'\t'+str(posmin-cleavage_site)+'\t'+miRNAseq+'\t'+str(len(miRNAseq))+'\t'+str(recordmarker)+'\n')
    resultw=open(outputfilename,'w')
    resultw.writelines(result)
    print("Analysis completed! %s records saved!"%(len(result)-1))



if __name__ == '__main__':
    main()