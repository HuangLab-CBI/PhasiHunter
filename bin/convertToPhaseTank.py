import sys
from customeDic import *
from Bio import SeqIO

def convertFa(inp):
    count = 1
    for query in SeqIO.parse(inp,"fasta"):
        name = query.id
        seq = str(query.seq)
        ele = name.split('@')
        abun = ele[1]
        print(f'>t{count}_x{abun}\n{seq}')
        count += 1

def convertgDNA(inp):
    for query in SeqIO.parse(inp, 'fasta'):
        name = query.id
        seq = query.seq
        ele = name.split(' ')
        print(f'>{ele[0]}\n{seq}')

def main():
    help = '''
    phase usage:
        option:
            # necessary options:
            -i: file  --  fasta file
            -m: int  --  function
                        [1]: convertFa, >SRR062265#5@97.21\\nGTCGTTGTAGTATAGTGGT to >t1_x97.21\\nGTCGTTGTAGTATAGTGGT
                        [2]: convertgDNA
                        [3]: integration phasiRNA

            # other
            -v:       --  print version information
            -h:       --  print help information
    '''
    version = '''
    '''
    for i in range(1, len(sys.argv)):
        if sys.argv[i] == '-i':
            inp = sys.argv[i+1]
        elif sys.argv[i] == '-v':
            print(version)
            sys.exit()
        elif sys.argv[i] == '-m':
            method = int(sys.argv[i+1])
        elif sys.argv[i] == '-h':
            print(help)
            sys.exit()

    rDic = {
        1: convertFa,
        2: convertgDNA,
        3: integrationPhasiRNA,
    }

    if method in rDic:
        rDic[method](inp)

def integrationPhasiRNA(inp):
    dic = OneDepDic()
    for query in SeqIO.parse(inp, 'fasta'):
        descript = query.description
        n_descript = descript.replace(' ', '_')
        seq = str(query.seq)
        dic[seq].append(n_descript)
    for i in dic:
        tmp = ";".join(dic[i])
        print(f'>{tmp}\n{seq}')

if __name__ == "__main__":
    main()