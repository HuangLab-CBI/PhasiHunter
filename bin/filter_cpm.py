import sys
from Bio import SeqIO

help = '''
phase usage:
    option:
        # necessary options:
        -i:  file  --  format output fasta
        -c:  float --  cpm cutoff

        # other
        -v:       --  print version information
        -h:       --  print help information
'''
version = '''
'''
for i in range(1, len(sys.argv)):
    if sys.argv[i] == '-i':
        inp = sys.argv[i+1]
    elif sys.argv[i] == '-c':
        cpm_cutoff = float(sys.argv[i+1])
    elif sys.argv[i] == '-v':
        print(version)
        sys.exit()
    elif sys.argv[i] == '-h':
        print(help)
        sys.exit()

def main():
    for query in SeqIO.parse(inp,"fasta"):
        name = query.id
        ele = name.split("@")
        cpm = float(ele[1])
        if cpm < cpm_cutoff:
            continue

        seq = str(query.seq)
        Id = f'{ele[0]}@{ele[1]}'
        print(f'>{Id}\n{seq}')



if __name__ == "__main__":
    main()