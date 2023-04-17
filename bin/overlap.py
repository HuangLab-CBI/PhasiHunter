# %%
import sys
import os


def main():
    help = '''
    usage:
        option:
            # necessary options:
            -a:  file  --  bed format
            -b:  file  --  bed format

            # options with default value

            # optional options

            # other
            -v:       --  print version information
            -h:       --  print help information
    '''
    version = '''
    '''
    for i in range(1, len(sys.argv)):
        if sys.argv[i] == '-a':
            afile = sys.argv[i+1]
        elif sys.argv[i] == '-b':
            bfile = sys.argv[i+1]
        elif sys.argv[i] == '-v':
            print(version)
            sys.exit()
        elif sys.argv[i] == '-h':
            print(help)
            sys.exit()


    # specific a
    a_specific = annotate(afile, bfile)
    b_specific = annotate(bfile, afile)
    overlap = intersect(afile, bfile)
    print(f'a specific: {a_specific}\na b overlap: {overlap}\nb specific: {b_specific}')


def annotate(afile, bfile):
    cmd = f"bedtools annotate -i {afile} -files {bfile}"  + " | awk '$4==0{print $0}' | wc -l"
    a_specific_io = os.popen(cmd)
    for i in a_specific_io:
        a_specific = i.strip()
    return a_specific

def intersect(afile, bfile):
    cmd = f"bedtools intersect -a {afile} -b {bfile} | wc -l"
    overlap_io = os.popen(cmd)
    for i in overlap_io:
        overlap = i.strip()
    return overlap

if __name__ == "__main__":
    main()