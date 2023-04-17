import sys
import re

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


def main():
    help = '''
    option:
        -i: file  --  tab format chr relationship
        -j: file  --  file need to format
        -o: out  --  chr converted file
    '''
    version = '''
    '''
    for i in range(1, len(sys.argv)):
        if sys.argv[i] == '-i':
            inp = sys.argv[i+1]
        elif sys.argv[i] == '-j':
            jnp = sys.argv[i+1]
        elif sys.argv[i] == '-o':
            out = sys.argv[i+1]
        elif sys.argv[i] == '-v':
            print(version)
            sys.exit()
        elif sys.argv[i] == '-h':
            print(help)
            sys.exit()
        
    relationDic = generatingRelationDci(inp)
    fo = open(out, 'w')
    with open(jnp, 'r') as fn:
        for line in fn:
            l = line.strip().split('\t')
            if l[0] in relationDic:
                fo.write(relationDic[l[0]] + '\t' + "\t".join(l[1:]) + "\n")

    fo.close()

def generatingRelationDci(file):
    dic = {}
    with open(file, 'r') as fn:
        for line in fn:
            l = reSplit('\s+', line.strip())
            dic[l[0]] = l[1]
            dic[l[1]] = l[0]
    return dic

if __name__ == "__main__":
    main()