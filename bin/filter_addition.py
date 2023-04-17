import sys
'''
@File    :   filter_addition.py
@Time    :   2022/12/05 18:00:14
@Author  :   zrf 
@Version :   0.0
@Contact :   zrfeng1@gmail.com
@License :   (C)Copyright 2022-, NJAU-CBI
@Desc    :   None
'''
# ==============================================================> declare variable and function --> <================================================================ #
def FilterStringForList(list_: list, string: str):
    """return filterd list

    Arguments:
        list_ {list} -- list
        string {str} -- filtered string from list 

    Returns:
        list -- filtered list
    """
    return list(map(lambda x: x.replace(string, ''), list_))

def GenerateListFromFile(file: str):
    """generate a list from a file

    Arguments:
        file {str} -- _description_

    Returns:
        list -- _description_
    """
    tmpList = []
    with open(file, 'r') as fn:
        for line in fn:
            line1 = line.strip()
            tmpList.append(line1)
    return tmpList

help = '''
    this script is designed to combine libararies according to the sequence
    option:
        -i: <file>  --  A file containing the filename to be merged
        -t: [str]   --  verbose | compact; default is compact
        -o: <file>  --  combined outputfile
'''
version = '''
'''
type_ = 'compact'
# ==============================================================> declare variable and funtion <-- <================================================================ #
for i in range(1, len(sys.argv)):
    if sys.argv[i] == '-i':
        inp = sys.argv[i+1]
    elif sys.argv[i] == '-o':
        out = sys.argv[i+1]
    elif sys.argv[i] == '-t':
        type_ = sys.argv[i+1]
    elif sys.argv[i] == '-v':
        print(version)
        sys.exit()
    elif sys.argv[i] == '-h':
        print(help)
        sys.exit()

inpList = GenerateListFromFile(inp)
tmpDic = {}

for j in inpList:
    with open(j, 'r') as fn:
        lines = fn.readlines()
        lines = FilterStringForList(lines, "\n")
        for i in list(range(1, len(lines), 2)):
            srr = lines[i-1]
            srnaId = srr.split("@")[0]
            seq = lines[i]
            abun = float(srr.split("@")[-1])
            if seq not in tmpDic:
                tmpDic[seq] = [abun, srnaId]
            else:
                tmpDic[seq][0] += abun
                tmpDic[seq][1] = tmpDic[seq][1] + '_' + srnaId
        print(f'{j} finished')

if type_ == 'verbose':
    with open(out, 'w') as fo:
        for i in tmpDic:
            fo.write(tmpDic[i][1] + '_' + str(tmpDic[i][0]) + "\n")
            fo.writelines(i + "\n")
elif type_ == 'compact':
    id_prefix = '>t'
    count = 0
    with open(out, 'w') as fo:
        for i in tmpDic:
            id = id_prefix + str(count)
            fo.write(id + '_' + str(tmpDic[i][0]) + "\n")
            fo.writelines(i + "\n")
            count += 1