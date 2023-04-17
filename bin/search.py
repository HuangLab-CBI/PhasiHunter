import os


def fullWordSearch(query, pool):
    cmd = f'cat {query} | rush ' + '\'grep {} ' + pool + '\''
    os.system(cmd)

def main():
    import sys
    help = '''
        phasiHunter search <query> <pool>
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

    query = sys.argv[1]
    pool = sys.argv[2]

    fullWordSearch(query, pool)

if __name__ == '__main__':
    main()