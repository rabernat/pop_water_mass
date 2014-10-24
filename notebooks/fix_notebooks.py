import sys
import re

def main():
    in_head = False
    in_bad = False
    fname = sys.argv[1]
    with open(fname, 'r') as f:
        fout = open(fname + '.FIXED', 'w')
        for l in f:
            if re.match('<<<<<<< HEAD', l):
                in_head = True
            elif re.match('=======', l):
                in_head = False
                in_bad = True
            elif re.match('>>>>>>>', l):
                in_bad = False
            else:
                if not in_bad:
                    fout.write(l)
        fout.close()
        
if __name__ == "__main__":
    main()