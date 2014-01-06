#!/usr/bin/python3
import re,sys, getopt

###################
# 01/03/2014 [JM] #
###################

def usage():
    print ("""
    #####################################################
    #  filterOrthoAll.py -g groups.grp -n 3 -o new.grp  #
    #####################################################
   
    options:
    -g, --grp=FILE.grp     
    -o, --output=OUTFILE.grp output file [default: FILE.allspec.grp]
    [-r, --rename discard old orthogroup ids and start numbering from 0]    
    [-n, --numspec=N    min number of species ] 
    -h, --help         prints this    
    
    ------------
    .grp format: GroupID\tgeneIDiSp1\tgeneIDjSp2\t...geneIDkSpn
    """)
    sys.exit(2)



def filterRedundancy(grp, out):
    seen = set()
    ok = 0
    redundant = 0
    isOK = True
    try:
        g = open(grp,'r')
    except IOError as e:
        print(e)
    for l in g:
            l = l.strip()
            isOK = True
            tmp = l.split("\t")
            for t in tmp:
                #print(t)
                if t in seen:
                    isOK = False
                
                seen.add(t)
            if isOK:
                ok+=1
            else:
                    #print(tmp)
                redundant+=1
    print(ok)
    print(redundant)
def filterGroups(grp, out, numspec=None,rename=False):
    try:
        g = open(grp,'r')
    except IOError as e:
        print(e)
        usage()
    
    if not numspec:
        max = 0
        for l in g:
            l = l.strip()
            tmp = len(l.split("\t"))
            print(tmp)
            if tmp > max:
                max = tmp
        g.close()
    numspec=max
    try:
        g = open(grp,'r')
    except IOError as e:
        print(e)
        usage()
        
    try:
        print(out)
        out = open(out,'w')
        
    except IOError as e:
        print(e)
        usage()
    
    cnt=0
    for l in g:
        ln = l.strip()
        tmp = len(ln.split("\t"))
        if tmp == numspec:
            if not rename:
                out.write(l)
            else:
                out.write(str(cnt)+"\t"+"\t".join(l.split("\t")[1:]))
                cnt+=1
    out.close()

def main():
    out = None
    numspec = None
    rename=False
    ###################################
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "g:n:o:hr", ["grp=","numspec=","output=", "help","rename"])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
    for o, a in opts:
        if o in ("-g", "--grp"):
            grp= a
        elif o in ("-n", "--numspec"):
            numspec = int(a)
        elif o in ("-h", "--help"):
            usage()
        elif o in ("-o", "--output"):
            out = a     
        elif o in ("-r", "--rename"):
            rename = True
        else:
            assert False, "unhandled option"
    ########################################
    if grp is None:
        usage()
    if out is None:
        if grp.endswith(".grp"):
            tmp = grp[:-3]
            out = tmp+"allspec.grp"
    ########################################
    print("# Output: ", out)
   
    #out = open(out, "w")
    #out.write(res)
    filterGroups(grp, out, numspec, rename)
    filterRedundancy(out,"tmp")
    filterRedundancy(grp,"tmp")
if __name__=="__main__":
    main()