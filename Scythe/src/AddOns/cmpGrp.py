#!/usr/bin/python
import re,sys, getopt

###################
# 04/15/2013 [JM] #
# v0.1            #
###################

def usage():
    print ("""
    ################################
    #  Scythe_cmpGrp.py  #
    ################################
   
    general options:
    -a, --grp_a=FILE_A.grp     
    -b, --grp_b=FILE_B.grp
    -o, --output=OUTFILE     output file [default: FILE_A.FILE_B.cmp]
    -h, --help               prints this    
    
    ------------
    .grp format: GroupID\tgeneIDiSp1\tgeneIDjSp2\t...geneIDkSpn
    """)
    sys.exit(2)


def cpmGrp(grpA, grpB):
    res=""
    dctA = {}
    dctB = {}
    dctAt = {}
    dctBt = {}
    grpA = open(grpA, 'r')
    cntA = 0
    for ln in grpA:
        cntA+=1
        ln = ln.rstrip()
        kA = ln.split("\t")[0]
        dctA[kA] = set(ln.split("\t")[1:]) 
    grpA.close()
    
    cntB = 0
    grpB = open(grpB, 'r')
    for ln in grpB:
        cntB+=1
        kB = ln.split("\t")[0]
        ln = ln.rstrip()
        dctB[kB] = set(ln.split("\t")[1:]) 
    grpB.close()
    
    dctAt = dctA.copy()
    dctBt = dctB.copy()
    
    idA = 0
    idB = 0
    subA = 0
    subB = 0
    superA = 0
    superB = 0
    for ka, va in dctAt.items():
       for kb, vb in dctBt.items():
           if va == vb:
               idA += 1
               idB += 1
               dctA.pop(ka)
               dctB.pop(kb)
           elif (va.issubset(vb)):
                subA +=1
                superB +=1
                
                try:
                    dctA.pop(ka)
                    dctB.pop(kb)
                    print("# ",ka, " subset of ",kb)
                except KeyError as k:
                    print("# ",ka, "also subset of ",kb)
           elif (va.issuperset(vb)):
                subB +=1
                superA +=1
                try:
                    dctA.pop(ka)
                    dctB.pop(kb)
                    print("# ",ka, " superset of ",kb)
                except KeyError as k:
                    print("# ",ka, "also superset of ",kb)
    #unrelated sets left in dctA and dctB
    unrelatedA = len(dctA.keys())
    unrelatedB = len(dctB.keys())
    #print(unrelatedA,  superA,  subA, idA, cntA)
    #print(unrelatedB,  superB,  subB, idB, cntB)
    
    #assert(unrelatedA + superA + subA + idA == cntA)
    #assert(unrelatedB + superB + subB +idB == cntB)
    res += "groups"+"\t"+str(cntA)+"\t"+str(cntB)+"\n"
    res += "identical"+"\t"+str(idA)+"\t"+str(idB)+"\n"
    res += "% identical (relative)"+"\t"+"{0:.2%}".format(idA/cntA)+"\t"+"{0:.2%}".format(idB/cntB)+"\n"
    res += "% identical (total)"+"\t"+"{0:.2%}".format((idA+idB)/(cntA+cntB))+"\t"+"{0:.2%}".format((idA+idB)/(cntA+cntB))+"\n"
    res += "subset"+"\t"+str(subA)+"\t"+str(subB)+"\n"
    res += "superset"+"\t"+str(superA)+"\t"+str(superB)+"\n"
    res += "other\t"+str(unrelatedA)+"\t"+str(unrelatedB)
    
    return(res)
def main():    
    ###################################
    grpA = None
    grpB = None
    out = None
    ###################################
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "a:b:o:h", ["grp_a=","grp_b=","output=", "help"])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
    for o, a in opts:
        if o in ("-a", "--grp_a"):
            grpA= a
        elif o in ("-b", "--grp_b"):
            grpB = a
        elif o in ("-h", "--help"):
            usage()
        elif o in ("-o", "--output"):
            out = a            
        else:
            assert False, "unhandled option"
    ########################################
    if grpA is None:
        usage()
    if grpB is None:
        usage()
    if out is None:
        out = grpA+"."+grpB+".cmp"
    ########################################
    print("#  Output: ", out)
    res = "files\t"+grpA+"\t"+grpB+"\n"
    res +=cpmGrp(grpA, grpB)
    print(res)
    out = open(out, "w")
    out.write(res)

if __name__ == "__main__":
    main()