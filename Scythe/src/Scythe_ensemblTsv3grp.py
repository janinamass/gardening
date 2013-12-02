import sys, getopt
#####################################
# last update 23/10/2013 by J. Mass #
# version = '0.1'                   #
#####################################

def usage():
    print ("""
    ######################################
    #  Scythe_ensembleTsv2grp.py (v0.1)  #
    ######################################
    -f, --files=list of ensembl tsv files (eg sA.tsv,sB.tsv,sC.tsv)
    -o, --output=FILE         output file
    -h, --help                prints this
    #----------------------------------#
    
    """)
    sys.exit(2)
def main():
    ###################################
    outfile = None
    infiles = None
    ###################################
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "f:ho:l:", ["files=","help", "output=","lengths_file="])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
    for o, a in opts:
        if o in ("-f", "--files"):
            infiles = a
        elif o in ("-h", "--help"):
            usage()
        elif o in ("-o", "--output"):
            outfile = a
        else:
            assert False, "unhandled option"

    if not infiles:
        usage()
    if not outfile:
        outfile = "out.grp"
    infiles = infiles.split(",")
    readTsvFiles(infiles, outfile)
########################################
#def pfloat(float):
#    return "%0.2f" % float
#def printInfo(numLoc, numTr, outfile=outfile, verb=False):
#    if verb:
#        print("# Success\n# Formatted "+str(numLoc)+" loci and "+str(numTr)+" transcripts.")
#        print("# That's an average of "+pfloat(float(numTr)/float(numLoc))+" Transcripts per Locus")
#        print("# Saved to file "+outfile)

def readTsvFiles(listoftsv, outfile):
    print(outfile)
    if listoftsv is None:
        return(-1)
    ortho=dict()#
    seen = set()
    done = set()
    res = ""
    for t in listoftsv:
        #print(t)
        infile = open(t, "r")
        for l in infile:
            l = l.strip()
            l = l.split("\t")
           # print(l[0], " - ", l[1])
            seen.add(l[0])
            seen.add(l[1])
            
            if l[0] not in ortho:
                ortho[l[0]] = [l[1]]
                #print("new")
            else:
                ortho[l[0]].append(l[1])
               # print("add")
            if l[1] not in ortho:
                ortho[l[1]] = [l[0]]
               # print("new2")
            else:
                ortho[l[1]].append(l[0])
               # print("add2")
        
    
    #for o in ortho.keys():
    #    print(o)
        #print(o,"+",ortho[o])
        #if len(ortho[o])>1:
            #print(o,ortho[o])         
    #print(len(ortho.keys()))
    #print(len(seen))
    cntr = 0
    for s in seen:
        if s not in done:
            res+=str(cntr)+"\t"+s+"\t"+"\t".join(ortho[s])
            res+="\n"
            cntr+=1
            done.add(s)
            for d in ortho[s]:
                done.add(s)
        
    #print(res)
    out=open(outfile,"w")
    #print(res)
    out.write(res)
if __name__ == "__main__":
    main()