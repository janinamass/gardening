import sys, getopt
VERB=True
#####################################
# last update 03/05/2013 by J. Mass #
# version = '0.1'                   #
#####################################

def usage():
    print ("""
    ###################################
    #  Scythe_ensemble2loc.py (v0.1)  #
    ###################################
    -f, --file=ENSEMBL BioMart input
                              format: 1st column: gene id, 2nd column:transcript id;
                                      gene ids can occur multiple times
    -l, --lengths_file        format: gene id\ttranscript id\ttranscript length
    -o, --output=FILE         output file [default: ENSEMBLE input.loc]
    -h, --help                prints this
    #----------------------------------#
    .loc format:SP0G0\tSP0G0T0\tSP0G0T1\t...\tSP0G0TN
    """)
    sys.exit(2)
###################################
outfile = None
ensemblLoc = None
lengthsfile = None
###################################
try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "f:ho:l:", ["file=","help", "output=","lengths_file="])
except getopt.GetoptError as err:
    print (str(err))
    usage()
for o, a in opts:
    if o in ("-f", "--file"):
        ensemblLoc = a
    elif o in ("-h", "--help"):
        usage()
    elif o in ("-o", "--output"):
        outfile = a
    elif o in ("-l", "--lengths_file"):
        lengthsfile = a
    else:
        assert False, "unhandled option"
def main():
    if not ensemblLoc:
        usage()
    if not outfile:
        outfile = ensemblLoc+".loc"
########################################
def pfloat(float):
    return "%0.2f" % float
def printInfo(numLoc, numTr, outfile=outfile, verb=False):
    if verb:
        print("# Success\n# Formatted "+str(numLoc)+" loci and "+str(numTr)+" transcripts.")
        print("# That's an average of "+pfloat(float(numTr)/float(numLoc))+" Transcripts per Locus")
        print("# Saved to file "+outfile)
def readLength(lengthfile):
    lenDct={}
    lfile = open(lengthfile,'r')
    for l in lfile:
        if "Ensembl" in l:
            continue
        l=l.rstrip()
        try:
            geneid,trid,length = l.split("\t")
            #print(l)
            lenDct[trid]=int(length)
        except ValueError as ve:
           # print("# No length given ",ve,l)
            continue
        except IndexError:
            trid,length = l.strip().split("\t")
            print("Warning,less than 3 columns is that right:transcript: ",trid," length:"+length+" ?\n")
            lenDct[trid]=int(length)
    return (lenDct)


def readEnsemblLoc(ensembleLoc, outfile, lenDct=None):
    out=open(outfile,"w")
    infile = open(ensembleLoc,"r")
    numLoc = 0
    numTr = 0
    locDct = {}
    maxDct = {} #max length of transcript at locus
    longestTr = {}
    missing = False
    for l in infile:
        if "Ensembl" in l:
            continue
        l = l.strip()
        tmp = l.split("\t")
        if tmp[0] not in locDct:
            numLoc+=1
            numTr+=1
            locDct[tmp[0]]=[tmp[1]]
        else:
            numTr+=1
            locDct[tmp[0]].append(tmp[1])
    
    if lenDct:
        for l,tr in locDct.items():
            try:
                maxDct[l] = max([lenDct[i] for i in tr])
            except KeyError as ke:
                tmp=[]
                missing = True
                for i in tr:
                    if i in lenDct:
                        tmp.append(i)
                    else: 
                        lenDct[i]=0
                        tmp.append(i)
                maxDct[l] = max([lenDct[i] for i in tmp])
                if not maxDct[l]:
                    #last resort
                    maxDct[l]=tr[0]
            longestTr[l]=[t for t in tr if lenDct[t] == maxDct[l]]
            #print(longestTr[l])
            longestTr[l] = longestTr[l][0]
        locDct2 = locDct.copy()
        for loc,tr in locDct2.items():#rm
            locDct[loc] = [t for t in tr if t != longestTr[loc]]
        for loc,tr in locDct.items():
            out.write(loc+"\t"+longestTr[loc]+"\t"+"\t".join(tr)+"\n")
        if missing:
            print("# There were missing lengths\n")

    else: #no lengths given 
        for loc,tr in locDct.items():
            out.write(loc+"\t"+"\t".join(tr)+"\n")
    out.close()
    printInfo(numLoc=numLoc, numTr=numTr, outfile=outfile, verb=VERB)

    if lengthsfile:
       ldct = readLength(lengthsfile)
       readEnsemblLoc(ensemblLoc,outfile,ldct)
    else:
        readEnsemblLoc(ensemblLoc,outfile)
if __name__ == "__main__":
    main()