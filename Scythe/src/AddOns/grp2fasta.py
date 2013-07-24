import getopt
import imp
import sys
import os.path
import re
gffFastaTools = imp.load_source('gffFastaTools', '../BioHelpers/gffFastaTools.py')
#clusteringTools = imp.load_source('clusteringTools', '../BioHelpers/clusteringTools.py')

from gffFastaTools import *
VERBOSE = False
#####################################
# [last update 04/11/13 by JM]      #
# DEV                               #
#####################################
#todo: throw error if id not in fasta
"""
Read .grp, loc and fasta file (all species, concatenated) 
and return one fasta file for each orthogroup with the default sequence.
"""

def usage():
    print ("""
    ##########################################################################
    # python grp2fasta.py 
    #     -g file.grp -l concatenated.loc -f concatenatedfasta.fa [-i ID] [-o outdir] 
    #########################################################################
    
    general options:

    -g, --grp=file.grp                  grp as input
    -f, --fasta=CONCAT.fa               concatenated fasta file 
                                        (of all species in PROTEINORTHO.out)
    -l, --loc=CONCAT.loc
    
    optional:
    -n, --numid=CLUSTERNUM0,CLUSTERNUM1,...  only write files for orthogroup CLUSTERNUM
    -o, --outdir=OUTPUTDIR            write files to OUTPUTDIR instead of ./
    
    -v, --verbose 
    -h, --help  prints this
    
    advanced options:
    
    -d, --drop_prefix=SEPARATOR    
                 ignore everything before SEPARATOR (incl)
    """)
    sys.exit(2)
    
def annoy(*args):
    if VERBOSE:
        print("# ", *args)
    else:
        pass

def main():
    global VERBOSE
    fasta = None
    grp = None
    regex = None
    id = None
    fastaid = None
    ofInterest = None
    headers = dict()
    outdir = "./"
    sep = None
    cleared = set()
    ###################################
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "f:g:l:n:o:d:hv", ["fasta=","grp=","loc=", "numid=","outdir=","drop_prefix=","help","verbose"])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
    for o, a in opts:
        if o in ("-f", "--fasta"):
            fasta = a
        elif o in ("-g", "--grp"):
            grp = a
        elif o in ("-l", "--loc"):
            loc = a
        elif o in ("-n", "--numid"):
            id = a
        elif o in ("-o", "--outdir"):
            outdir = a
        elif o in ("-h", "--help"):
            usage()
        elif o in ("-v", "--verbose"):
            VERBOSE = True    
        elif o in ("-d", "--drop_prefix"):
            sep = a
        else:
            assert False, "unhandled option"
    
    ########################################
    if fasta is None:
        usage()
    if grp is None:
        usage()
    if outdir: 
        if not os.path.exists(outdir):
            os.makedirs(outdir)
    ########################################
    grpfile = open(grp,"r")
    locfile = open(loc,"r")
    locdct = {}
    headers = {}
    for l in locfile:
        l=l.rstrip()
        locdct[l.split("\t")[0]] = l.split("\t")[1]
    
    grpdct = {}
    for l in grpfile:
        l=l.rstrip()
        try:
            tmp = [locdct[i] for i in l.split("\t")[1:]]
        except KeyError as ke:
            print(ke)
        grpdct[l.split("\t")[0]] = tmp

    headers = grpdct
    for a in gffFastaTools.FastaParser().read_fasta(fasta):
        tmp = a[0]
        annoy(tmp, fasta)
        if sep:
            try:
                tmp = a[0].split(sep)[0]
            except IndexError as ie:
                annoy("#", a[0])
            annoy("->",tmp, fasta)

        for k in headers:
            if outdir:
                filek = outdir+"/"+str(k)+".ortho.fasta"
            if tmp in headers[k]:
                    if not k in cleared:
                        cleared.add(k)
                        #print("cleared", cleared)
                        out = open(filek,'w')
                        out.close()
                    out = open(filek,'a+')
                    out.write(">"+tmp+"\n"+a[1]+"\n")
                    out.close()
            else:
                if tmp in headers[k]:
                        print(tmp)
    #annoy("done with "+fasta)
    
 
 
if __name__ == "__main__": 
    main()

