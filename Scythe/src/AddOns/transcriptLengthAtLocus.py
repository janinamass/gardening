import getopt
import imp
import sys
import os.path
import re
gffFastaTools = imp.load_source('gffFastaTools', '../BioHelpers/gffFastaTools.py')
clusteringTools = imp.load_source('clusteringTools', '../BioHelpers/clusteringTools.py')

from gffFastaTools import *
VERBOSE = False
#####################################
# [last update 4/05/13 by JM]       #
# DEV                               #
#####################################
"""
Read loc and fasta file(s) 
and return table: loc\tlen(tr_fa1)\tlen(tr_fa2)\t...\tlen(tr_faN) 
"""

def usage():
    print ("""
    ##########################################################################
    # python transcriptLengthAthLocus.py 
    #     -l locFile -f listOfFastaFiles [-o outfile] 
    #########################################################################
    
    general options:

    -l, --portho=PROTEINORTHO.out    proteinortho output as input
    -f, --fasta=FA1,FA2,...,FAn      ','-separated list of fasta files 
                                        (of all species in PROTEINORTHO.out)
    optional:
    -n, --numid=CLUSTERNUM0,CLUSTERNUM1,...  only write files for orthogroup CLUSTERNUM
    
    -i, --id=FASTAHEADER0,FASTAHEADER1,...   only write files for orthogroup including FASTAHEADER
    -r, --regex=PATTERN               like -i but with one RegEx pattern
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
    portho= None
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
        opts, args = getopt.gnu_getopt(sys.argv[1:], "f:p:r:n:i:o:d:hv", ["fasta=","portho=","regex=","numid=","id=","outdir=","drop_prefix=","help","verbose"])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
    for o, a in opts:
        if o in ("-f", "--fasta"):
            fasta = a
        elif o in ("-p", "--portho"):
            portho = a
        elif o in ("-n", "--numid"):
            id = a
        elif o in ("-r", "--regex"):
            regex = a
        elif o in ("-i", "--id"):
            fastaid = a
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
    if portho is None:
        usage()
    if outdir: 
        if not os.path.exists(outdir):
            os.makedirs(outdir)
    ########################################
    # look for orthogroup with id 'id'
    if id:
        id = id.split(",")
        print(id)
        for g in clusteringTools.ProteinOrthoParser().readGroups(portho):
            if (str(g[0]) in id):
                tmp = [h for h in g[1] if h !="*" ]
                headers[g[0]]=tmp
                print(headers[g[0]])
    # fasta header
    if fastaid:
        fastaid = fastaid.split(",")
        for g in clusteringTools.ProteinOrthoParser().readGroups(portho):
            tmpheader = [f for f in fastaid if f in g[1]]
            if (tmpheader):
                tmp = [h for h in g[1] if h !="*" ]
                headers[g[0]]=tmp
    # regex
    if regex:
        p = re.compile(regex)
        print(regex)
        print(p)
        for g in clusteringTools.ProteinOrthoParser().readGroups(portho):
            tmpheader = [g for g in g[1] if p.match(g)]
            if (tmpheader):
                tmp = [h for h in g[1] if h !="*" ]
                headers[g[0]]=tmp
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
                            print("cleared", cleared)
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

