import getopt
import imp
import sys
import os.path
import re
gffFastaTools = imp.load_source('gffFastaTools', 'BioHelpers/gffFastaTools.py')
clusteringTools = imp.load_source('clusteringTools', 'BioHelpers/clusteringTools.py')

from gffFastaTools import *
VERBOSE = False
#####################################
# [last update 04/16/13 by JM]      #
# DEV                               #
#####################################
#todo: throw error if id not in fasta
"""
Read Proteinortho output and fasta file (all species, concatenated) 
and return one fasta file for each orthogroup.
Read additional .loc file to filter out single/multi model loci.
"""

def usage():
    print ("""
    ######################################################################################
    # python orthogroup2fasta_multi.py 
    #     -p proteinOrtho.out -f concatenatedfasta.fa [-i ID] [-o outdir] [-l loc locfile]
    ######################################################################################
    
    general options:

    -p, --portho=PROTEINORTHO.out    proteinortho output as input
    -f, --fasta=CONCAT.FA            concatenated fasta file 
                                        (of all species in PROTEINORTHO.out)
    optional:
    -n, --numid=CLUSTERNUM0,CLUSTERNUM1,...  only write files for orthogroup CLUSTERNUM
    
    -R, --range				     use range from numid_0 to numid_1 instead (eg -n0,3 -R prints 0,1,2)
    -i, --id=FASTAHEADER0,FASTAHEADER1,...   only write files for orthogroup including FASTAHEADER
    -r, --regex=PATTERN               like -i but with one RegEx pattern
    -o, --outdir=OUTPUTDIR            write files to OUTPUTDIR instead of ./
    
    -l, --loc=LOCFILE                 read .loc file (concatenated) and split output 
                                      into groups with only single model loci, 
                                      only multi model loci and the rest
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
    locfile = None
    multiModels = set()
    allMulti = dict()
    mixedMulti = dict()
    nonMulti = dict()
    is_range = False
    ###################################
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "f:p:r:n:i:o:l:d:hvR", ["fasta=","portho=","regex=","numid=","id=","outdir=","loc=","drop_prefix=","help","verbose", "range"])
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
        elif o in ("-l", "--loc"):
            locfile = a
        elif o in ("-o", "--outdir"):
            outdir = a
        elif o in ("-h", "--help"):
            usage()
        elif o in ("-v", "--verbose"):
            VERBOSE = True    
        elif o in ("-d", "--drop_prefix"):
            sep = a
        elif o in ("-R", "--range"):
            is_range=True
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
    if (locfile):
        locfile = open(locfile,'r')
        for ln in locfile:
            ln = ln.rstrip()
            loc, gm = ln.split("\t")[0], ln.split("\t")[1:]
            if len(gm)>1:
                for g in gm:
                    multiModels.add(g)
        #print(multiModels) 
    # look for orthogroup with id 'id'
    if id:
        id = id.split(",")
        #print(id)
        if is_range: 
            id = [str(a) for a in list(range(int(id[0]),int(id[1])))]
            #print(id)
        for g in clusteringTools.ProteinOrthoParser().readGroups(portho):
            if (str(g[0]) in id):
                #print(g[0],id)
                tmp = [h for h in g[1] if h !="*" ]
                tmpAllMulti = [t for t in tmp if t in multiModels]
                print("allMuti ",tmpAllMulti)
                if len(tmpAllMulti) == len(tmp): 
                    print("all multi")
                    allMulti[g[0]]=tmp
                elif len(tmpAllMulti) < 1:
                    print("non multi")
                    nonMulti[g[0]]=tmp
                else:
                    print("mixedMulti")
                    mixedMulti[g[0]]=tmp
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
                filek_am = filek+".all_multi"
                filek_mm = filek+".mixed_multi"
                filek_nm = filek+".all_single"
                if tmp in headers[k]:
                        if not k in cleared:
                            #cleared.add(k)
                            #print("cleared", cleared)
                            if k in allMulti:
                                out2n = filek_am
                                #out2 = open(filek_am,'w')
                            elif k in mixedMulti:
                                out2n = filek_mm
                            elif k in nonMulti:
                                print(k, "all single model -> skip")
                                out2n = filek_nm
                            else:
                                out2n = ".tmpfile"
                            #out2 = open(out2n,'w')
                      #      out = open(filek,'w')
                      #      out.close()
                            #out2.close()
                        out2 = open(out2n,'a+')
                      #  out = open(filek,'a+')
                      #  print("out", tmp,filek)
                      #  out.write(">"+tmp+"\n"+a[1]+"\n")
                      #  out.close()
                        print("out: ",tmp,out2n)
                        out2.write(">"+tmp+"\n"+a[1]+"\n")
                        out2.close()
            else:
                if tmp in headers[k]:
                        print(tmp)
    #annoy("done with "+fasta)
    
 
 
if __name__ == "__main__": 
    main()

