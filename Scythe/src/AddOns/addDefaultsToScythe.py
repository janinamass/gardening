import imp
import sys
import subprocess
import getopt
gffFastaTools = imp.load_source('gffFastaTools', '../BioHelpers/gffFastaTools.py')
from gffFastaTools import *
VERBOSE = False
#####################################
# [last update 4/1/13 by JM]       #
# DEV                               #
#####################################


def default_from_loc(loc):
    loc = open(loc, 'r')
    # if no default was specified during the .loc file creation
    # this will yield artefacts
    default_models = set()
    models2locus = {}
    for ln in loc:
        nondefault_models = None
        ln = ln.rstrip()
        #print(ln)
        try:
            tmp = ln.split("\t")
            locus = tmp[0]
            default_model = tmp[1]
            nondefault_models = tmp[2:]
        except IndexError as ie:
            try:
                tmp =ln.split("\t")
                locus = tmp[0]
                default_model =tmp[1]
            
                #locus,default_model = ln.split("\t")[0][1]
            except ValueError as ie:
                print("Error: ",ln," does not have any gene model. Please check your .loc file\n")
                print(ie)
                exit(2)
        default_models.add(default_model)
        models2locus[default_model] = locus
        for m in nondefault_models:
            models2locus[m]=locus
    return (default_models, models2locus)

def get_fasta_headers(fasta, sep):
    fasta_headers = set()
    for fa in gffFastaTools.FastaParser().read_fasta(fasta, delim=sep, asID=1):
        header= fa[0]
        fasta_headers.add(header)
    return (fasta_headers)
    
#deprecated, use default_from_loc instead
def read_gff2loc(infile):
    infile = open(infile, 'r')
    loci = {}
    longest = set()
    rawstr = r"""(Name=)(.*);pacid.*(longest=)(.*);(Parent=)(.*)"""
    for ln in infile:
        s =ln
        m = re.findall(rawstr, s)
        if len(m)  >0:
            name =  m[0][1]
            isLongest = m[0][3]
            parent = m[0][5]
            if isLongest == str(1):
                longest.add(name)
            if parent in loci:
                loci[parent].append(name)
            else:
                loci[parent]=[name]
    s = sorted(loci.keys())
    return loci,longest


#gff = sys.argv[1]
#fa1 = sys.argv[2]
#fa2 = sys.argv[3]
#findMax = True
#findMin = False
#
#notCovered = 0
#coveredButSame= 0
#total = 0
#fromFasta1 = set()
#fromFasta2 = set()
#fasta1 = {}
#fasta2 = {}
#test,longest = read_gff2loc(gff)
#
#for a in GFFFastaTools.FastaParser().read_fasta(fa1, delim="|", asID=0):
#    fromFasta1.add(a[0])
#    fasta1[a[0]] =a[1]
#for a in GFFFastaTools.FastaParser().read_fasta(fa2, delim="|", asID=1):
#    fromFasta2.add(a[0])
#    fasta2[a[0]] =a[1]
#
#for k in test:
#    total+=1
#    tmplist = [x for x in fromFasta2 if x in test[k]]
#    if not tmplist:
#        tmplong = [x for x in test[k] if x in longest]
#        notCovered+=1
#        try:
#            print(">"+tmplong[0]+"\n"+fasta1[tmplong[0]])
#        except IndexError as e:
#            #sys.stderr.write(e)
#            sys.stderr.write(" ".join(test[k]))
#    else:
#         print(">"+tmplist[0]+"\n"+fasta2[tmplist[0]])
#         if tmplist[0] in longest:
#             coveredButSame+=1
#
#print("##","\t".join(["total", "not_covered", "covered", "concurs"]))
#
#print("#",total, notCovered, total-notCovered,coveredButSame )



def usage():
    print ("""
    ##########################################################
    # python3  xxx.py -l file.loc  -s res.fasta -f full.fasta#
    ##########################################################
    
    general options:
    
    -l, --loc=FILE.LOC 
                (see Tools/inputFormatsReadme.txt)
    -s, --scythe_result=RESULT.FASTA    
                single output FASTA file from scythe run
    -f, --full_fasta=ALL_MODELS.FASTA
                the FASTA file with all gene models
    -v, --verbose 
    -h, --help  prints this
    
    advanced options:
    
    -p, --drop_prefix=SEPARATOR    
                 ignore everything before SEPARATOR (incl)
    -t, --table  print tabular summary
    """)
    sys.exit(2)
    
def annoy(*args):
    if VERBOSE:
        print("# ", *args)
    else:
        pass


def main():
    global VERBOSE
    loc = None
    scytheFasta = None
    sep = None
    fullFasta = None
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "l:s:f:p:thv", ["loc=","scythe_result=","full_fasta=","drop_prefix=","table","help","verbose"])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
    for o, a in opts:
        if o in ("-l", "--loc"):
            loc=a
        elif o in ("-s", "--scythe_result"):
            scytheFasta = a
        elif o in ("-f", "--full_fasta"):
             fullFasta = a
        elif o in ("-t", "--table"):
            table = True
        elif o in ("-h", "--help"):
            usage()
        elif o in ("-v", "--verbose"):
            VERBOSE = True    
        elif o in ("-p", "--drop_prefix"):
            sep = a
        else:
            assert False, "unhandled option"
    ########################################
    if loc is None:
        usage()
    if scytheFasta is None:
        usage()
    if fullFasta is None:
        usage()
    ########################################
    default_models, mod2loc = default_from_loc(loc)
    scythe_headers = get_fasta_headers(scytheFasta, sep)
    covered = set()
    for sh in scythe_headers:
        covered.add(mod2loc[sh])
        #print(mod2loc[sh]," covered\n")
    
    additional_default_models = [x for x in default_models if mod2loc[x] not in covered]
    #print("additionally, "+str(len(additional_default_models)))
    #print(len(covered)+len(additional_default_models))
    outfile = scytheFasta+".merged"
    outfile = open(outfile, "w")
    for ln in open(scytheFasta,"r"):
        #rm prefix
        if (ln.startswith(">")):
            ln =  ">"+ln.split("|")[1]
        outfile.write(ln)
    
    for tmp in gffFastaTools.FastaParser().read_fasta(fullFasta, delim=sep, asID=0):
        header, seq= tmp[0], tmp[1]
        if (header in additional_default_models):
            outfile.write(">"+header+"\n"+seq+"\n")
            additional_default_models.remove(header)
    #    fasta_headers.add(header)
    #return (fasta_headers)
    #print(default_models)
    outfile.close()
if __name__ == "__main__": 
    main()
