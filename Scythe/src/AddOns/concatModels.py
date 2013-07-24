#!/usr/bin/python
import re,sys, getopt

#####################################
# last update 10/26/2012 by J. Mass #
# version = '0.2'                   #
#####################################

def usage():
    print ("""
    ###################################
    #  concatModels.py v0.1    #
    ###################################
    
    -g, --gff=gff3_FILE      (tested w/ .gff3 from phytozome)
    -f, --fasta=fasta_FILE
    -o, --output=FILE        output file [default: fasta_FILE.concat]
    -n, --name=STRING        the concatenated ID gets this as suffix
    -d, --delim=STRING       split fasta headers at 'delim'
    -i, --asID=INT           use 'asID' part of fasta header as ID         
    -h, --help               prints this
    """)
    sys.exit(2)

def read_fasta(fasta, delim = None, asID = 0):
        """read from fasta fasta file 'fasta' 
        and split sequence id at 'delim' (if set)\n
        example:\n
        >idpart1|idpart2\n
        ATGTGA\n
        and 'delim="|"', 'asID'=0 returns ("idpart1", "ATGTGA")
        """
        name = ""
        fasta = open(fasta, "r")
        while True:
            line = name or fasta.readline()
            if not line:
                break
            seq = []
            while True:
                name = fasta.readline()
                name = name.rstrip()
                if not name or name.startswith(">"):
                    break
                else:
                    seq.append(name)
            joinedSeq = "".join(seq)
            line = line[1:]
            if delim:
                line = line.split(delim)[asID]
            yield (line.rstrip(), joinedSeq.rstrip())    
        fasta.close()

def concat(fasta, gff, outfile, suffix, delim, asID):
    
    gff = open(gff, 'r')
    outfile = open(outfile, 'w')
    parents = {}
    cds = {}
    seqs = {}
    rawstr = r"""(Name=)(.*);pacid.*(longest=)(.*);(Parent=)(.*)"""
    for ln in gff:
        s = ln
        m = re.findall(rawstr, s)
        if len(m) >0:
            name =  m[0][1]
            parent = m[0][5]
            parents[name] = parent
            if parent in cds:
                cds[parent].append(name)
            else:
                cds[parent]=[name]
    gff.close()
    
    for f in read_fasta(fasta, delim = delim, asID=asID):
        seqs[f[0]]=f[1]
    new = {}
    for loc in cds:
        new[loc]=""
        for c in cds[loc]:
            if c in seqs:
                new[loc]+=seqs[c]
       
    keys = new.keys()
    for k in sorted(new.keys()):
        str = ""
        if suffix:
           str =">"+k+suffix+"\n"+new[k]+"\n"
           print(str)
        else:
            str = ">"+" ".join(cds[k])+"\n"+new[k]+"\n"
            print(str)

    outfile.close()
    return None

###################################
outfile = None
fasta = None
gff = None
suffix = None
delim = None
asID = None
###################################
try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "f:g:o:s:d:i:h", ["fasta=","gff=", "output=", "suffix=", "delim=", "asID=","help"])
except getopt.GetoptError as err:
    print (str(err))
    usage()
for o, a in opts:
    if o in ("-f", "--fasta"):
        fasta = a
    elif o in ("-g", "--gff"):
        gff = a
    elif o in ("-o", "--output"):
        outfile = a
    elif o in ("-s", "--suffix"):
        suffix = a
    elif o in ("-d", "--delim"):
        delim = a
    elif o in ("-i", "--asID"):
        asID = int(a)
    elif o in ("-h", "--help"):
        usage()
    else:
        assert False, "unhandled option"

########################################
if fasta is None or gff is None:
    usage()
if outfile is None:
    outfile = fasta+".cc"
########################################
concat(fasta, gff, outfile, suffix, delim, asID)
