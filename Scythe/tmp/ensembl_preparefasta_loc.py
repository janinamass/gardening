from gffFastaTools import FastaParser
import sys, getopt


class Pep(object):
    def __init__(self,pep,gene,length ):
        self.gene = gene
        self.pep  = pep
        self.length = length
        self.isLongest = None
def prepareLocFromFasta(fasta, outfile):
    genes = {}
    longest = {}
    if not outfile:
        out = fasta+".loc"
    else:
        out = outfile
    fp = FastaParser()
    out = open(out,'w')
    for i in fp.read_fasta(fasta):
        tmp = i[0].split(" ")
        geneid = [t for t in tmp if "gene:" in t]
        geneid = geneid[0].split("gene:")[-1]
        proteinid = tmp[0]
        protlen=len(i[1])
        peptide= Pep(proteinid,geneid,protlen)
        try:
            genes[geneid].append(peptide)
            if peptide.length > longest[geneid].length:
                longest[geneid].isLongest =False
                longest[geneid]=peptide
                peptide.isLongest = True
        except KeyError as e:
            genes[geneid] = [peptide]
            longest[geneid] = peptide
            peptide.isLongest=True
    for g in genes:
        firstcol = [w for w in genes[g] if w.isLongest==True ]
        restcol = [w for w in genes[g] if w.isLongest!=True ]
        s = firstcol[0].gene+"\t"+firstcol[0].pep+"\t"
        s+="\t".join([v.pep for v in restcol])
        s+="\n"
        out.write(s)
def q(s):
    return ('"'+s+'"')
def usage():
    print ("""
    ##############################
    #
    ##############################
    
    -i, --infile=FASTA
    -o, --output=FILE        
    -h, --help
    """)
    sys.exit(2)
def main():

    ###################################
    infile = None
    outfile = None
    ###################################
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "i:o:h", ["infile=","outfile=","help"])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
    for o, a in opts:
        if o in ("-i", "--infile"):
            infile = a
        elif o in ("-h", "--help"):
            usage()
        elif o in ("-o", "--output"):
            outfile = a
        else:
            assert False, "unhandled option"

    if not infile:
        usage()
    
    
    prepareLocFromFasta(infile, outfile)

if __name__ == "__main__":
    main()