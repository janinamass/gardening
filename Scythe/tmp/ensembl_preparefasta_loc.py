from gffFastaTools import FastaParser


class Pep(object):
    def __init__(self,pep,gene,length ):
        self.gene = gene
        self.pep  = pep
        self.length = length
        self.isLongest = None
def prepareLocFromFasta(fasta):
    print("PL")
    genes = {}
    longest = {}
    out = fasta+".tmp"
    fp = FastaParser()
    #fp.read_fasta(fasta)
    open(out,'w')
    for i in fp.read_fasta(fasta):
        tmp = i[0].split(" ")
        geneid = [t for t in tmp if "gene:" in t]
        geneid = geneid[0].split("gene:")[-1]
        proteinid = tmp[0]
        protlen=len(i[1])
        #print("\t".join([geneid, proteinid,str(protlen)]))
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
        #for v in genes[g]:
        print(firstcol[0].gene, firstcol[0].pep, firstcol[0].length)
        for v in restcol:    
            print(v.gene,v.pep,v.length)
        print("#")
    
def q(s):
    return ('"'+s+'"')


prepareLocFromFasta("homo_sapiens.fa")
    