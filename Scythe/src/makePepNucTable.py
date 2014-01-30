import getopt
import imp
import sys
import os.path
import re
import glob
global BIOHELPERS
BIOHELPERS = "/home/janina/BioHelpers/"
gffFastaTools = imp.load_source('gffFastaTools', BIOHELPERS+'gffFastaTools.py')

from gffFastaTools import *
"""
Read fasta from ensembl with header lines like:
>ENSP00000451042 pep:known chromosome:GRCh37:14:22907539:22907546:1 
gene:ENSG00000223997 transcript:ENST00000415118 gene_biotype:TR_D_gene 
transcript_biotype:TR_D_gene
"""

def usage():
    print ("""
    #########################################
    # python3 
    #########################################
    
    general options:
    -f, --fasta=FASTA input file
    -D, --fastaDir=DIR input file dir
    -o, --out=OUTFILE    file result will be written to
    -h, --help  prints this
    [-F, --fasta_list=F1,F2,F3 ]
    [-a, alignment= MSA.fasta ]
    """)
    sys.exit(2)

def mkHeader():
    header = "PepID\tNucID\n"
    return(header)
def getAlignmentHeaders(fasta):
    res =[]
    fp = FastaParser()
    for i in fp.read_fasta(fasta):
        res.append(">"+i[0].strip()) 
    return(res)  
        
def nucFastaSeqAndIds(fastalist, idmapping = None, relevantHeaders=None):
    nucSeqDct = {}
    fp = FastaParser()
    fastalist = fastalist.split(",")
    for fa in fastalist:
        for i in fp.read_fasta(fa):
            tmp = i[0].split(" ")[0].strip()
            try:
                if idmapping[tmp]:
                #nucSeqDct[i[0].strip()]= i[1]
                    if relevantHeaders:
                        if idmapping[tmp] in relevantHeaders:
                            #print("relevant")
                            nucSeqDct[idmapping[tmp]] = i[1]
                        else:
                            pass
                            #print("tmp",tmp)
                    else:
                        nucSeqDct[idmapping[tmp]] = i[1]
            except KeyError:
                pass
    return(nucSeqDct.copy())
def readEnsemblHeaders(fasta):
    fastah = open(fasta,'r')
    for l in fastah:
        m = re.findall('>(.*?) pep.*transcript:(.*?) gene', l)
        if (m):
            #print(m)
            #res = "\t".join([m[0],m[1]])
            yield(m)#res)
def main():
    ###################################
    fasta = None
    alignmentdir = None
    out = None
    alignmentMode = False
    alignment = None
    alignmentHeaders = None
    fastalistMode = False
    fastalist = None
    idmapping = {}
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "f:o:F:a:d:h", ["fasta=","out=","fasta_list=", "alignment=","alignment_dir=","help"])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
    for o, a in opts:
        if o in ("-f", "--fasta"):
            fasta = a  
        elif o in ("-h", "--help"):
            usage()
        elif o in ("-o", "--out"):
            out=a
        elif o in ("-F", "--fasta_list"):
            fastalistMode =True
            fastalist = a
        elif o in ("-a", "--alignment"):
            alignmentMode =True
            alignment = a
        elif o in ("-d", "--alignment_dir=DIR"):
            alignmentMode=True
            alignmentdir =a
        else:
            assert False, "unhandled option"
    
    ########################################
    if fasta is None:
        usage()
    if out is None:
        out=fasta+".pepNuc"
    ########################################
    outh = open(out,"w")
    outh.write(mkHeader())
    for r in readEnsemblHeaders(fasta):
        s = "\t".join([r[0][0],r[0][1]])+"\n"
        outh.write(s)
    outh.close()
   
    filesInDir = []
    if alignmentdir:
        dir = alignmentdir+"/*.fas"
        print (glob.glob(dir))
        for fname in glob.glob(dir):
            print(fname)
            filesInDir.append(fname)
    else:
        filesInDir=[alignment]
    
    for r in readEnsemblHeaders(fasta):
                idmapping[r[0][0]] = r[0][1]
                idmapping[r[0][1]] = r[0][0]

    relevantHeaders =[]
    for alignment in filesInDir:   
        if fastalistMode:
            if alignmentMode:
                tmp = getAlignmentHeaders(alignment)
                tmp = [h.split("|")[1] for h in tmp]
                #for t in tmp:
                #    relevantHeaders.append(t)
                
    dct =nucFastaSeqAndIds(fastalist, idmapping)
    #dct =nucFastaSeqAndIds(fastalist, idmapping, relevantHeaders)
    for alignment in filesInDir:  
        if alignmentMode:
            alignmentHeaders = getAlignmentHeaders(alignment)
            print(alignmentHeaders)
            alignout = open(alignment+".nuc",'w')
            for a in alignmentHeaders:
                tmp = a.split("|")[1]
                try:
                    alignout.write(a+"\n"+dct[a]+"\n")
                except KeyError:
                    try:
                        alignout.write(a+"\n"+dct[tmp]+"\n")
                    except KeyError as k:
                        print("ERROR",k)
                        
            alignout.close()
if __name__ == "__main__": 
    main()

