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
#
#####################################
#todo: throw error if id not in fasta
"""
rm sequence unavailable  and save in new folder
"""

def usage():
    print ("""
    ##########################################################################
    # python 
    #########################################################################
    
    general options:

    -f, --fasta	                      fasta file 
    -o, --outdir=OUTPUTDIR            write files to OUTPUTDIR instead of ./
    
    -v, --verbose 
    -h, --help  prints this
    
    """)
    sys.exit(2)
    

def main():
    fasta = None
    outdir = "./"
    sep = None
    ###################################
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "f:o:hv", ["fasta=","outdir=","help","verbose"])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
    for o, a in opts:
        if o in ("-f", "--fasta"):
            fasta = a
        elif o in ("-o", "--outdir"):
            outdir = a
        elif o in ("-h", "--help"):
            usage()
        elif o in ("-v", "--verbose"):
            VERBOSE = True    
        else:
            assert False, "unhandled option"
    
    ########################################
    if fasta is None:
        usage()
    if outdir:
        print(outdir) 
        if not os.path.exists(outdir):
            os.makedirs(outdir)
    ########################################
    outfile = open(outdir+"/"+os.path.basename(fasta)+".available.fa","w")
    for a in gffFastaTools.FastaParser().read_fasta(fasta):
        tmp = a[0]
        if not "unavailable" in a[1]:
            outfile.write(">"+tmp+"\n")
            outfile.write(a[1]+"\n")
        else:
            print("## ",a[0],a[1])
        
 
    outfile.close()
if __name__ == "__main__": 
    main()

