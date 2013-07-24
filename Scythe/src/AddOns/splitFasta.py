import getopt
import os
import glob
import imp
import sys
path = os.path.join(os.path.dirname(__file__), 'BioHelpers/*.py')
#-- import from BioHelpers --#
for infile in glob.glob(path):
    fullbasename = os.path.basename(infile)
    basename = fullbasename[:-3]
    imp.load_source(basename, infile)
from gffFastaTools import *
#-- usage --#
def usage():
    print ("""
    ######################################
    #  splitFasta.py   [03/25/2013 JM]   #
    ######################################
  usage:
     splitFasta.py -i INFILE -n NUM_SEQ
  options:
  -i INFILE, --infile=INFILE  input file in FASTA format
  -n NUM_SEQ, --num_seq=NUM_SEQ  infile will be split every NUM_SEQ sequences
    """)
    sys.exit(2)
#-----------#



def main():
  infile = None
  outfile = None
  num = 5000
  cnt = 0
  suff = 0
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "i:n:h", ["infile=", "num_seq=","help"])
  except getopt.GetoptError as err:
        print (str(err))
        usage()
  for o, a in opts:
        if o in ("-i", "--infile"):
            infile = a
        elif o in ("-n", "--num_files"):
            num = int(a)
        elif o in ("-h", "--help"):
            usage()
        else:
            assert False, "unhandled option"
    
  if not (infile):
        usage()             

  fh = FastaParser().read_fasta(infile)

  tmp = open(infile+"."+str(suff), 'w')
  for h in fh:
      cnt+=1
      tmp.write(">"+h[0]+"\n"+h[1]+"\n")
      if(cnt%num==0):
          print(cnt, h[0])
          suff+=1
          tmp.close()
          tmp = open(infile+"."+str(suff), 'w')

          
  tmp.close()

if __name__ == "__main__":
    main()