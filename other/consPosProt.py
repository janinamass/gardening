import sys
import getopt

##########
# dev
# J.Mass
# mod: 26Feb14
##########

class Alignment(object):
    def __init__(self, id=None, fasta = None, members = []):
        self.id = id
        self.fasta = fasta
        self.members = []
        self.foregroundConserved = []
        self.scores = []
    def calcForegroundConserved(self,ignoreGaps = False, strictly_different = False):
        cons = []
        alignmentLength = len(self.members[0].sequence)
        #check if all labelled seq have the same character
        foreground = [m for m in self.members if m.isForeground == True]
        print("FORE",foreground, self.members)
        background = [m for m in self.members if m.isForeground == False]
        print("BACK",foreground, self.members[3].isForeground)
        for i in range(0, alignmentLength):
            #list comp
            fgpos = [fg.sequence[i] for fg in foreground]
            #print("FPO",fgpos)
            #just check if all the same:
            tmpcons = True
            for f in fgpos:
                #print(f,fgpos)
                if len([tmp for tmp in fgpos if tmp == f]) != len(fgpos):
                    #print("NC")
                    #not conserved
                    tmpcons=False
                if ignoreGaps==True:
                    if fgpos[0]=="-":
                        tmpcons = False
            #now check if background is different:
            cons_diff = False
            if tmpcons==True:
                bgpos = [bg.sequence[i] for bg in background]
                for b in bgpos:
                    # strict: no background sequence is allowed to have 
                    # the foreground character
                    if b != fgpos[0]:
                            cons_diff = True
                    if strictly_different:
                        if len([f for f in fgpos if f == b]) >0:
                            cons_diff = False
                            break
            cons.append(cons_diff)
        return(cons)
            
    def calcScores(self, ignoreGaps = False, strictly_different = False):
        cons = self.calcForegroundConserved(ignoreGaps,strictly_different )
        numMem = len(self.members)
        alignmentLength = len(self.members[0].sequence)
        print("Alignment length ",alignmentLength)
        for i in range(0,alignmentLength):
            s = ""
            pos = ""
            for j in range(0,numMem):
                s+= "" + self.members[j].sequence[i]
                pos +=self.members[j].sequence[i]
            ps = calcPosScore(pos)
            pss = ""
            for k,v in sorted(ps.items(), key=lambda x:x[1], reverse=True):
                pss += ":".join((k,str(v)))
                pss +=" "
            print(s,pss, cons[i],i)
def calcPosScore(str):
    posScore = {}
    for s in str:
        posScore[s] = len([t for t in str if t == s])/len(str)
    return(posScore)
class Sequence():
    def __init__(self, id = "", sequence = None, isForeground = False):
        self.id = id
        self.sequence = sequence
        self.isForeground = isForeground
    
    def setForeground(self, bool = True):
        self.isForeground = bool

class SeqTranslator(object):
    RNAmap = {
           "UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
           "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
           "UAU":"Y", "UAC":"Y", "UAA":"*", "UAG":"*",
           "UGU":"C", "UGC":"C", "UGA":"*", "UGG":"W",
           "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
           "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
           "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
           "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
           "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
           "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
           "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
           "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
           "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
           "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
           "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
           "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"
    }
    DNAmap = {
           "TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
           "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
           "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",
           "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",
           "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
           "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
           "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
           "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
           "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
           "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
           "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
           "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
           "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
           "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
           "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
           "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"
    }
    def triplets(self, string, frameshift = 0):
        #chop string into 3 char slices 
        for i in range(0, int(len(string)/3)):
            #print(string[i*3+frameshift:i*3+3+frameshift])
            yield string[i*3+frameshift:i*3+3+frameshift]
    def triplets1(self, string, frameshift = 0):
        if frameshift == 3:
            frameshift=0
        #chop string into 3 char slices 
        for i in range(0, int(len(string)/3)):
            #print(string[i*3+frameshift:i*3+3+frameshift])
            if frameshift==0:
                yield string[i*3+frameshift:i*3+3+frameshift]
            elif frameshift==1:
                yield string[1+ i*3+frameshift:1+ i*3+3+frameshift]
            elif frameshift ==2:
                yield string[2+ i*3+frameshift:2+ i*3+3+frameshift]
            else:
                print("\t??")
    def dna2prot(self, string, frameshift=0):
        res = ""
        for a in self.triplets(string, frameshift):
            try:
                res+=self.DNAmap[a]
            except KeyError as e:
                print(e, "But don't panic.")
                res+=""
        return res



class FastaParser(object):     
    def read_fasta(self, fasta, delim = None, asID = 0):
        """read from fasta fasta file 'fasta' 
        and split sequence id at 'delim' (if set)\n
        example:\n
        >idpart1|idpart2\n
        ATGTGA\n
        and 'delim="|"' returns ("idpart1", "ATGTGA")
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
        
###########################################
def usage():
    print ("""
    ######################################
    #  consPos.py       [JM]   
    ######################################
    usage:
        consPos.py -f multifasta alignment -l foreground labels
    options:
        -f, --fasta    multifasta alignment    align.fa
        -l, --foreground_labels                LAB1,LAB2,...
        -i, --ignore_gaps    gaps do not count as shared
        -s, --strictly_different  no background seq is allowed to share fg character at given pos
        -d, --dir=DIR   TODO...
        -h, --help      prints this
        -T, --translate translate into protein
        
    """)
    sys.exit(2)
############################################
def main():
  fasta = None
  fgLabels = None
  ignoreGaps = False
  strictDiff = False
  toProt = False
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "f:l:siht", ["fasta=","foreground_labels=","strictly_different","ignore_gaps","help","translate"])
  except getopt.GetoptError as err:
        print (str(err))
        usage()
  for o, a in opts:
        if o in ("-f", "--fasta"):
            fasta = a
        elif o in ("-l","--foreground_labels"):
            fgLabels = a
        elif o in ("-i","--ignore_gaps"):
            ignoreGaps = True
        elif o in ("-s","--strictly_different"):
            strictDiff = True
        elif o in ("-h","--help"):
            usage()
        elif o in ("-t","--translate"):
            toProt = True
        else:
            assert False, "unhandled option"
  if not fasta:
      usage()
  if not fgLabels:
      usage()
  else:
      if toProt:
          st = SeqTranslator()
          fp = FastaParser()
          newAlignment = Alignment(id = fasta)
          if fgLabels:
              fgLabels = fgLabels.split(",")
              fp = FastaParser()
          for f in fp.read_fasta(fasta):
              if fgLabels:
                  if f[0].startswith(tuple(fgLabels)):
                      tmpIsForeground = True
                      print("#",f[0],"is foreground")
                  else:
                      tmpIsForeground = False
              else:
                  tmpIsForeground = False
              newSeq = Sequence(id = f[0], sequence = st.dna2prot(f[1],0), isForeground = tmpIsForeground)
          #print(newSeq.isForeground,"@@", tmpIsForeground)
              newAlignment.members.append(newSeq)
          
          
          
      else:
          newAlignment = Alignment(id = fasta)
          if fgLabels:
              fgLabels = fgLabels.split(",")
          fp = FastaParser()
          for f in fp.read_fasta(fasta):
              if fgLabels:
                  if f[0].startswith(tuple(fgLabels)):
                      tmpIsForeground = True
                      print("#",f[0],"is foreground")
                  else:
                       tmpIsForeground = False
              else:
                  tmpIsForeground = False
              newSeq = Sequence(id = f[0], sequence = f[1], isForeground = tmpIsForeground)
              #print(newSeq.isForeground,"@@", tmpIsForeground)
              newAlignment.members.append(newSeq)
      print(strictDiff,"SD", ignoreGaps, "IG")
      newAlignment.calcScores(ignoreGaps = ignoreGaps,strictly_different = strictDiff)


#############################################
if __name__ == "__main__":
    main()
