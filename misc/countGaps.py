import getopt
import re
import sys



def usage():
    print ("""
    ##############################
    #
    ##############################
    -i, --infile=ALIGNMENT
    -h, --help
    """)
    sys.exit(2)

def readNexus(nexus, stfu=True):
    alignment_length = None
    number_species = None
    number_gaps = None
    infile = open(nexus,'r')
    tmp = infile.read()
    res = re.findall("ntax=(.*)nchar=(.*);", tmp)
    number_species=int(res[0][0].strip())
    alignment_length=int(res[0][1].strip())
    alignment = re.findall("matrix(.*);", tmp, re.DOTALL)[0]
    number_gaps=alignment.count("-")
    perc_gaps=number_gaps/(alignment_length*number_species)*100
    if not stfu:
        print(number_gaps)
    #print(number_species,alignment_length)
    return(perc_gaps)

def readFasta(fasta,stfu=True):
    alignment = ""
    alignment_length = 0
    number_species = 0
    number_gaps = None
    infile = open(fasta,'r')
    tmp = infile.read()
    for ln in tmp:
        if ln.startswith(">"):
            number_species+=1
        else:
            ln = ln.strip()
            #print(len(ln))
            alignment += ln
    alignment_length= len(alignment)
    number_gaps = alignment.count("-")
    perc_gaps = number_gaps/alignment_length*100
    if not stfu:
        print(number_gaps)
    return(perc_gaps)
            
def guessFormat(infile, stfu=True):
    if not stfu:
        print("guessing format: believe it's... ", end="")
    infile = open(infile,'r')
    ln = infile.readline()
    infile.close()
    if ln.startswith("#NEXUS"):
        if not stfu:
            print("nexus.\n")
        return ("nexus")
    elif ln.startswith(">"):
        if not stfu:
            print("fasta.\n")
        return("fasta")
    else:
        sys.stderr.write("unknown format\n")
        return(None)
   
    
def main():

    ###################################
    infile = None
    format = None
    stfu = True
    verbose=False
    ###################################
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "i:vVh", ["infile=","verbose","Verbose","help"])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
    for o, a in opts:
        if o in ("-i", "--infile"):
            infile = a
        elif o in ("-h", "--help"):
            usage()
        elif o in ("-v", "--verbose"):
            stfu=False
        elif o in ("-V", "--Verbose"):
            verbose=True
        else:
            assert False, "unhandled option"

    if not infile:
        usage()
    else:
        if format == None:
            format = guessFormat(infile,stfu)
        #check format
        if format=="nexus":
            pg = readNexus(infile,stfu)
        elif format == "fasta":
            pg = readFasta(infile,stfu)
        else:
            exit(1)
        if verbose:
            print(infile+"\t"+"{:.2f}".format(pg))
        else:
            print(pg)
####################################
if __name__ == "__main__":
    main()
####################################