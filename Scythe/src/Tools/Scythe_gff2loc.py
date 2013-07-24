#!/usr/bin/python
import re,sys, getopt

#####################################
# last update 03/31/2013 by J. Mass #
# version = '0.1'                   #
#####################################

def usage():
    print ("""
    ##############################
    #  Scythe_gff2loc.py v0.1    #
    ##############################
    
    -f, --file=gff3_FILE     (tested w/ _gene.gff3 from phytozome)
    -o, --output=FILE        output file [default: gff3_FILE.loc]
    -h, --help               prints this
    """)
    sys.exit(2)

def read_gff2loc(infile, outfile):
	infile = open(infile, 'r')
	outfile = open(outfile, 'w')
	loci = {}
	longest = {}
	rawstr = r"""(Name=)(.*);pacid.*(longest=)(.*);(Parent=)(.*)"""
	cnt = 0
	for ln in infile:
		s =ln
		m = re.findall(rawstr, s)
		if len(m)  >0:
			name =  m[0][1]
			isLongest = m[0][3]
			parent = m[0][5]
			if isLongest == str(1):
				if parent in longest:
					print("#Warning "+parent+" has more than one default model\nCheck your gff -> ", longest[parent], name)
				longest[parent]=name #longest will be printed to 2nd col
			elif isLongest == str(0):
				if parent in loci:
					loci[parent].append(name)
				else:
					loci[parent]=[name]
	
	s_def =  sorted(longest.keys())
	for k_def in s_def:
		try:
			outfile.write(k_def+"\t"+longest[k_def]+"\t"+"\t".join(loci[k_def])+"\n")
		except KeyError as ke:
			outfile.write(k_def+"\t"+longest[k_def]+"\n")
		if k_def in loci: 
			del loci[k_def]
	s = sorted(loci.keys())
	for k in s:
		try:
			outfile.write(k+"\t"+longest[k]+"\t"+"\t".join(loci[k])+"\n")
		except KeyError as ke:
			print("#Warning "+k+" has no default model\n")
			outfile.write(k+"\t"+"\t".join(loci[k])+"\n")
	return loci

###################################
outfile = None
infile = None
###################################
try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "f:ho:", ["file=","help", "output="])
except getopt.GetoptError as err:
    print (str(err))
    usage()
for o, a in opts:
    if o in ("-f", "--file"):
        infile=a
    elif o in ("-h", "--help"):
        usage()
    elif o in ("-o", "--output"):
        outfile = a
    else:
        assert False, "unhandled option"

########################################
if infile is None:
	usage()
if outfile is None:
	outfile = infile+".loc"
########################################
read_gff2loc(infile, outfile)
