#!/usr/bin/env python
import sys
import getopt
#########################
# last update:
# Mo 19 Mai 2014 12:49:54 CEST
# [JMass]
#########################
def usage():
    print ("""
    ###########################
    #  scythe_loc_tsv.py      #
    ###########################
    -f, --file=ENSEMBLBioMart.tsv
                              format: 1st column: gene id, 2nd column:transcript id,
                                      3rd column: peptide id, 4th column: cds length;
                                      gene ids can occur multiple times
    -c, --custom=COLx,COLy,COLz,...   COLi in ["gene","transcript", "protein", "length"]
                                      Use this if your file is different from the described biomart output.
                                      "cds_length" is optional but recommended, at least one of
                                      ["transcript", "protein"] need to be included
    -o, --output=FILE         output file [default: ENSEMBLEBioMart.tsv.loc]
    -h, --help                prints this
    -H, --HELP                show help on format
    #----------------------------------#
    """)
    sys.exit(2)
def formatHelp():
    print("""
    #------------ loc output format ---------------------------#

    LOCUS0\tTRANSCRIPT0_0\tTRANSCRIPT0_1\t...\tTRANSCRIPT0_n
    LOCUS1\tTRANSCRIPT1_0\t...\tTRANSCRIPT1_m
    .
    .
    .
    LOCUSk\tTRANSCRIPTk_0\t...TRANSCRIPTk_l

    #----------------------------------------------------------#

    """)
    howTo()

def howTo():
    """
    http://www.ensembl.org/biomart/
    <?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
    <Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
    <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
    <Filter name = "biotype" value = "protein_coding"/>
    <Attribute name = "ensembl_gene_id" />
    <Attribute name = "ensembl_transcript_id" />
    <Attribute name = "ensembl_peptide_id" />
    <Attribute name = "cds_length" />
    </Dataset>
    </Query>
    """
class WrongNumberOfColumnsException(Exception):
    pass
def checkTsv(custom=None, infile = None):
    print(custom)
    with open(infile,'r') as f:
        h = f.readline().strip()
    print(h)
    header = h.split("\t")
    f = lambda x: x.lower().replace(" ","")
    header_tmp = [f(h) for h in header]
    if len(custom)!=len(header):
        raise WrongNumberOfColumnsException("Lengths of header and file don't match! ", custom, header)

    for c,h in zip(custom,header_tmp):
        if not c in h:
            raise Warning("Something wrong with the headers:"+",".join(header))

def writeLoc(genesDct, outfile):
    f = lambda x: "\t".join(genesDct[x])
    tmp = [g+"\t"+f(g) for g in sorted(genesDct.keys())]
    res = "\n".join(tmp)
    with open(outfile,'w') as out:
        out.write(res)

def readTsv(infile=None):
    genes = set()
    genest = {}
    genesp = {}
    maxlen = {}
    gene_index = None
    transcript_index = None
    protein_index = None
    length_index = None
    #first line should be header
    with open(infile,'r') as f:
        headers = f.readline().strip().split("\t")
        headers = [h.lower().replace(" ","") for h in headers]
        for i,h in enumerate(headers):
            if "gene" in h:
                gene_index = i
            elif "transcript" in h:
                transcript_index = i
            elif "protein" in h or "peptide" in h:
                protein_index = i
            elif "length" in h:
                length_index = i
        for l in f:
            tmp = l.strip().split("\t")
            tmpg = tmp[gene_index]
            if not tmpg in genes:
                genes.add(tmpg)
                if transcript_index:
                    genest[tmpg] = [tmp[transcript_index]]
                if protein_index:
                    try:
                        genesp[tmpg] = [tmp[protein_index]]
                    except IndexError as e:
                        sys.stderr.write(str(tmp)+" ignored (missing protein id field)\n")
                if length_index:
                    try:
                        maxlen[tmpg] = tmp[length_index]
                    except IndexError as e:
                        sys.stderr.write(str(tmp)+" ignored (missing length field)\n")
                        maxlen[tmpg] = 0
            else:
                if length_index:
                    try:
                        if  tmp[length_index] >  maxlen[tmpg]:
                            maxlen[tmpg] = tmp[length_index]
                            if transcript_index:
                                genest[tmpg].insert(0, tmp[transcript_index])
                            if protein_index:
                                try:
                                    genesp[tmpg].insert(0,tmp[protein_index])
                                except KeyError as e:
                                    genesp[tmpg] = [tmp[protein_index]]
                        else:
                            if transcript_index:
                                genest[tmpg].append(tmp[transcript_index])
                            if protein_index:
                                genesp[tmpg].append(tmp[protein_index])
                    except IndexError as e:
                        sys.stderr.write(str(tmp)+" ignored (missing length field)\n")
                else:
                    if transcript_index:
                        genest[tmpg].append(tmp[transcript_index])
                    if protein_index:
                        genesp[tmpg].append(tmp[protein_index])

    return(genest.copy(),genesp.copy())
def readEnsemblLoc(ensembleLoc, outfile, lenDct=None):
    out=open(outfile,"w")
    infile = open(ensembleLoc,"r")
    numLoc = 0
    numTr = 0
    locDct = {}
    maxDct = {} #max length of transcript at locus
    longestTr = {}
    missing = False
    for l in infile:
        if "Ensembl" in l:
            continue
        l = l.strip()
        tmp = l.split("\t")
        if tmp[0] not in locDct:
            numLoc+=1
            numTr+=1
            locDct[tmp[0]]=[tmp[1]]
        else:
            numTr+=1
            locDct[tmp[0]].append(tmp[1])

    if lenDct:
        for l,tr in locDct.items():
            try:
                maxDct[l] = max([lenDct[i] for i in tr])
            except KeyError as ke:
                tmp=[]
                missing = True
                for i in tr:
                    if i in lenDct:
                        tmp.append(i)
                    else:
                        lenDct[i]=0
                        tmp.append(i)
                maxDct[l] = max([lenDct[i] for i in tmp])
                if not maxDct[l]:
                    #last resort
                    maxDct[l]=tr[0]
            longestTr[l]=[t for t in tr if lenDct[t] == maxDct[l]]
            #print(longestTr[l])
            longestTr[l] = longestTr[l][0]
        locDct2 = locDct.copy()
        for loc,tr in locDct2.items():#rm
            locDct[loc] = [t for t in tr if t != longestTr[loc]]
        for loc,tr in locDct.items():
            out.write(loc+"\t"+longestTr[loc]+"\t"+"\t".join(tr)+"\n")
        if missing:
            print("# There were missing lengths\n")

    else: #no lengths given
        for loc,tr in locDct.items():
            out.write(loc+"\t"+"\t".join(tr)+"\n")
    out.close()
    printInfo(numLoc=numLoc, numTr=numTr, outfile=outfile, verb=VERB)


def main():
    custom = None
    outfile = None
    outfilep = None
    in_tsv = None
    default_header = ["gene","transcript","protein", "length"]
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "f:hHo:c:", ["file=","help","HELP","output=","custom="])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
    for o, a in opts:
        if o in ("-f", "--file"):
            in_tsv = a
        elif o in ("-h", "--help"):
            usage()
        elif o in ("-H", "--HELP"):
            formatHelp()
        elif o in ("-o", "--output"):
            outfile = a
            outfilep = a+"p"
        else:
            assert False, "unhandled option"

    if not in_tsv:
        usage()
    if not outfile:
            outfile = in_tsv+".loc"
            outfilep = in_tsv+".locp"
    if custom:
        header = custom.split(",")
        valid = [h for h in header if h in default_header]
        if not default_header[0] in valid:
            usage()
        elif not default_header[1] in valid and not default_header[2] in valid:
            usage()
    else:
        valid = default_header
    try:
        checkTsv(custom = valid , infile = in_tsv)
    except WrongNumberOfColumnsException as e:
        print(str(e))
        usage()

    genest, genesp = readTsv(in_tsv)
    if genest:
        writeLoc(genest, outfile)
    if genesp:
        writeLoc(genesp, outfilep)
########################################
if __name__ == "__main__":
    main()
