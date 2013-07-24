import getopt
import imp
import sys
import os.path

gffFastaTools = imp.load_source('gffFastaTools', '../BioHelpers/gffFastaTools.py')
from gffFastaTools import *
VERBOSE = False
#####################################
# [last update 3/29/13 by JM]       #
# DEV                               #
#####################################

def usage():
    print ("""
    ##########################################################
    # python compare_loc_scythe.py -l file.loc  -s res.fasta #
    ##########################################################
    
    general options:
    
    -l, --loc=FILE.LOC 
 
    -s, --scythe_result=RESULT.FASTA    
                single output FASTA file from scythe run
    
    -v, --verbose 
    -h, --help  prints this
    
    advanced options:
    
    -p, --drop_prefix=SEPARATOR    
                 ignore everything before SEPARATOR (incl)
    -t, --table  print tabular summary
    """)
    sys.exit(2)
    
def annoy(*args):
    if VERBOSE:
        print("# ", *args)
    else:
        pass

def main():
    global VERBOSE
    loc = None
    scytheFasta = None
    sep = None
    table = False
    ###################################
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "l:s:p:thv", ["loc=","scythe_result=","drop_prefix=","table","help","verbose"])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
    for o, a in opts:
        if o in ("-l", "--loc"):
            loc=a
        elif o in ("-s", "--scythe_result"):
            scytheFasta = a
        elif o in ("-t", "--table"):
            table = True
        elif o in ("-h", "--help"):
            usage()
        elif o in ("-v", "--verbose"):
            VERBOSE = True    
        elif o in ("-p", "--drop_prefix"):
            sep = a
        else:
            assert False, "unhandled option"
    
    ########################################
    if loc is None:
        usage()
    if scytheFasta is None:
        usage()
    ########################################
    scytheFastaIDs = set() 
    locDefaultIDs = set()
    locNonDefaultIDs = set()
    locMultiModelLoci = set()
    locMultiModelLociDef = set()
    locMultiModelLociNoDef = set()
    locTotalmRNA = set()
    locTotalGenes = set()
    annoy("reading "+scytheFasta)
    for a in gffFastaTools.FastaParser().read_fasta(scytheFasta):
        tmp = a[0]
        annoy(tmp, scytheFasta)
        if sep:
            try:
                tmp = a[0].split(sep)[1]
            except IndexError as ie:
                annoy("#", a[0])
            annoy("->",tmp, scytheFasta)
        scytheFastaIDs.add(tmp)
    annoy("done with "+scytheFasta)
    
    annoy("reading "+loc)
    for ln in open(loc, 'r'):
        #print(ln)
        ln = ln.rstrip()
        locus, genemodels = ln.split("\t")[0],ln.split("\t")[1:]
        locTotalGenes.add(locus)
        #print(genemodels)
        locDefaultIDs.add(genemodels[0])
        for gm in genemodels:
            locTotalmRNA.add(gm)
        if (len(genemodels) >1):
            locMultiModelLoci.add(locus)
            locMultiModelLociDef.add(genemodels[0])
            for gm in genemodels[1:]:
                locMultiModelLociNoDef.add(gm)
                locNonDefaultIDs.add(gm)
        
    ##########
    numTotalGenes = len(locTotalGenes)
    numTotalGeneModels = len(locTotalmRNA)
    numTotalDefaultModels = len(locDefaultIDs)
    numMultiGenes = len(locMultiModelLoci)
    numMultiGenesDefaultModels = len(locMultiModelLociDef)
    
    numScytheGenes = len(scytheFastaIDs)
    numScytheGenesMulti = len(scytheFastaIDs.intersection(locMultiModelLociNoDef.union(locMultiModelLociDef)))
    assert(numTotalDefaultModels == numTotalGenes)
    assert(numMultiGenesDefaultModels == numMultiGenes)
    numMultiGeneAgree = len(scytheFastaIDs.intersection(locMultiModelLociDef))
    numMultiGeneDisagree = len(scytheFastaIDs.intersection(locMultiModelLociNoDef))
    
    #######
    loci_covered = "{:.0%}".format(len(scytheFastaIDs)/len(locTotalGenes))
    consensus = len(scytheFastaIDs.intersection(locDefaultIDs))
    nonConsensus = len(scytheFastaIDs)-len(scytheFastaIDs.intersection(locDefaultIDs))
    agreeRatio = "{:.0%}".format(consensus/len(scytheFastaIDs))
    disagreeRatio = "{:.0%}".format(nonConsensus/len(scytheFastaIDs))
    
    consensusMulti=len(scytheFastaIDs.intersection(locMultiModelLociDef))
    nonConsensusMulti=len((scytheFastaIDs.intersection(locMultiModelLociNoDef)))
    #print(consensusMulti/len(locMultiModelLoci)*100)#, nonConsensusMulti)
    #print(len(locTotalGenes)-len(locMultiModelLoci))
    ##???
    try:
        agreeRatioMulti = "{:.0%}".format(numMultiGeneAgree/numScytheGenesMulti)
        disagreeRatioMulti = "{:.0%}".format(numMultiGeneDisagree/numScytheGenesMulti)
    except ZeroDivisionError as zde:
        print("Error: No multi model genes in result fasta file. Please check your input.")
        print(zde)
        exit(2)
    
   # print(len(locMultiModelLoci))
   # print(len(locMultiModelLociDef))
   # print(len(locMultiModelLociNoDef), "nodef")

    models_per_gene ="{0:.2f}".format(len(locTotalmRNA)/len(locTotalGenes))
    perc_multi_models =" {:.0%}".format(numMultiGenes/numTotalGenes)
    multiModelLoci = len(locMultiModelLoci)
    multiModelLociPerc = "{:.0%}".format(multiModelLoci/len(locTotalGenes))
    #-->
    multiModelLociCovered = locMultiModelLociDef.union(locMultiModelLociNoDef).intersection(scytheFastaIDs)
    
    fileNameLoc= os.path.basename(loc)
    fileNameFa= os.path.basename(scytheFasta)
    if table:
        print("# "+fileNameFa+" "+fileNameLoc)
        print("\t".join([str(len(locTotalGenes)), "total genes (loc)"]))
        print("\t".join([str(numMultiGenes), "multi loci genes (loc)"]))
        print("\t".join([str(perc_multi_models), "multi loci genes % (loc)"]))
        print("\t".join([str(len(locTotalmRNA)), "total models (loc)"]))
        print("\t".join([str(len(locDefaultIDs)), "default models (loc)"]))
        print("\t".join([str(len(locNonDefaultIDs)), "non-default models (loc)"]))
        print("\t".join([models_per_gene, "models per gene (loc)"]))
        print("\t".join([str(len(scytheFastaIDs)), "genes processed (fasta)"]))
        print(loci_covered+"\tpercent genes covered (fasta vs loc)")
        print(str(numScytheGenesMulti)+"\tmulti model genes covered (fasta vs loc)")
        print("{:.0%}".format(numScytheGenesMulti/numMultiGenes)+"\tmulti model genes covered % (fasta vs loc)")
        print("\t".join([str(consensus),"consensus total (fasta)" ]))
        print("\t".join([str(nonConsensus),"non consensus total  (fasta)" ]))
        print("\t".join([agreeRatio,"agree ratio (fasta) total" ]))
        print("\t".join([disagreeRatio,"disagree ratio (fasta) total" ]))
        print("\t".join([agreeRatioMulti,"agree ratio (fasta) multi model genes" ]))
        print("\t".join([disagreeRatioMulti,"disagree ratio (fasta) multi model genes" ]))
    else:
        tmp= fileNameLoc + " has a total of "+str(numTotalGenes)+ " gene loci\n"
        tmp+=fileNameFa + " coveres "+ str(numScytheGenes)+ " of those (={0:.0%}".format(numScytheGenes/numTotalGenes) +").\n"
        tmp+=str(numMultiGenes)+" of the total gene loci have more than one gene model (={0:.0%}".format(numMultiGenes/numTotalGenes)+").\n"
        tmp+=fileNameFa + " coveres "+ str(numScytheGenesMulti)+ " of those (={0:.0%}".format(numScytheGenesMulti/numMultiGenes)+").\n"
        tmp+="Among the "+str(numScytheGenesMulti)+" covered multi model loci, Scythe agrees in "+str(numMultiGeneAgree)+" cases with the default (="+agreeRatioMulti+")\n"
        tmp+="and disagrees in "+str(numMultiGeneDisagree)+" cases  (="+disagreeRatioMulti+").\n"
        
        print(tmp)

if __name__ == "__main__": 
    main()

