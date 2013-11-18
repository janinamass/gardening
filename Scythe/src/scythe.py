#!/usr/bin/python3
import imp
import sys
import subprocess
import os
import glob
import getopt
from itertools import chain
import configparser

path = os.path.join(os.path.dirname(__file__), 'BioHelpers/*.py')
#-- import from BioHelpers --#
for infile in glob.glob(path):
    fullbasename = os.path.basename(infile)
    basename = fullbasename[:-3]
    imp.load_source(basename, infile)
from gffFastaTools import *
from clusteringTools import *
from scythe_classes import *
#----/import------------------#
logo = """
          _____            __  __       
         / ___/_______  __/ /_/ /_  ___ 
         \__ \/ ___/ / / / __/ __ \/ _ \\
        ___/ / /__/ /_/ / /_/ / / /  __/
       /____/\___/\__, /\__/_/ /_/\___/ 
                 /____/"""
#-- usage --#
def usage():
    """Print help."""
    print ("""
    ######################################
    #  scythe.py v0.5b   [10/31/2013 JM] #
    ######################################
  usage:
     scythe.py -i DIR -g .grpFILE --cleanup
     
  usage with configuration file:
     scythe.py --config configuration.scy

  general options:  
    -C, --config                     use configuration file instead of
                                     command line parameters
    -c, --cleanup                    remove temporary files when done
    -h, --help                       prints this
    -i, --in_dir=DIR                 folder w/ subfolder "fa" and "loc"
    -m, --global_max                 always start with best pair
    -M, --global_sum		     calculate sum of best scores and choose gene model accordingly
    -o, --out_dir=DIR                output directory [default:./]
    -v, --verbose                    be wordy

  alignment options:
     -O, --gap_open=FLOAT           needleall gap opening cost [default 10]
     -E, --gap_extend=FLOAT         needleall gap extension cost  

  fasta options:
    -d, --delim=STRING               split fasta headers at STRING
    -a, --asID=INT                   use INTth part of fasta header as transcript-ID
                                     (default:0)
  debug options:
    -s, --stop_after=NUM             stop after NUM groups
  
  further help:
    Please see documentation in the 'doc'-directory.
    """)
    sys.exit(2)
#-----/usage------#

#----parse configuration----#
def parseConfig(pathconfig):
    VERBOSE =False
    CF_MODE = "Mode"
    CF_MODE_use_ensembl = "use_ensembl_api"
    CF_MODE_use_local_files = "use_local_files"
    CF_PATHS = "Paths"
    CF_PATHS_fasta_directory = "fasta_directory"
    CF_PATHS_loc_directory = "loc_directory"
    CF_PATHS_grp_file = "grp_file"
    CF_PATHS_output_directory = "output_directory"
    CF_OUTPUT ="Output"
    CF_OUTPUT_attach_output_prefix="attach_output_prefix"
    CF_OUTPUT_output_prefix="output_prefix"
    CF_OUTPUT_output_orthogroups = "output_orthogroups"
    CF_OUTPUT_output_species_fasta = "output_species_fasta"
    CF_OUTPUT_merge_species_fasta_with_defaults = "merge_species_fasta_with_defaults"
    CF_CLEANUP = "Cleanup"
    CF_CLEANUP_clean_up_directories = "clean_up_directories"
    CF_RUN="Run_options"
    CF_RUN_max_threads ="max_threads"
    CF_RUN_split_input="split_input"
    CF_RUN_split_each="split_each"
    CF_RUN_use_seqan="use_seqan"
    CF_PENALTIES = "Penalties"
    CF_PENALTIES_gap_open_cost = "gap_open_cost"
    CF_PENALTIES_gap_extend_cost="gap_extend_cost"
    CF_PENALTIES_substitution_matrix="substitution_matrix"
    CF_ALGORITHM = "Algorithm"
    CF_ALGORITHM_use_global_max ="use_global_max"
    CF_ALGORITHM_use_default="use_default"
    CF_ALGORITHM_use_global_sum="use_global_sum"
    CF_PARALOGS="Paralogs"
    CF_PARALOGS_include_paralogs = "include_paralogs"

    
# SECTIONS
# CF_MODE: CF_MODE_use_ensembl,CF_MODE_use_local_files
# CF_PATHS: CF_PATHS_fasta_directory, CF_PATHS_loc_directory,CF_PATHS_grp_file,
#       CF_PATHS_output_directory
# CF_OUTPUT: CF_OUTPUT_attach_output_prefix, CF_OUTPUT_output_prefix,
#     CF_OUTPUT_output_orthogroups, CF_OUTPUT_output_species_fasta
#     CF_OUTPUT_merge_species_fasta_with_defaults
# CF_CLEANUP: CF_CLEANUP_clean_up_directories
# CF_RUN: CF_RUN_max_threads, CF_RUN_split_input,CF_RUN_split_each,
#    CF_RUN_use_seqan
# CF_PENALTIES: CF_PENALTIES_gap_open_cost, CF_PENALTIES_gap_extend_cost
#    CF_PENALTIES_substitution_matrix
# CF_ALGORITHM: CF_ALGORITHM_use_global_max,CF_ALGORITHM_use_default,
#     CF_ALGORITHM_use_global_sum 
# CF_PARALOGS:[CF_PARALOGS_include_paralogs  ]:
  
    config = configparser.ConfigParser()
    print(config.sections())
    print(pathconfig)
    config.read(pathconfig)
    print(config.sections())
    print(config.sections())
    for c in config:
        print(c)
   # print(config)
    #global VERBOSE
    global GLOBMAX
    global GLOBSUM 
    #VERBOSE = None
    if config.get(CF_ALGORITHM,CF_ALGORITHM_use_global_max)!="yes":
        GLOBMAX = False
    else:
         GLOBMAX = True
    if config.get(CF_ALGORITHM,CF_ALGORITHM_use_global_sum)!="yes":
        GLOBSUM = False
    else:
        GLOBSUM = True
    
    if config.get(CF_CLEANUP,CF_CLEANUP_clean_up_directories) !="yes":
        cleanUp = False
    else:
        cleanUp = True
    
    
    groups= config.get(CF_PATHS,CF_PATHS_grp_file)
    namesList = None
    faDir = config.get(CF_PATHS,CF_PATHS_fasta_directory)
    inDir = faDir
    outDir = config.get(CF_PATHS,CF_PATHS_output_directory)
    locDir = config.get(CF_PATHS,CF_PATHS_loc_directory)
    fastaList = os.listdir(faDir)
    #gffList = None
    delim = None
    asID = 0
    stopAfter = False
    gapOpen= config.get(CF_PENALTIES,CF_PENALTIES_gap_open_cost)
    gapExtend =config.get(CF_PENALTIES,CF_PENALTIES_gap_extend_cost)
    faFileList = os.listdir(faDir)
    namesList = os.listdir(faDir)
    namesList = [n[0:3] for n in namesList]
    
    
    print(groups)
    print(namesList)
    print(gapOpen,gapExtend )
    print(faDir, faFileList)
    print("LODIR", locDir)
    runScythe(groups=groups, delim=delim, 
              asID=asID, faFileList=faFileList, 
              namesList=namesList, cleanUp=cleanUp, 
              stopAfter=stopAfter, inDir=inDir, outDir=outDir,
              gapOpen=gapOpen, gapExtend=gapExtend,
              locDir=locDir,faDir=faDir)


#----/parse configuration----#



#-- be verbose --#
def annoy(*args):
    """Prints if VERBOSE is set to True. Use for debugging."""
    if VERBOSE==True:
        print(*args)
#----------------#
def append_add(dct, key, val, unique = False):
    """Append new (unique) value to the list of 'key'. 
    If 'key' does not exist, create new list.
    """
    if key in dct:
        if unique:
            if val in dct[key]:
                pass
            else:
                dct[key].append(val)
        else:
            dct[key].append(val)
    else:
        dct[key] = [val]

def initCollProc(avd, seqDct):
    """Initialize sets of species. """
    species2id = {}
    unprocessed = set()
    processed =set()
    coll = set()
    uncoll = set()
    """Set of species that have already been processed"""
    for fkey in avd.keys():
        uncoll.add(fkey) #not yet collected
        unprocessed.add(seqDct[fkey].species) #species on to do list
        append_add(species2id, seqDct[fkey].species, fkey)
        for skey in avd[fkey].keys():
            uncoll.add(skey)
            unprocessed.add(seqDct[skey].species)
            append_add(species2id, seqDct[skey].species, skey, unique=True)
    return(processed, unprocessed, coll, uncoll, species2id)

def checkSingle(avd, seqDct, singles, all):
    species2id = {}
    unprocessed = set()
    processed =set()
    coll = set()
    uncoll = set()
    for fkey in avd.keys():
        tmpspec = seqDct[fkey].species
        uncoll.add(fkey)
        unprocessed.add(tmpspec)
        append_add(species2id, tmpspec, fkey)
        for skey in avd[fkey].keys():
            tmpspec = seqDct[skey].species
            uncoll.add(skey)
            unprocessed.add(tmpspec)
            append_add(species2id, tmpspec, skey, unique=True)        
    # adding to coll, uncoll, proc, unproc
    for fkey in avd.keys():
        # fkey key for first dimension
        tmpspec = seqDct[fkey].species
        if fkey in singles and tmpspec in unprocessed:
            uncoll.remove(fkey)
            coll.add(fkey)
        # if fkey in singles and tmpspec not in processed:
            unprocessed.remove(tmpspec)
            processed.add(tmpspec)    
        skl=list(avd[fkey].keys())
        
        for o in [x for x in skl if x in singles and seqDct[x].species in unprocessed]:
            unprocessed.remove(seqDct[o].species)
            processed.add(seqDct[o].species)
            uncoll.remove(o)
            coll.add(o)
    return(processed, unprocessed, coll, uncoll, species2id)


def getPairwiseAsTuples(avd):
    annoy("pairwise")
    firstdim = avd.keys()
    secdim = []
    pwdict={}
    listofkeys =[]
    for f in firstdim:
         if f not in listofkeys:
             listofkeys.append(f)
         for s in avd[f].keys():
            if s not in listofkeys: 
                 listofkeys.append(s)
            pwdict[(f,s)]=  avd[f][s]
    return(pwdict.copy(), listofkeys)

def adddyn(listoftuples, pairwisedist, actualdict, allkeys, sequencesdct):
    tmdct = actualdict.copy()
    lenot = [len(t) for t  in listoftuples]
    maxlen = max(lenot)
    listoftuples = [t for t in listoftuples if len(t) == maxlen] 
    for t in listoftuples:#start with pairwise dist (i,j)
         #print("TUPLE ",t)
         specdone=[sequencesdct[e].species for e in t]
         for k in allkeys: #should be single keys
            #print("ALLKEYS k",k)
            if sequencesdct[k].species not in specdone:
                #if k not in t: #and species not in etc
                newtup=tuple(chain.from_iterable([t,[k]]))
                if newtup not in tmdct:
                    #print("NEWTUP",newtup)
                    specdone.append(sequencesdct[k].species)
                    addscore=0
                #print(newtup,t)
                    for l in t:
                        #print("WILL BE added to ",l)
                        try:
                            addscore = addscore+pairwisedist[(l,k)]
                            #annoy("added:",(l,k), pairwisedist[(l,k)]) 
                        except KeyError as ke:
                            try:
                                addscore = addscore+pairwisedist[(k,l)]#both di i,j an j,i
                                annoy("k, added:",(k,l), pairwisedist[(k,l)])
                            except KeyError:
                                pass
                                print("KEYERROR",ke)
                    tmdct[newtup] = tmdct[t]+addscore
    return(tmdct.copy())

def algo_globsum(avd, seqDct, defaults):
    processed, unprocessed, coll, uncoll, species2id  = initCollProc(avd, seqDct)
    pairwise,allkeys = getPairwiseAsTuples(avd)
    #annoy("#@\t",pairwise, allkeys)
    actual = pairwise.copy()
    lot = actual.keys()
    keylengths = [len(key) for key in lot]
    while(max(keylengths) < len(unprocessed)): #proxy for num of spec
        newdct = adddyn(lot, pairwise, actual, allkeys, seqDct)
        actual = newdct 
        lot = newdct.keys()
        keylengths = [len(key) for key in lot]
    tupkeys = []
    scores = []
    for k,v in actual.items():
        tupkeys.append(k)
        scores.append(v)
    tmax=max(scores)
    tieList = [(tupkeys[n],e,n) for (n, e) in enumerate(scores) if e == tmax]
    #prefer defaults
    tiedef = {}
    for tie in tieList:
        tiedef[tie[0]]=sum([1 for x in tie[0] if x in defaults])
    maxDefaults = max(tiedef.values())
    proc = [k for k, v in tiedef.items() if v == maxDefaults]
    first = proc[0]
    #proc = [f for f in first[0]]
    return(seqDct, first, species2id)
def findGlobalMax(avd, seqDct, singles, all, coll, uncoll, defaultForms):
    """Find the best scoring pair. In case of a tie, prefer default models. """
    globMax = -1
    globMaxIds = ""
    globMaxSps = ""
    #print(avd)
    for u in uncoll: #uncollected species
        tmpspec = seqDct[u].species
        
        if avd[u]:
            skl=list(avd[u].keys())
            scorel=list(avd[u].values())
            tmpMax = max(scorel)
            if tmpMax > globMax:
                globMax = tmpMax
                globMaxIds = (skl[scorel.index(tmpMax)], u)
                tieList = [n for (n, e) in enumerate(scorel) if e == tmpMax]
                #print(u,"tmpmax: ", [(n,e) for (n, e) in enumerate(scorel) if e == tmpMax],tieList)
                if len(tieList) >1:
                    annoy("tmpMax ties", [(n,e,skl[n],u ) for (n, e) in enumerate(scorel) if e == tmpMax])
                #favoring defaults:
                fav = [(skl[n],u)  for n in tieList if (skl[n] in defaultForms and u in defaultForms)]
                fav1 = [(skl[n],u)  for n in tieList if (skl[n] in defaultForms or u in defaultForms)]
                if fav:
                    #print("fav ", fav)
                    globMaxSps = (seqDct[fav[0][0]].species, seqDct[fav[0][1]].species)
                    globMaxIds = (fav[0][0],fav[0][1])
                elif fav1:
                    #print("fav1", fav1)
                    fav1tmp = fav1[0] #tuple
                    globMaxSps = (seqDct[fav1tmp[0]].species, seqDct[fav1tmp[1]].species)
                   #################what about globmaxids??
                    globMaxIds =  (fav1tmp[0], fav1tmp[1])
                else:
                    globMaxSps = (seqDct[globMaxIds[0]].species, seqDct[u].species)
                    globMaxIds = (globMaxIds[0], u)

            else:
                pass
    return(globMaxIds, globMaxSps)

def algo(avd, seqDct, singles, all, defaultForms):
    print("AVD",avd)
    """Process distance matrix"""
    # avd is a two dimensional dictionary with sequence identifiers as indices
    # it represents the distance / scoring matrix
    # the seqDct maps identifiers to the actual ScytheSequence objects 
    # so that checking for their species etc becomes more convenient
    processed, unprocessed, coll, uncoll, species2id  = initCollProc(avd, seqDct)#, singles, all) 
    if not GLOBMAX:
        processed, unprocessed, coll, uncoll, species2id = checkSingle(avd, seqDct, singles, all)  
    #singles round
        annoy("# After singles check:","\n#\t unprocessed",unprocessed,
              "\n#\t processed",processed,"\n#\t collected",coll,"\n#\t uncollected",
              uncoll,"\n# ", species2id)
    
    if not coll or GLOBMAX:
        globMaxIds, globMaxSps = findGlobalMax(avd, seqDct, singles, all, coll, uncoll, defaultForms)
        annoy("# global max pair:", globMaxIds, globMaxSps)
        print("GLM",globMaxIds)
        coll.add(globMaxIds[0])
        coll.add(globMaxIds[1])
        processed.add(globMaxSps[0])
        processed.add(globMaxSps[1])
        uncoll.remove(globMaxIds[0])
        uncoll.remove(globMaxIds[1])
        unprocessed.remove(globMaxSps[0])
        unprocessed.remove(globMaxSps[1])
   
    while(unprocessed):
        if coll: #already collected
            print(avd)
            for c in coll:
                cmax = -1
                cmaxid = ""
                cmaxsp = ""
                for up in unprocessed:#not yet processed
                    for uc in species2id[up]:
                        print(uc)
                        if uc in avd: # is in distance dict
                            try: 
                                tmp = int(avd[uc][c])
                            except TypeError as e:
                                tmp = -1
                            #print("TMP ",tmp, "uc ",uc,"c", c  )
                            if (tmp >cmax or (tmp == cmax and uc in defaultForms)) :
                                # tie resolved
                                cmax = int(avd[uc][c])
                                print("new max ", avd[uc][c], uc, c)
                                cmaxid = uc
                                cmaxsp = up
                                cmax=tmp
                        if c in avd:
                            try:
                                tmp = int(avd[c][uc])
                            except TypeError as e:
                                    tmp = -1
                            if (tmp >cmax or (tmp == cmax and uc in defaultForms)):
                                cmax = int(avd[c][uc])
                                cmaxid = uc
                                cmaxsp = up
                                cmax=tmp
            uncoll.remove(cmaxid)
            unprocessed.remove(cmaxsp)
            coll.add(cmaxid)
            processed.add(cmaxsp)
        else:
            print("?")
    annoy("# to do: ",unprocessed)
    annoy("# processed: ",processed)
    annoy("# trash: ",uncoll)
    annoy("# collected: ",coll)
    annoy("# \t#------------------------#")
    return(seqDct, coll, species2id)

def makeFasta(listofspecies, group, frame, stopAfter, gapOpen, gapExtend):
    singles = set()
    allSpec = set()
    seqDct = {}
    pattern  = re.compile(r"""(.*)\s+([a-zA-Z0-9_.]*)\s+[a-zA-Z0-9_.]*\s+\((.*)\)""")
    outfile = ""
    singleGrp = {}
    #!ToDo check types etc
    sp = {}
    defaultForms = set()
    for l in listofspecies:
        sp[l.name] = l
        for i in l._defForm:
            defaultForms.add(i)
    debug_stop = 0
    for g in group.groups:
        if stopAfter and g > stopAfter:
            break
        spl = list(group.groups[g])
        allSpec = set(spl)
        
        for s in spl:
            #print("FAT, "+frame._fat)
            outfile = frame._fat+".".join([str(g),s,"fa"])
            out = open(outfile, "w")
            spa = group.groups[g][s]
            for locus in spa:
                if len(spa)==1:
                    singles.add(locus)
                try:
                    out.write(sp[s].sequences[locus].toFasta())
                    seqDct[sp[s].sequences[locus].name]=sp[s].sequences[locus]
                except KeyError as ke:
                    print ("Are all gene models in your fasta files? - KeyError for ",ke)
            out.close()
    test = "" 
    deb=0
    for g in group.groups:
        avd = None
        avd = AutoViviDict()
        if stopAfter and g> stopAfter:
            break
        annoy("# \t#-- Processing group",g,"--#")
        spl = list(group.groups[g])
        for i in range(0,len(spl)-1):
            for j in range(i+1,len(spl)):
                outfile = frame._fat+".".join([str(g),spl[i],"fa"])
                fileA = outfile
                outfile = frame._fat+".".join([str(g),spl[j],"fa"])
                fileB = outfile
                outfile=frame._sr+".".join([str(g),spl[i],spl[i+1],"needle"])
                try: 
                    task = frame.callNeedleAll(fileA, fileB, outfile = outfile,stdout=True, gapOpen=gapOpen, gapExtend=gapExtend)
                    fulldata = task.stdout.read()
                    assert task.wait() == 0
                except AssertionError as ae:
                    print(ae)
                    data="#"
                data  =  fulldata.decode("utf-8")
                for l in data.split("\n"):
                    if l.startswith("#"):
                        break
                    else:
                        tmp = pattern.findall(l)
                        if tmp:
                            print(tmp)
                            res = tmp[0]
                            score = int(float(res[2])*10)
                            avd[res[0]][res[1]]=score
                            avd[res[1]][res[0]]=score
        if GLOBSUM:
            yield(algo_globsum(avd, seqDct, defaultForms),str(g))       

        else:
            yield(algo(avd, seqDct, singles, allSpec, defaultForms),str(g))

def runScythe(groups, delim, asID, namesList, cleanUp, stopAfter, faFileList, inDir, outDir, gapOpen, gapExtend, locDir=None, faDir=None):
    stopAfter=int(stopAfter)
    specsList = []
    grpMapList = []
    if locDir:
        locFileList = os.listdir(locDir)
        print("locdir set", locFileList)
    else:
        locFileList = os.listdir(inDir+"/loc/")
    if groups:
        if locDir:
            locfl = [locDir+x for x in locFileList]
            print(locfl,"locdl")
        else:
            locfl = [inDir+"/loc/"+x for x in locFileList]
        dct = GrpParser().groupDct(groups, locf=locfl)
        #print(dct)
    else:
        usage()
   # locFileList = os.listdir(inDir+"/loc/")
    for n,f in zip(namesList,faFileList):
        print(namesList, faFileList,"match")
        """Find matching .loc and .fa files"""
        #locFileList = os.listdir(inDir+"/loc/")
        if locDir:
            locFileList = os.listdir(locDir)
        else:
            locFileList = os.listdir(inDir+"/loc/")
        pf = ".".join(f.split(".")[:-1])
        print("PF",pf)
        pf = pf.split("_")[0]
        print("PF2",pf)
        print("LOC",locFileList)
        locFileList = [x for x in locFileList if x.startswith(pf)]
        print("lfl",locFileList)
        #try:
        #    locFileList[0] == ""
        #except IndexError as ke:
        #    print("Something is wrong. File names don't match, please try renaming.",ke)
        #    exit(2)
        n = n.strip()
        if locDir:
            specsList.append(ScytheSpec(name = n, format = "loc", source = locDir+locFileList[0], fasta = faDir+f))

        else:
            specsList.append(ScytheSpec(name = n, format = "loc", source = inDir+"/loc/"+locFileList[0], fasta = inDir+"/fa/"+f))
        for a in specsList:
            annoy("#\t",a.name,a.fasta)
        if locDir:
            grpMapList.append(ScytheGroupMap(name=n, locfile=locDir+locFileList[0], dct=dct, separator=delim, asID=asID))
        else:
            grpMapList.append(ScytheGroupMap(name=n, locfile =inDir+"/loc/"+locFileList[0], dct=dct, separator=delim, asID=asID))
    for g in grpMapList:
        g.free()
    for sp in specsList:
        sp.fillDefForm()
        sp.fillLociCDS()
        sp.fillSequences(sep = delim, asID=asID)
    
    grp = ScytheGroup("tmpgrp", grpMapList)
    frame = ScytheFrame(path=outDir)
    frame.mkAllDirs()
   
    outfiles = {}
    outfilesGroups = {}
    cnt = 0
    for s in specsList:
        outfile = frame._srfa+".".join([s.name,"fa"])
        outfiles[s.name] = open(outfile, 'a')
    for r,R in makeFasta(specsList, grp, frame, stopAfter, gapOpen,gapExtend):
        outfileGroup = frame._srofa+".".join([R,"fa"])
        outfilesGroups[R] = open(outfileGroup, 'a')
        for s in specsList:   
            tmp = r[1]
            ok  = [x for x in tmp if x in s.cds]
            if ok:
                ok = ok[0]
                outfiles[s.name].write(r[0][ok].toFasta())
                outfilesGroups[R].write(r[0][ok].toFasta())
        cnt+=1
    if cleanUp:
        annoy("# Cleaning up...")
        frame.cleanUp()
        annoy("# ...done.")
        
    for out in outfiles.values():
        out.close()
    print("\ndone.")
    ##########################################################
def main():
    global VERBOSE
    global GLOBMAX
    global GLOBSUM 
    VERBOSE = None
    GLOBMAX = False
    GLOBSUM = False 
    cleanUp = False
    groups = None
    namesList = None
    inDir = None
    outDir = "./"
    fastaList = None
    gffList = None
    delim = None
    asID = 0
    stopAfter = False
    gapOpen=str(10)
    gapExtend =str(0.5)
    isUsingConfig = False
    ##################################
    
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "C:i:o:s:d:a:g:O:E:mcvhM", 
                                        ["config=","in_dir=","out_dir=",
                                         "stop_after=",
                                         "delim=", "asID=",
                                         "groups=","gap_open=", "gap_extend=","global_max","cleanup", 
                                         "verbose", "help","global_sum"]
                                       )
    except getopt.GetoptError as err:
        print (str(err))
        print(logo)
        usage()
    for o, a in opts:
        if o in ("-C", "--config"):
            print(logo)
            parseConfig(a)
            
            isUsingConfig=True
        elif o in ("-i", "--in_dir"):
            inDir = a
        elif o in ("-o", "--outdir"):
            outDir = a
        elif o in ("-s", "--stop_after"):
            stopAfter = int(a)
        elif o in ("-d", "--delim"):
            delim = a
        elif o in ("-a", "--asID"):
            asID = int(a)
        elif o in ("-g", "--groups"):
            groups = a
        elif o in ("-v", "--verbose"):
            VERBOSE = True
        elif o in ("-c", "--cleanup"):
            cleanUp = True
        elif o in ("-m", "--global_max"):
            GLOBMAX = True 
        elif o in ("-M", "--global_sum"):
            GLOBSUM = True
        elif o in ("-O", "--gap_open"):
            gapOpen = a
        elif o in ("-E", "--gap_extend"):
            gapExtend = a
         
        elif o in ("-h", "--help"):
            usage()
        else:
            assert False, "unhandled option"
    if not isUsingConfig:
        
        if VERBOSE:
            for o, a in opts:
                if o in ("-v", "--verbose","-c", "--cleanup"):
                    annoy(o, "set")
                else:
                    annoy(o, "set to",a)
    
        if not (inDir and groups):
            print(logo)
            usage()
        
        try:
            annoy(os.listdir(inDir+"/fa/"))
            annoy(os.listdir(inDir+"/loc/"))
        except OSError as e:
            print(e)
            print("Please provide a directory containing folders 'fa' with fasta files and 'loc' with .loc files.")
            usage()
            
        faFileList = os.listdir(inDir+"/fa/")
        namesList = os.listdir(inDir+"/fa/")
        namesList = [n[0:3] for n in namesList]
        #print (namesList)        
    
        if (len(faFileList)!=len(namesList)):
            print("Names:", namesList)
            print("Fasta files:", faFileList)
            usage()
        
        runScythe(groups=groups, delim=delim, 
                  asID=asID, faFileList=faFileList, 
                  namesList=namesList, cleanUp=cleanUp, 
                  stopAfter=stopAfter, inDir=inDir, outDir=outDir, gapOpen=gapOpen, gapExtend=gapExtend)

#----------------------------------------------------------------#
   
if __name__ == "__main__":
    main()
