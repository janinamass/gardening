#!/usr/bin/env python

import imp
import sys
import subprocess
import os
import glob
import getopt
import re
from itertools import chain
import configparser
import datetime

from helpers.grphelper import GrpParser
from helpers.scythecore import ScytheSpec
from helpers.scythecore import ScytheGroupMap, ScytheGroup
from helpers.scythecore import ScytheFrame
#!to get rid of that:
from helpers.scythecore import AutoViviDict

from helpers.fastahelper import FastaParser
from algo.algomod import AlgoHandler
from algo.algomod import EmptyScoringDctException,EmptySequenceDctException

import threading
import time
import queue
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
    print ("""
    ######################################
    # scythe.py v0.1                     #
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
    -i, --in_dir=DIR                 folder w/ subfolders "fa" and "loc"

    -o, --out_dir=DIR                output directory [default:./]
    -v, --verbose                    be wordy

 algorithm options:
    -R, --sl_ref                     find best matches to reference
    -G, --sl_glob                    best scoring pair as seed
    -M, --mx_sum                     optimize sum of pairwise scores


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
    Please see documentation in the 'docs'-directory.
    """)
    sys.exit(2)
#-----/usage------#

#----parse configuration----#
def parseConfig(pathconfig):
    VERBOSE =False
    CF_MODE = "Mode"
    CF_MODE_use_ensembl = "use_ensembl"
    CF_MODE_use_local_files = "use_local_files"
    CF_PATHS = "Paths"
    CF_PATHS_fasta_directory = "fasta_directory"
    CF_PATHS_loc_directory = "loc_directory"
    CF_PATHS_grp_file = "grp_file"
    CF_PATHS_output_directory = "output_directory"
    CF_CLEANUP = "Cleanup"
    CF_CLEANUP_clean_up_directories = "clean_up_directories"
    CF_RUN="Run_options"
    CF_RUN_max_threads ="max_threads"
    CF_RUN_split_input="split_input"
    CF_PENALTIES = "Penalties"
    CF_PENALTIES_gap_open_cost = "gap_open_cost"
    CF_PENALTIES_gap_extend_cost="gap_extend_cost"
    CF_PENALTIES_substitution_matrix="substitution_matrix"
    CF_ALGORITHM = "Algorithm"
    CF_ALGORITHM_use_global_max ="use_sl_glob"
    CF_ALGORITHM_use_default="use_sl_ref"
    CF_ALGORITHM_use_global_sum="use_mx_sum"
    CF_FASTAHEADER="Fasta_header"
    CF_FASTAHEADER_delimiter = "fasta_header_delimiter"
    CF_FASTAHEADER_part = "fasta_header_part"

    config = configparser.ConfigParser()
    config.read(pathconfig)
    ### GLOBMAX = sl_max GLOBSUM = mx_sum
    global GLOBMAX
    global GLOBSUM

    # check for multiple algos in config
    if config.get(CF_ALGORITHM,CF_ALGORITHM_use_global_max)!="yes":
        GLOBMAX = False
    else:
        GLOBMAX = True
    if config.get(CF_ALGORITHM,CF_ALGORITHM_use_global_sum)!="yes":
        GLOBSUM = False
    else:
        GLOBSUM = True
    if config.get(CF_ALGORITHM,CF_ALGORITHM_use_default)!="yes":
        SL_REF = False
    else:
        SL_REF = True
    if sum([GLOBMAX, GLOBSUM, SL_REF])>1:
        sys.stderr.write("Problem with you config file. Please select one algorithm.\n")
        sys.exit(1)
    elif sum([GLOBMAX, GLOBSUM, SL_REF])<1:
        sys.stderr.write('Problem with you config file. Please select one algorithm ( ...  = "yes").\n')
        sys.exit(1)
    else:
        pass
########################################
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
    delim = config.get(CF_FASTAHEADER,CF_FASTAHEADER_delimiter)
    asID = int(config.get(CF_FASTAHEADER,CF_FASTAHEADER_part))
    stopAfter = False
    gapOpen= config.get(CF_PENALTIES,CF_PENALTIES_gap_open_cost)
    gapExtend =config.get(CF_PENALTIES,CF_PENALTIES_gap_extend_cost)
    faFileList = os.listdir(faDir)
    namesList = os.listdir(faDir)
    namesList = [n[0:3] for n in namesList]
#TODO fix numThreads
    runScythe(groups=groups, delim=delim.strip('"'),
              asID=asID, faFileList=faFileList,
              namesList=namesList, cleanUp=cleanUp,
              stopAfter=stopAfter, inDir=inDir, outDir=outDir,
              gapOpen=gapOpen, gapExtend=gapExtend,
              locDir=locDir,faDir=faDir, numThreads = None)

#----/parse configuration----#

def getPairwiseAsTuples(avd):
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
                        except KeyError as ke:
                            try:
                                addscore = addscore+pairwisedist[(k,l)]#both di i,j an j,i
                            except KeyError:
                                #FixThis
                                pass
                                print("KEYERROR",ke)
                    tmdct[newtup] = tmdct[t]+addscore
    return(tmdct.copy())

#def algo_globsum(avd, seqDct, defaults):
#    processed, unprocessed, coll, uncoll, species2id  = initCollProc(avd, seqDct)
#    pairwise,allkeys = getPairwiseAsTuples(avd)
#    actual = pairwise.copy()
#    lot = actual.keys()
#    keylengths = [len(key) for key in lot]
#    while(max(keylengths) < len(unprocessed)): #proxy for num of spec
#        newdct = adddyn(lot, pairwise, actual, allkeys, seqDct)
#        actual = newdct
#        lot = newdct.keys()
#        keylengths = [len(key) for key in lot]
#    tupkeys = []
#    scores = []
#    for k,v in actual.items():
#        tupkeys.append(k)
#        scores.append(v)
#    tmax=max(scores)
#    tieList = [(tupkeys[n],e,n) for (n, e) in enumerate(scores) if e == tmax]
#    #prefer defaults
#    tiedef = {}
#    for tie in tieList:
#        tiedef[tie[0]]=sum([1 for x in tie[0] if x in defaults])
#    maxDefaults = max(tiedef.values())
#    proc = [k for k, v in tiedef.items() if v == maxDefaults]
#    first = proc[0]
#    #proc = [f for f in first[0]]
#    return(seqDct, first, species2id)
#def findGlobalMax(avd, seqDct, singles, all, coll, uncoll, defaultForms):
#    """Find the best scoring pair. In case of a tie, prefer default models. """
#    globMax = -1
#    globMaxIds = ""
#    globMaxSps = ""
#    if avd =={}:
#        #print("avd empty")
#        return(None,None)
#
#    for u in uncoll: #uncollected species
#        tmpspec = seqDct[u].species
#
#        if avd[u]:
#            skl=list(avd[u].keys())
#            scorel=list(avd[u].values())
#            tmpMax = max(scorel)
#            if tmpMax > globMax:
#                globMax = tmpMax
#                globMaxIds = (skl[scorel.index(tmpMax)], u)
#                tieList = [n for (n, e) in enumerate(scorel) if e == tmpMax]
#                #print(u,"tmpmax: ", [(n,e) for (n, e) in enumerate(scorel) if e == tmpMax],tieList)
#                if len(tieList) >1:
#                    pass
#                    #print("tmpMax ties", [(n,e,skl[n],u ) for (n, e) in enumerate(scorel) if e == tmpMax])
#                #favoring defaults:
#                fav = [(skl[n],u)  for n in tieList if (skl[n] in defaultForms and u in defaultForms)]
#                fav1 = [(skl[n],u)  for n in tieList if (skl[n] in defaultForms or u in defaultForms)]
#                if fav:
#                    #print("fav ", fav)
#                    globMaxSps = (seqDct[fav[0][0]].species, seqDct[fav[0][1]].species)
#                    globMaxIds = (fav[0][0],fav[0][1])
#                elif fav1:
#                    #print("fav1", fav1)
#                    fav1tmp = fav1[0] #tuple
#                    globMaxSps = (seqDct[fav1tmp[0]].species, seqDct[fav1tmp[1]].species)
#                   #################what about globmaxids??
#                    globMaxIds =  (fav1tmp[0], fav1tmp[1])
#                else:
#                    globMaxSps = (seqDct[globMaxIds[0]].species, seqDct[u].species)
#                    globMaxIds = (globMaxIds[0], u)
#
#            else:
#                pass
#    return(globMaxIds, globMaxSps)
#
#def algo(avd, seqDct, singles, all, defaultForms):
#    #print("AVD",avd)
#    if avd =={}:
#        return(None)
#    """Process distance matrix"""
#    # avd is a two dimensional dictionary with sequence identifiers as indices
#    # it represents the distance / scoring matrix
#    # the seqDct maps identifiers to the actual ScytheSequence objects
#    # so that checking for their species etc becomes more convenient
#    processed, unprocessed, coll, uncoll, species2id  = initCollProc(avd, seqDct)#, singles, all)
#    if not GLOBMAX:
#        processed, unprocessed, coll, uncoll, species2id = checkSingle(avd, seqDct, singles, all)
#    #singles round
#        #frame.writeLog("debug","# After singles check:","\n#\t unprocessed "+unprocessed+
#              #"\n#\t processed "+processed+"\n#\t collected "+coll+"\n#\t uncollected "+uncoll+"\n# "+species2id)
#
#    if not coll or GLOBMAX:
#        globMaxIds, globMaxSps = findGlobalMax(avd, seqDct, singles, all, coll, uncoll, defaultForms)
#        #print("# global max pair:", globMaxIds, globMaxSps)
#        #print("GLM",globMaxIds)
#        coll.add(globMaxIds[0])
#        coll.add(globMaxIds[1])
#        processed.add(globMaxSps[0])
#        processed.add(globMaxSps[1])
#        uncoll.remove(globMaxIds[0])
#        uncoll.remove(globMaxIds[1])
#        unprocessed.remove(globMaxSps[0])
#        unprocessed.remove(globMaxSps[1])
#
#    while(unprocessed):
#        if coll: #already collected
#            #print(avd)
#            for c in coll:
#                cmax = -1
#                cmaxid = ""
#                cmaxsp = ""
#                for up in unprocessed:#not yet processed
#                    for uc in species2id[up]:
#                        #print(uc)
#                        if uc in avd: # is in distance dict
#                            try:
#                                tmp = int(avd[uc][c])
#                            except TypeError as e:
#                                tmp = -1
#                            #print("TMP ",tmp, "uc ",uc,"c", c  )
#                            if (tmp >cmax or (tmp == cmax and uc in defaultForms)) :
#                                # tie resolved
#                                cmax = int(avd[uc][c])
#                                #frame.writeLog("debug","new max "+avd[uc][c]+" "+uc+" "+c)
#                                cmaxid = uc
#                                cmaxsp = up
#                                cmax=tmp
#                        if c in avd:
#                            try:
#                                tmp = int(avd[c][uc])
#                            except TypeError as e:
#                                    tmp = -1
#                            if (tmp >cmax or (tmp == cmax and uc in defaultForms)):
#                                cmax = int(avd[c][uc])
#                                cmaxid = uc
#                                cmaxsp = up
#                                cmax=tmp
#            uncoll.remove(cmaxid)
#            unprocessed.remove(cmaxsp)
#            coll.add(cmaxid)
#            processed.add(cmaxsp)
#        else:
#            print("?")
#    #frame.writeLog("debug","#to do: "+unprocessed+"\n# processed: "+processed\
#    #            +"\n# trash: "+uncoll+"\n# collected: "+coll)
#    return(seqDct, coll, species2id)

def makeFasta(listofspecies, group, frame, stopAfter, gapOpen, gapExtend,task,  startAt = None, queue = None):
    print("DEBUg makeFasta")
    singles = {}
    skip = {}
    allSpec = set()
    pattern  = re.compile(r"""(.*)\s+([a-zA-Z0-9_.]*)\s+[a-zA-Z0-9_.]*\s+\((.*)\)""")
    outfile = None
    sp = {}
    for l in listofspecies:
        sp[l.name] = l
    for g in group.groups:
        seqDct = {}
        if stopAfter and g > stopAfter:
            break
        if startAt and g < startAt:
            print("continue", g)
            continue
        spl = list(group.groups[g])
        allSpec = set(spl)
        singles[g] = set()
#collect single model species for this group so it could be skipped if all species are in this set
        for s in spl:
#tmp files
            outfile = frame._fat+".".join([str(g),s,"fa"])
            out = open(outfile, "w")
            spa = group.groups[g][s]
            if len(spa)==1:
                singles[g].add(s)
            for locus in spa:
                    try:
                        out.write(sp[s].sequences[locus].toFasta())
                        seqDct[sp[s].sequences[locus].name]=sp[s].sequences[locus]
                    except KeyError as ke:
                        print ("Are all gene models in your fasta files? - KeyError for ",ke)
            out.close()
        if len(singles[g]) == len(spl):
            print("SKIP ",g)
            skip[g] = True
        else:
            skip[g] = False
        #SKIP
        if skip[g]:
            frame.writeLog("debug","#-- Skipping group "+str(g)+"--#")
            #yield((seqDct,set([x.name for x in seqDct.values()]),"SKIP"),str(g))
            queue.put(((seqDct,set([x.name for x in seqDct.values()]),"SKIP"),str(g)))
        else:
            avd = None
            avd = AutoViviDict()
            if stopAfter and g > stopAfter:
                break
            frame.writeLog("debug","#-- Processing group "+str(g)+"--#")
            spl = list(group.groups[g])
            ah = AlgoHandler()
            for i in range(0,len(spl)-1):
                for j in range(i+1,len(spl)):
                    outfile = frame._fat+".".join([str(g),spl[i],"fa"])
                    fileA = outfile
                    outfile = frame._fat+".".join([str(g),spl[j],"fa"])
                    fileB = outfile
                    outfile=frame._sr+".".join([str(g),spl[i],spl[i+1],"needle"])
                    try:
                        print("CALLING NEEDLE ",g)
                        task = frame.callNeedleAll(fileA, fileB, outfile = outfile,stdout=True, gapOpen=gapOpen, gapExtend=gapExtend)
                        fulldata = task.stdout.read()
                        #print(fulldata)
                        assert task.wait() == 0
                        task.stdout.close()
                    except AssertionError as ae:
                        sys.stderr.write(ae)
                        data="#"
                        sys.stderr.write("WARNING:", fileA, fileB,"excluded")
                        frame.writeLog("error","WARNING:"+fileA+" "+fileB+" excluded, AssertionError")
                        #yield((None,None))
                        queue.put((None,None))
                    data  =  fulldata.decode("utf-8")
                    for l in data.split("\n"):
                        if l.startswith("#"):
                            break
                        else:
                            tmp = pattern.findall(l)
                            if tmp:
                                res = tmp[0]
                                score = int(float(res[2])*10)
                                avd[res[0]][res[1]]=score
                                avd[res[1]][res[0]]=score
            if GLOBSUM:
                #yield(algo_globsum(avd, seqDct, defaultForms),str(g))
                queue.put((algo_globsum(avd, seqDct, defaultForms),str(g)))
            else:
                #yield(ah.sl_ref(scoringDct = avd, sequenceDct = seqDct), str(g))
                queue.put((ah.sl_ref(scoringDct = avd, sequenceDct = seqDct), str(g)))
class MakeFastaThread(threading.Thread):
    def __init__(self,listofspecies, group, frame, stopAfter, gapOpen, gapExtend,task,startAt, queue):
        super(MakeFastaThread, self).__init__()
        self.listofspecies = listofspecies
        self.group = group
        self.frame = frame
        self.stopAfter = stopAfter
        self.gapOpen = gapOpen
        self.gapExtend = gapExtend
        self.task = task
        self.startAt = startAt
        self.queue = queue
    def run(self):
        print("DEBUG, RUUN")
        SEMAPHORE.acquire()
        try:
            makeFasta(listofspecies = self.listofspecies, group = self.group, frame = self.frame,
                    stopAfter = self.stopAfter,  gapOpen = self.gapOpen,  gapExtend = self.gapExtend,
                    task = self.task, startAt =self.startAt, queue = self.queue)
        except Exception as e:
            sys.stderr.write(e)
            sys.exit(1)
        print("semaphore release")
        SEMAPHORE.release()

def runScythe(groups, delim, asID, namesList, cleanUp, stopAfter, faFileList, inDir, outDir, gapOpen, gapExtend, locDir=None, faDir=None, numThreads=None, startAt =0):
    global SEMAPHORE
    print("NUMTHREADS", numThreads)
    SEMAPHORE=threading.BoundedSemaphore(numThreads)
    print(delim, asID, locDir, faDir,inDir)
    stopAfter=int(stopAfter)
    specsList = []
    grpMapList = []
    if locDir:
        if not locDir.endswith(os.sep):
            locDir=locDir+os.sep
    if faDir:
        if not faDir.endswith(os.sep):
            faDir=faDir+os.sep
    if locDir:
        locFileList = os.listdir(locDir)
    else:
        locFileList = os.listdir(inDir+os.sep+"loc")
    if groups:
        if locDir:
            locfl = [locDir+x for x in locFileList]
        else:
            locfl = [inDir+os.sep+"loc"+os.sep+x for x in locFileList]
        dct = GrpParser().groupDct(groups, locf=locfl)

    else:
        usage()
    for n,f in zip(namesList,faFileList):
        print(namesList, faFileList,"match")
        """Find matching .loc and .fa files"""
        if locDir:
            locFileList = os.listdir(locDir)
        else:
            locFileList = os.listdir(inDir+os.sep+"loc"+os.sep)
        locFileListtmp = locFileList
        pf = ".".join(f.split(".")[:-1])
        pf = pf.split("_")[0]
        locFileList = [x for x in locFileList if x.startswith(pf)]

        if len(locFileList) < len(faFileList):
            """less stringend matching"""
            print("match again", locFileListtmp, locFileList, pf, pf[0:3])
            locFileList = [x for x in locFileListtmp if x.startswith(pf[0:3])]

        n = n.strip()
        if locDir:
            specsList.append(ScytheSpec(name = n, format = "loc", source = locDir+locFileList[0], fasta = faDir+f))

        else:
            specsList.append(ScytheSpec(name = n, format = "loc", source = inDir+os.sep+"loc"+os.sep+locFileList[0], fasta = inDir+os.sep+"fa"+os.sep+f))
        if locDir:
            grpMapList.append(ScytheGroupMap(name=n, locfile=locDir+locFileList[0], dct=dct, separator=delim, asID=asID))
        else:
            grpMapList.append(ScytheGroupMap(name=n, locfile =inDir+os.sep+"loc"+os.sep+locFileList[0], dct=dct, separator=delim, asID=asID))

    for g in grpMapList:
        g.free()

    for sp in specsList:
        sp.fillLociCDS()
        sp.fillSequences(sep = delim, asID=asID)
        sp.fillDefForm()

    grp = ScytheGroup("tmpgrp", grpMapList)
    frame = ScytheFrame(path=outDir)
    frame.mkAllDirs()

    outfiles = {}
    outfilesGroups = {}
###june
    outDctSp = {}#key output filename, value scytheSeq
    outDctOg = {}#same but for orthogroup
###
    cnt = 0
    for s in specsList:
#### species output ####

        outfile = frame._srfa+".".join([s.name,"fa"])
        outfiles[s.name] = open(outfile, 'a')
        outfile = frame._srfa+".".join([s.name,"skipped.fa"])
        outfiles[s.name+".skipped"] = open(outfile, 'a')
######################################### parallel
    if stopAfter:
        maxNumGrp = stopAfter
    else:
        print(grp, len(grp.groups))
        maxNumGrp = len(grp.groups)
    if startAt:
        minNumGrp = startAt
    else:
        minNumGrp = 0
    numGrps = maxNumGrp - minNumGrp
#todo keep track of what has already started
    numPerThread = round((float(numGrps)/float(numThreads))+0.5)
    print("debug", startAt,stopAfter, numGrps, numPerThread, numThreads, "NT" )
    qList = []
    tList = []
    for i in range(0, numThreads):
        qList.append(queue.Queue())
        print(i)
        tList.append(MakeFastaThread(listofspecies = specsList, group = grp, frame = frame, stopAfter=stopAfter,gapOpen =  gapOpen,gapExtend = gapExtend, task="needleall", startAt = startAt+i*numPerThread-1, queue = qList[i]))
    for t in tList:
        t.run()

    for k in qList:
        while True:
            print("true")
            try:
                tmp = k.get(True, 15)
                R = tmp[1]
                r = tmp[0]
                print(r, "r\n")
    #for r,R in MakeFastaThread(listofspecies = specsList, group = grp, frame = frame, stopAfter=stopAfter,gapOpen =  gapOpen,gapExtend = gapExtend, task="needleall", startAt = startAt+i*numPerThread).run():
        #########
                if r is None:
                    sys.stderr.write("Failed to process group {}\n".format(str(R)))
                    continue
                else:
                    if r[2] == "SKIP":
                            outfileGroup = frame._srofa+".".join([R,"skipped","fa"])
                            outfilesGroups[R] = open(outfileGroup, 'a')
                    else:
                        outfileGroup = frame._srofa+".".join([R,"fa"])
                        outfilesGroups[R] = open(outfileGroup, 'a')
                    for s in specsList:
                        tmp = r[1]
                        ok  = [x for x in tmp if x in s.cds]
                        if ok:
                            ok = ok[0]
                            if not r[2] =="SKIP":
                                outfiles[s.name].write(r[0][ok].toFasta())
                            else:
                                outfiles[s.name+".skipped"].write(r[0][ok].toFasta())
                            outfilesGroups[R].write(r[0][ok].toFasta())
                    cnt+=1
                    outfilesGroups[R].close()
            except queue.Empty:
                break
###########################################
    #for t in tList:
    #    t.join()

    if cleanUp:
        frame.cleanUp()

    for out in outfiles.values():
        out.close()
    print("\ndone.")

def main():
    global VERBOSE
    global GLOBMAX
    global GLOBSUM
    global SL_REF
    VERBOSE = None
    GLOBMAX = False
    GLOBSUM = False
    SL_REF = False
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
    numThreads = None
    ##################################

    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "C:i:g:o:s:d:a:O:E:N:RGMcvh",
                                        ["config=",
                                        "in_dir=",
                                        "groups=",
                                        "out_dir=",
                                        "stop_after=",
                                        "delim=", "asID=",
                                        "gap_open=",
                                        "gap_extend=",
                                        "num_cores = ",
                                        "sl_ref",
                                        "sl_glob",
                                        "mx_sum",
                                        "cleanup",
                                        "verbose",
                                        "help"]
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
            if not inDir.endswith(os.sep):
                inDir = inDir+os.sep
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
        elif o in ("-R", "--sl_ref"):
            SL_REF = True
        elif o in ("-G", "--sl_glob"):
            GLOBMAX = True
        elif o in ("-M", "--mx_sum"):
            GLOBSUM = True
        elif o in ("-O", "--gap_open"):
            gapOpen = a
        elif o in ("-E", "--gap_extend"):
            gapExtend = a
        elif o in ("-N", "--num_cores"):
            numThreads = int(a)
        elif o in ("-h", "--help"):
            usage()
        else:
            assert False, "unhandled option"
    if not isUsingConfig:

        if VERBOSE:
            for o, a in opts:
                if o in ("-v", "--verbose","-c", "--cleanup"):
                    print(o, "set")
                else:
                    print(o, "set to",a)

        if not (inDir and groups):
            print(logo)
            usage()

        try:
            print(os.listdir(inDir+"fa"))
            print(os.listdir(inDir+"loc"))
        except OSError as e:
            sys.stderr.write(str(e))
            print("Please provide a directory containing folders 'fa' with fasta files and 'loc' with .loc files.\nAlternatively, use Scythe with gui or configuration file\n")
            usage()

        faFileList = os.listdir(inDir+"fa")
        namesList = os.listdir(inDir+"fa")
        namesList = [n[0:3] for n in namesList]
        locDir = inDir+"loc"+os.sep
        locFileList = os.listdir(locDir)
        faDir = inDir+"fa"+os.sep
        faFileList = os.listdir(faDir)

        print("debug", locDir, locFileList)

        if (len(faFileList)!=len(namesList)) or (len(locFileList)!=len(namesList)):
            sys.stderr.write("Number of files doesn't match. Please check {} and {}\n".format(locDir, faDir))
            usage()


        runScythe(groups=groups, delim=delim,
                  asID=asID, faFileList=faFileList,
                  namesList=namesList, cleanUp=cleanUp,
                  stopAfter=stopAfter, inDir=faDir, outDir=outDir, gapOpen=gapOpen, gapExtend=gapExtend, locDir = locDir, faDir=faDir, numThreads = numThreads)

#----------------------------------------------------------------#

class ThreadedScythe(threading.Thread):
    def __init__(self, queue, argdct):
        threading.Thread.__init__(self)
        self.queue = queue
        self.argdct = argdct

    def run(self):
        try:
            runScythe(**self.argdct)
            self.queue.put(1)
        except Exception as e:
            print(e)
            sys.exit(1)

if __name__ == "__main__":
    main()
