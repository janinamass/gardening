class AlgoHandler(object):
    def __init__(self):
        pass

    def getSingleGMS(self, sequenceDct):
        """Return list of species that only have one gene model"""
        res = [s.species for s in sequenceDct.values() if s.isSingle]
        return(res)

    def initCollProc(self, scoringDct, sequenceDct):
        """Initialize sets of species. """
        "speciesID -> [gmIDs]"
        species2id = {}

        #species:
        processed = set()
        unprocessed = set()

        #seq:
        coll = set()
        uncoll = set()

        for fkey in scoringDct:
            uncoll.add(fkey)
            unprocessed.add(sequenceDct[fkey].species)
            append_add(species2id, sequenceDct[fkey].species, fkey, unique = True)
            for skey in scoringDct[fkey]:
                uncoll.add(skey)
                unprocessed.add(sequenceDct[skey].species)
                append_add(species2id, sequenceDct[skey].species, skey, unique=True)
        return(processed, unprocessed, coll, uncoll, species2id)

    def getMaxSeed(self, scoringDct, sequenceDct, uncoll):
        """Find the best scoring pair. In case of a tie, prefer default models. """
        globMax = -1
        globMaxIds = None
        globMaxSps = None
        if scoringDct =={}:
            print("scoringDct empty")
            raise Warning("ScoringDct empty")
            return(None,None)
        defaultForms = [s.name for s in sequenceDct.values() if s.isReference]
        for u in uncoll: #uncollected species
            tmpspec = sequenceDct[u].species
            if scoringDct[u]:
                skl=list(scoringDct[u].keys())
                scorel=list(scoringDct[u].values())
                tmpMax = max(scorel)
                if tmpMax > globMax:
                    globMax = tmpMax
                    globMaxIds = (skl[scorel.index(tmpMax)], u)
                    tieList = [n for (n, e) in enumerate(scorel) if e == tmpMax]
                    if len(tieList) >1:
                        pass
                    #favoring defaults:
                    fav = [(skl[n],u)  for n in tieList if (skl[n] in defaultForms and u in defaultForms)]
                    fav1 = [(skl[n],u)  for n in tieList if (skl[n] in defaultForms or u in defaultForms)]
                    if fav:
                        globMaxSps = (sequenceDct[fav[0][0]].species, sequenceDct[fav[0][1]].species)
                        globMaxIds = (fav[0][0],fav[0][1])
                    elif fav1:
                        fav1tmp = fav1[0] #tuple
                        globMaxSps = (sequenceDct[fav1tmp[0]].species, sequenceDct[fav1tmp[1]].species)
                        globMaxIds =  (fav1tmp[0], fav1tmp[1])
                    else:
                        globMaxSps = (sequenceDct[globMaxIds[0]].species, sequenceDct[u].species)
                        globMaxIds = (globMaxIds[0], u)

                else:
                    pass
            else:
                sys.stderr.write(scoringDct, u)
                raise Exception("Sth wrong with ScoringDct")
        return(globMaxIds, globMaxSps)

    def sl_ref(self, scoringDct = {}, sequenceDct = {}, referenceAlgo = True):
        """Reference algorithm"""
        print("debug", "sl_ref")
        if not scoringDct:
            raise EmptyScoringDctException("sth wrong during sl_ref")
        if not sequenceDct:
            raise EmptyScoringDctException("sth wrong during sl_ref")
        singleGMSpec = self.getSingleGMS(sequenceDct)
        processed, unprocessed, coll, uncoll, species2id  = self.initCollProc(scoringDct, sequenceDct)
        if referenceAlgo :
        ### add single GM ###
            for sgms in set(singleGMSpec):
                print("ASSERT",species2id[sgms],sgms)
                assert len(species2id[sgms]) == 1
                processed.add(sgms)
                unprocessed.remove(sgms)
                coll.add(species2id[sgms][0])
                uncoll.remove(species2id[sgms][0])
        #####################

        while(unprocessed):
            if not coll:
                max_seqid,max_specid = self.getMaxSeed(scoringDct = scoringDct, sequenceDct= sequenceDct, uncoll = uncoll)
                for j,k in  zip(max_seqid,max_specid):
                    coll.add(j)
                    uncoll.remove(j)
                    processed.add(k)
                    unprocessed.remove(k)
            for c in coll:
                cmax = -1
                print("cmax",cmax)
                cmaxid = None
                cmaxsp = None
                for up in unprocessed:#spec
                    for uc in species2id[up]:#seq
                        print(uc,"UUUUUUCCCCCC", species2id, up, sequenceDct[uc])
                        if uc in scoringDct: # is in distance dict
                            try:
                                tmp = int(scoringDct[uc][c])
                            except TypeError as e:
                                tmp = -1
                                #raise Warning("tmp -1")
                            if (tmp >cmax or (tmp == cmax and sequenceDct[uc].isReference)) :
                                # tie resolved
                                cmax = int(scoringDct[uc][c])
                                cmaxid = uc
                                cmaxsp = up
                                cmax=tmp
                                print("cmax",cmax)
                        elif c in scoringDct:
                            try:
                                tmp = int(scoringDct[c][uc])
                            except TypeError as e:
                                    tmp = -1
                                 #   raise Warning("TMP -1")
                            if (tmp >cmax or (tmp == cmax and sequenceDct[uc].isReference)):
                                cmax = int(avd[c][uc])
                                cmaxid = uc
                                cmaxsp = up
                                cmax=tmp
                                print("Cmax",cmax)

            uncoll.remove(cmaxid)
            unprocessed.remove(cmaxsp)
            coll.add(cmaxid)
            processed.add(cmaxsp)
        return(sequenceDct, coll, species2id)

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




#    if not coll or GLOBMAX:
#        globMaxIds, globMaxSps = findGlobalMax(avd, seqDct, singles, all, coll, uncoll, defaultForms)
        #print("# global max pair:", globMaxIds, globMaxSps)
        #print("GLM",globMaxIds)
#        coll.add(globMaxIds[0])
#        coll.add(globMaxIds[1])
#        processed.add(globMaxSps[0])
#        processed.add(globMaxSps[1])
#        uncoll.remove(globMaxIds[0])
#        uncoll.remove(globMaxIds[1])
#        unprocessed.remove(globMaxSps[0])
#        unprocessed.remove(globMaxSps[1])







    def sl_glob(scoringDct = {},sequenceDct = {}):
        print("debug", "sl_glob")
        if not scoringDct:
            raise EmptyScoringDctException("sth wrong during sl_glob")
        if not sequenceDct:
            raise EmptyScoringDctException("sth wrong during sl_glob")
        self.sl_ref(scoringDct = scoringDct, sequenceDct = sequenceDct, referenceAlgo=True)

    def mx_sum(scoringDct = {},sequenceDct = {}):
        print("debug", "mx_sum")
        if not scoringDct:
            raise EmptyScoringDctException("sth wrong during mx_sum")
        if not sequenceDct:
            raise EmptyScoringDctException("sth wrong during mx_sum")


class EmptyScoringDctException(Exception):
    pass

class EmptySequenceDctException(Exception):
    pass


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
#
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
    def findGlobalMax(self, scoringDct, sequenceDct,  all, coll, uncoll, defaultForms):
        """Find the best scoring pair. In case of a tie, prefer default models. """
        print("\n\nFGM\n\n")
        globMax = -1
        globMaxIds = ""
        globMaxSps = ""
        if scoringDct =={}:
            raise Exception("No Scores")
            return(None,None)

        for u in uncoll: #uncollected seq
            tmpspec = sequenceDct[u].species
            print(scoringDct,"SCOOOO")
            if scoringDct[u]:
                skl=list(scoringDct[u].keys())
                scorel=list(scoringDct[u].values())
                tmpMax = max(scorel)
                if tmpMax > globMax:
                    globMax = tmpMax
                    globMaxIds = (skl[scorel.index(tmpMax)], u)
                    tieList = [n for (n, e) in enumerate(scorel) if e == tmpMax]
                    #print(u,"tmpmax: ", [(n,e) for (n, e) in enumerate(scorel) if e == tmpMax],tieList)
                    if len(tieList) >1:
                        pass
                        #print("tmpMax ties", [(n,e,skl[n],u ) for (n, e) in enumerate(scorel) if e == tmpMax])
                    #favoring defaults:
                    fav = [(skl[n],u)  for n in tieList if (skl[n] in defaultForms and u in defaultForms)]
                    fav1 = [(skl[n],u)  for n in tieList if (skl[n] in defaultForms or u in defaultForms)]
                    if fav:
                        globMaxSps = (sequenceDct[fav[0][0]].species, sequenceDct[fav[0][1]].species)
                        globMaxIds = (fav[0][0],fav[0][1])
                    elif fav1:
                        #print("fav1", fav1)
                        fav1tmp = fav1[0] #tuple
                        globMaxSps = (sequenceDct[fav1tmp[0]].species, sequenceDct[fav1tmp[1]].species)
                        globMaxIds =  (fav1tmp[0], fav1tmp[1])
                    else:
                        globMaxSps = (sequenceDct[globMaxIds[0]].species, sequenceDct[u].species)
                        globMaxIds = (globMaxIds[0], u)

                else:
                    pass
        return(globMaxIds, globMaxSps)

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
#
#
