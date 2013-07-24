class GrpParser(object):
    """read .loc files. Return dictionary with locus ID as key and transcript list as value"""
    def readLoc(self, locfiles, speciesNames = None):
        firstLine=True
        res = {}
        for lf in locfiles:
            tmpfile = open(lf,'r')        
            for ln in tmpfile:
                if firstLine:
                    firstLine = False
                ln = ln.rstrip()
                tmp = ln.split("\t")
                k,v = str(tmp[0]),tmp[1:]
                res[k]=v
            tmpfile.close()
        return(res)
    
    def readGroups(self, grpf, locf, speciesNames = None):
            """Iterate over grp return numbered list with orthologous transcripts."""
            grp = open(grpf, 'r')
            locDct = self.readLoc(locfiles =locf)
            traDct ={}
            #locDct has loci as keys, lists of transcripts as values
            for l in grp:
                tmp = l.strip().split("\t")
                grpID = tmp[0]
                tmpLoci=tmp[1:]
                try:                    
                    tr = [locDct[x] for x in tmpLoci]
                except KeyError as ke:
                    print("key error",ke," @readGroups. Locus does not appear to have transcripts.")
                traDct[grpID] = tr
                yield int(grpID), traDct[grpID]
            grp.close()
            
    def groupDct(self, grpf, locf,speciesNames = None):
        #This will need to return the actual CDS identifiers (from 
        # the .loc file
            res = {}
            for sth in self.readGroups(grpf, locf,speciesNames):
                res[sth[0]]=sth[1:][0]
            return res 
    

###







class ProteinOrthoParser(object):
    def readInfo(self, po):
        po = open(po, 'r')
        res = ""
        for l in po:
            if l.startswith("#"):
                res+=l
        po.close()
        return(res)
    
    def readGroups(self, po, speciesNames = None):
        """Iterate over proteinOrtho output and return numbered list with orthogroups."""
        po = open(po, 'r')
        cnt = 0
        #if columnNames aren't provided, use first line as names
        l = po.readline()
        if not speciesNames:
            sn = l.strip().split("\t")[3:]
            #print(sn)
        else:
            if len(speciesNames) == len(l.strip().split("\t")[3:]):
                sn = speciesNames
            else:
                raise ScytheError("Number of speciesNames must equal the number of species in the file.")
        for l in po: 
            if l.startswith("#"):
                pass
            else:
                tmp_list=l.strip().split("\t")[3:]
                tmp_list = [l.split(",") for l in tmp_list]
                tmp_list = [val for tmp in tmp_list for val in tmp]
                yield cnt, tmp_list, sn
                cnt +=1
        po.close()
    
    def groupDct(self, po, speciesNames = None):
        res = {}
        for gr,ids,spec in self.readGroups(po, speciesNames):
            res[gr] = ids
        return res 
    
    def findOrthologs(self, po, geneID):
        #return geneID's friends
        res = []
        for a in self.readGroups(po):
            if geneID in a[1]:
                res.append(a)
        return res
