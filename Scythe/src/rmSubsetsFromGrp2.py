#!/usr/bin/python3
import re,sys, getopt

###################
# 01/03/2014 [JM] #
###################

def usage():
    print ("""
    #####################################################
    #  rmSubSetsFromGrp.py -g groups.grp -o new.grp -r  #
    #####################################################
   
    options:
    -g, --grp=FILE.grp     
    -o, --output=OUTFILE.grp output file [default: FILE.allspec.grp]
    [-r, --rename discard old orthogroup ids and start numbering from 0]    
    -h, --help         prints this    
    [-n, --numspec=N    min number of species ] 
    ------------
    .grp format: GroupID\tgeneIDiSp1\tgeneIDjSp2\t...geneIDkSpn
    """)
    sys.exit(2)


class Orthogroup(object):
    def __init__(self,list):
        self._list = list[1:]
        self._isSubset = None
        self._size = len(self._list)
        self._checked =False 
        self._id = list[0]
    def setIsSubset(self,bool):
        self._isSubset=bool
    def compareToOrthogroup(self,other):
        self._checked = True
        other._checked = True
        notInSelf = [l for l in self._list if l not in other._list]
        notInOther = [l for l in other._list if l not in self._list]
        #print("Other unique: ",len(notInSelf), notInSelf)
        #print("Self unique:", len(notInOther), notInOther)
        if len(notInOther) == 0 and len(notInSelf)==0:
            pass
            #print("Same thing.")
        elif len(notInOther)==0 and len(notInSelf)>0:
            self._isSubset=True
        elif len(notInSelf)==0 and len(notInOther)>0:
            other._isSubset = True
        return((notInSelf, notInOther))
    def compareToOrthogroup2(self,other):
        self._checked = True
        other._checked = True
        selfSet = set(self._list)
        otherSet = set(other._list)
        if self._size > other._size:
            bigger=self
            smaller = other
        else:
            bigger = other 
            smaller = self
        intersect = selfSet.intersection(otherSet)
        if len(intersect)>=1: #something is in common
            #bigger is better
            bigger._isSubset=False
            smaller._isSubset=True
            
    def compareToOrthogroup3(self,other):
        #print("cmp", self._id, other._id)
        self._checked = True
        other._checked = True
        selfSet = set(self._list)
        otherSet = set(other._list)
        if self._size > other._size:
            bigger=self
            smaller = other
        else:
            bigger = other 
            smaller = self
       # print("SZ",bigger._size, smaller._size)
        intersect = selfSet.intersection(otherSet)
        if len(intersect)==self._size==other._size:
            #print("same", bigger._id, smaller._id)
            return(bigger)
        elif len(intersect)==smaller._size != bigger._size:
            #print("bigger", bigger._id)
            return(bigger)
        
        elif len(intersect)>=1 and self._id != other._id: #something is in common
            new = self.expand(other)
            #print("merge", new._id)
            other._isSubset=True
            self._isSubset=True
            return(new)
        else:
            return(None)
    def expand(self, other):
        tmpset = set()
        tmplist = []
        for i in self._list:
            tmplist.append(i)
        for j in other._list:
            #if j in tmplist:
                #pass
                #print("wtf?")
            #else:
            tmplist.append(j)
        #print(tmplist)
        tmpset = set(tmplist)
        
        #for l in self._list:
        #    tmplist.append(l)
        #for l in other._list:
        #    tmplist.append(l)
        #print(tmpset)
        newlist=list(tmpset)
        newid= str(self._id)+"_"+str(other._id)
        initlist=[newid]
        for i in newlist:
            initlist.append(i)
        new = Orthogroup(initlist)
        
        #new._checked=True
        #new._isSubset=False
        #other._isSubset = True
        #self._isSubset = True
        new._id = str(self._id)+"_"+str(other._id)
        #self=new
        return(new)
    def toString(self):
        s = "\t".join(self._list)
        return(s)
    
    
def unionOrthogroups(grp, out,numspec, rename=False):
    #print(numspec)
    seen = []
    idDct = dict()
    knownIDs = set()
    ok = 0
    redundant = 0
    isOK = True
    ml =0
    try:
        g = open(grp,'r')
    except IOError as e:
        print(e)
    for l in g:
            l = l.rstrip()
            checkThis = False
            tmp = l.split("\t")
            newOrtho = Orthogroup(tmp)
            #print(tmp,newOrtho._list,newOrtho._size, numspec)
            todelete =None
            for k in tmp[1:]:
                #print(k,tmp)
                if k in knownIDs:
                    checkThis = True
                #print(checkThis)
                knownIDs.add(k)
                if not k in idDct:
                    idDct[k] = [newOrtho]
                else:
                    idDct[k].append(newOrtho)
            if checkThis ==True:
                #print("checkthis", checkThis)
                if newOrtho._size==numspec: #all species
                    print("maxlength")
                    ml+=1
                    newOrtho._isSubset=False
                    newOrtho._checked =True
                    for k in tmp[1:]:
                        for l in idDct[k]:
                            if l._id==newOrtho._id:
                                pass
                            #print(l._id)
                            elif l._checked==False:
                                l._isSubset=True
                                l._checked=True
                                print(l._id,"--")
                            elif l._checked ==True and l._isSubset==False:
                                l._isSubset=True
                                print(l._id)
            
                #else: #smaller
                #    for k in tmp:
                #        tmpids = [j._id for j in idDct[k] ]
                #        #print ("tmpids", tmpids)
                #        for s in idDct[k]:
                #            if not s._checked:
                #            #print(idDct[k])
                #                new = s.compareToOrthogroup3(newOrtho)
                #                #print(new, new._id)
                #                seen.append(new)
                #                if s._isSubset:
                #                    todelete =s
                #                elif newOrtho._isSubset:
                #                    todelete=newOrtho
                #        if todelete:
                #            try:
                #                seen.remove(todelete)
                #            except ValueError as e:
                #                #print("#",e, todelete.toString())
                #                pass
                        #idDct[k].append(new)
            #if newOrtho._isSubset==False:
            #if newOrtho._size==numspec: # has not been seen before but is full size
            #    newOrtho._checked=True
            #    newOrtho._isSubset=False
            
            seen.append(newOrtho)
                #print(newOrtho.toString())
                #print(len(seen))
    #for s in range(0,int(len(seen)+1/2)):
    #    for t in range(int(len(seen)-1/2),len(seen)):
    #        seen[s].compareToOrthogroup(seen[t])
    #add non-redundant ones
    print("allspecies ",ml)
    g.close()
    for s in seen:
        if s._checked ==False:
            if s._isSubset==None:
                #print("ok")
                s._checked=True
                s._isSubset=False
            
            
    warn = 0
    proc = 0
    cnt=0
    redcnt=0
    try:
        red = open(out+".rm",'w')
        out = open(out,"w")
        
    except IOError as e:
        print(e)
        exit(1)
    
    for s in seen:
        if s._checked == False:
            warn+=1
           # print (s.toString())
        if s._checked == True:
            proc +=1
        
            if s._isSubset==False:
                #print(str(s._id)+"\t"+str(cnt)+"\t"+s.toString())
                if rename==True:
                    out.write(str(cnt)+"\t"+s.toString())
                else:
                    out.write(str(s._id)+"\t"+s.toString())
                out.write("\n")
                cnt+=1
            if s._isSubset==True:
                #print(str(s._id)+"\t"+str(redcnt)+"\t"+s.toString())
                red.write(str(s._id)+"\t"+s.toString())
                red.write("\n")
                redcnt+=1
            #print("WARNING", s, "not checked!")
    print("# warnings ", warn)
    print("# processed ", proc)

 
def checkConsistency(grp, outfile):
    knownGenes= set()
    knownOrthogroups = []
    godict = {}
    proc = set()
    try:
        g = open(grp,'r')
    except IOError as e:
        print(e)
    for l in g:
            l = l.rstrip()
            tmp = l.split("\t")
            #print("tmp",tmp)
            newOrtho= Orthogroup(tmp)
            #print("NO",newOrtho._id, newOrtho.toString())
            for t in tmp:                
                if t in knownGenes:
                    tmplist = []
                    for o in godict[t]:
                        #print(o._id)
                        #print("no",newOrtho.toString())
                        new = newOrtho.expand(o)
                        new._isSubset=False
                        #print(new._id, new.toString())
                        o._isSubset = True
                        newOrtho._isSubset = True
                        tmplist.append(new)
                    for l in tmplist:
                        godict[t].append(l)
                #tmplist=[]
                else:
                    #print("NO2",newOrtho._id, newOrtho.toString())
                    knownGenes.add(t)
                    godict[t] = [newOrtho]
    g.close()
    try:
        g = open(grp,'r')
        o = open(outfile,'w')
    except IOError as e:
        print(e)
    for l in g:
        l = l.rstrip()
        tmp = l.split("\t")
        for t in tmp[1:]:
            nm = [t._id for t in godict[t]]
            #print(t, nm)
            nmm = [t for t in godict[t] if  t._isSubset==False]
           # print(nmm)
            out=None
            if nmm !=[]:
                out = nmm[-1]
            #try:
                #for out in nmm:
                if out._id not in proc:
                    #print(out._id+"\t"+out.toString())
                    o.write(out._id+"\t"+out.toString())
                    o.write("\n")
                proc.add(out._id)
            #except IndexError as e:
            #    pass
            #if out ==[]:
            #    pass
            #else:
            #    out = out[0]
               # print(t, nmm[:-1], out._id)
    g.close()
    o.close()
                
def filterRedundancy(grp, out):
    seen = set()
    ok = 0
    redundant = 0
    isOK = True
    try:
        g = open(grp,'r')
        o = open(out,'w')
    except IOError as e:
        print(e)
    for l in g:
            l = l.rstrip()
            isOK = True
            tmp = l.split("\t")
            for t in tmp[1:]:
                #print(t)
                if t in seen:
                    isOK = False
                
                seen.add(t)
            if isOK:
                ok+=1
                o.write(l)
                o.write("\n")
            else:
                #print("RED",tmp)
                redundant+=1
    print("# non redundant ",ok)
    print("# redundant ",redundant)  
    g.close()
    o.close()
def filterGroups(grp, out, numspec=None,rename=False):
    try:
        g = open(grp,'r')
    except IOError as e:
        print(e)
        usage()
    
    if numspec:
        max = numspec
    if not numspec:
        max = 0
        for l in g:
            l = l.strip()
            tmp = len(l.split("\t")[1:])
            #print(tmp)
            if tmp > max:
                max = tmp
        g.close()
    numspec=max
    #print(numspec)
    try:
        g = open(grp,'r')
    except IOError as e:
        print(e)
        usage()
    try:
        #print(out)
        out = open(out,'w')
    except IOError as e:
        print(e)
        usage()
     
    cnt=0
    for l in g:
        ln = l.rstrip()
        tmp = len(ln.split("\t")[1:])
        if tmp == numspec:
            if not rename:
                out.write(l)
            else:
                out.write(str(cnt)+"\t"+"\t".join(l.split("\t")[1:]))
                cnt+=1
    out.close()
    return(numspec)
def main():
    out = None
    outNoSubsets = None
    rename=False
    numspec = None
    ###################################
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "g:n:o:hr", ["grp=","numspec=","output=", "help","rename"])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
    for o, a in opts:
        if o in ("-g", "--grp"):
            grp= a
        elif o in ("-h", "--help"):
            usage()
        elif o in ("-n", "--numspec"):
            numspec = int(a)
        elif o in ("-o", "--output"):
            out = a
            tmp=a    
        elif o in ("-r", "--rename"):
            rename = True
        else:
            assert False, "unhandled option"
    ########################################
    if grp is None:
        usage()
    if out is None:
        tmp="tmp."
        if grp.endswith(".grp"):
            tmp = grp[:-3]
    out = tmp+"allspec.grp"
    outNoSubsets = tmp+"nosubsets.grp"
    ########################################
    print("# Output: ", out)
    print("# Output: ", outNoSubsets)
   
    #out = open(out, "w")
    #out.write(res)
    #filterGroups(grp, out, numspec, rename)
    #filterRedundancy(out,"tmp")
    #filterRedundancy(grp,"tmp")
    #filterRedundancyPlusPlus(grp,outNoSubsets, rename)
    numspec= filterGroups(grp, out, numspec, rename)
    checkConsistency(grp,"tmp")
    filterRedundancy("tmp", outNoSubsets)
    #filterRedundancy(outNoSubsets,"tmp")

    #filterRedundancyPlusPlus(grp,outNoSubsets, rename)
    #unionOrthogroups(grp, outNoSubsets, numspec, rename)
    #filterGroups(outNoSubsets, out, numspec, rename=False)
    #filterRedundancy(outNoSubsets,"tmp")
    #filterRedundancy(out,"tmp")

    
if __name__=="__main__":
    main()