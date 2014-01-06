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
        self._isSubset = False
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
        if len(intersect)>1: #something is in common
            #bigger is better
            bigger._isSubset=False
            smaller._isSubset=True
            #self.expand(bigger)
        #print("Test: Is self bigger", bigger==self)
        #print("Test: Is other bigger", bigger==other)

        #notInSelf = [l for l in self._list if l not in other._list]
        #notInOther = [l for l in other._list if l not in self._list]
        #print("Other unique: ",len(notInSelf), notInSelf)
        #print("Self unique:", len(notInOther), notInOther)
        #if len(notInOther) == 0 and len(notInSelf)==0:
        #    pass
            #print("Same thing.")
        #elif len(notInOther)==0 and len(notInSelf)>0:
        #    self._isSubset=True
        #elif len(notInSelf)==0 and len(notInOther)>0:
        #    other._isSubset = True
        #return((notInSelf, notInOther))
    def compareToOrthogroup3(self,other):
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
        if len(intersect)>1: #something is in common
            self.expand(other)
            other._isSubset=True
    def expand(self, other):
        tmplist = []
        for l in self._list:
            tmplist.append(l)
        for l in other._list:
            tmplist.append(l)
        new = Orthogroup(tmplist)
        new._checked=True
        new._id = str(self._id)+"_"+str(other._id)
        self=new
    def toString(self):
        s = "\t".join(self._list)
        return(s)
    
    
def unionOrthogroups(grp, out,numspec, rename=False):
    seen = []
    idDct = dict()
    knownIDs = set()
    ok = 0
    redundant = 0
    isOK = True
    try:
        g = open(grp,'r')
    except IOError as e:
        print(e)
    for l in g:
            l = l.strip()
            checkThis = False
            tmp = l.split("\t")
            newOrtho = Orthogroup(tmp)
            todelete =None
            for k in tmp:
                if k in knownIDs:
                    checkThis = True
                #print(checkThis)
                knownIDs.add(k)
                if not k in idDct:
                    idDct[k] = [newOrtho]
                else:
                    idDct[k].append(newOrtho)
            if checkThis ==True:
                for s in idDct[k]:
                    s.compareToOrthogroup3(newOrtho)
                    if s._isSubset:
                        todelete =s
                    elif newOrtho._isSubset:
                        todelete=newOrtho
                #if todelete:
                #    try:
                #        seen.remove(todelete)
                #    except ValueError as e:
                        #print("#",e, todelete.toString())
                #        pass
            #if newOrtho._isSubset==False:
            if newOrtho._size==numspec:
                newOrtho._checked=True
                newOrtho._isSubset=False
            seen.append(newOrtho)
                #print(newOrtho.toString())
                #print(len(seen))
    #for s in range(0,int(len(seen)+1/2)):
    #    for t in range(int(len(seen)-1/2),len(seen)):
    #        seen[s].compareToOrthogroup(seen[t])
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
                print(str(s._id)+"\t"+str(redcnt)+"\t"+s.toString())
                red.write(str(s._id)+"\t"+s.toString())
                red.write("\n")
                redcnt+=1
            #print("WARNING", s, "not checked!")
    print("# warnings ", warn)
    print("# processed ", proc)

 


def filterRedundancyPlusPlus(grp, out, rename=False):
    seen = []
    idDct = dict()
    knownIDs = set()
    ok = 0
    redundant = 0
    isOK = True
    try:
        g = open(grp,'r')
    except IOError as e:
        print(e)
    for l in g:
            l = l.strip()
            checkThis = False
            tmp = l.split("\t")
            newOrtho = Orthogroup(tmp)
            todelete =None
            for k in tmp:
                if k in knownIDs:
                    checkThis = True
                #print(checkThis)
                knownIDs.add(k)
                if not k in idDct:
                    idDct[k] = [newOrtho]
                else:
                    idDct[k].append(newOrtho)
            if checkThis ==True:
                for s in idDct[k]:
                    s.compareToOrthogroup2(newOrtho)
                    if s._isSubset:
                        todelete =s
                    elif newOrtho._isSubset:
                        todelete=newOrtho
                if todelete:
                    try:
                        seen.remove(todelete)
                    except ValueError as e:
                        #print("#",e, todelete.toString())
                        pass
            if newOrtho._isSubset==False:
                seen.append(newOrtho)
                #print(newOrtho.toString())
                #print(len(seen))
    #for s in range(0,int(len(seen)+1/2)):
    #    for t in range(int(len(seen)-1/2),len(seen)):
    #        seen[s].compareToOrthogroup(seen[t])
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
                print(str(s._id)+"\t"+str(redcnt)+"\t"+s.toString())
                red.write(str(s._id)+"\t"+s.toString())
                red.write("\n")
                redcnt+=1
            #print("WARNING", s, "not checked!")
    print("# warnings ", warn)
    print("# processed ", proc)

 


def filterRedundancyPlus(grp, out):
    seen = []
    
    knownIDs = set()
    ok = 0
    redundant = 0
    isOK = True
    try:
        g = open(grp,'r')
    except IOError as e:
        print(e)
    for l in g:
            l = l.strip()
            checkThis = False
            tmp = l.split("\t")
            newOrtho = Orthogroup(tmp)
            todelete =None
            for k in tmp:
                if k in knownIDs:
                    checkThis = True
                print(checkThis)
                knownIDs.add(k)
            if checkThis ==True:
                for s in seen:
                    s.compareToOrthogroup2(newOrtho)
                    if s._isSubset:
                        todelete =s
                    elif newOrtho._isSubset:
                        todelete=newOrtho
                if todelete:
                    try:
                        seen.remove(todelete)
                    except ValueError as e:
                        print(e)
                        pass
            if newOrtho._isSubset==False:
                seen.append(newOrtho)
                print(newOrtho.toString())
                print(len(seen))
    #for s in range(0,int(len(seen)+1/2)):
    #    for t in range(int(len(seen)-1/2),len(seen)):
    #        seen[s].compareToOrthogroup(seen[t])
    warn = 0
    proc = 0
    for s in seen:
        if s._checked == False:
            warn+=1
           # print (s.toString())
        if s._checked == True:
            proc +=1
        
        if s._isSubset==False:
            print(s.toString())
    
            #print("WARNING", s, "not checked!")
    print("warnings ", warn)
    print("processed ", proc)
def filterRedundancy(grp, out):
    seen = set()
    ok = 0
    redundant = 0
    isOK = True
    try:
        g = open(grp,'r')
    except IOError as e:
        print(e)
    for l in g:
            l = l.strip()
            isOK = True
            tmp = l.split("\t")
            for t in tmp:
                #print(t)
                if t in seen:
                    isOK = False
                
                seen.add(t)
            if isOK:
                ok+=1
            else:
                    #print(tmp)
                redundant+=1
    print("# non redundant ",ok)
    print("# redundant ",redundant)  

def filterGroups(grp, out, numspec=None,rename=False):
    try:
        g = open(grp,'r')
    except IOError as e:
        print(e)
        usage()
     
    if not numspec:
        max = 0
        for l in g:
            l = l.strip()
            tmp = len(l.split("\t"))
            #print(tmp)
            if tmp > max:
                max = tmp
        g.close()
    numspec=max
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
        ln = l.strip()
        tmp = len(ln.split("\t"))
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
        elif o in ("-r", "--rename"):
            rename = True
        else:
            assert False, "unhandled option"
    ########################################
    if grp is None:
        usage()
    if out is None:
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
   # filterRedundancyPlusPlus(grp,outNoSubsets, rename)
    filterGroups(grp, "old", numspec, rename)
    unionOrthogroups(grp, outNoSubsets, 4, rename)
    filterGroups(outNoSubsets, out, numspec, rename=False)
    filterRedundancy(outNoSubsets,"tmp")
    #filterRedundancy(out,"tmp")

    
if __name__=="__main__":
    main()