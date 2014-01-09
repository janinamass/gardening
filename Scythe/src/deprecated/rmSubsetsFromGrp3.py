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
        self._ignore = False
    def expand(self, other):
        tmpset = set()
        tmplist = []
        dct = {}
        purge= False
        for i in self._list:
            tmplist.append(i)
            if i == "ENSBTAG00000000417":
                        print(i,dct)
                        #print( i != dct[i[0:6]])
            if i[0:6] in dct:
                if i != dct[i[0:6]]:
                    
                    print("Shit.",dct[i[0:6]],i)
                    purge=True
                    #return(None)
            else:
                dct[i[0:6]]=i
        for j in other._list:
            if j == "ENSBTAG00000000417":
                        print(j,dct)
                        #print( j != dct[j[0:6]])
            if j[0:6] in dct:
                if j != dct[j[0:6]]:
                    print("Shit.",dct[j[0:6]],j)
                    purge=True
                    #return(None)
            else:
                dct[j[0:6]]=j
            tmplist.append(j)
        tmpset = set(tmplist)
        newlist=list(tmpset)
        newid= str(self._id)+"_"+str(other._id)
        
        #print( newid, len(set(self._list))-len(set(other._list)), len(set(self._list).intersection(set(other._list)))) 
        #print( len(tmpset.intersection(set(self._list)))-len(set(other._list))) #last resort: compare first 7letters
        initlist=[newid]
        for i in newlist:
            initlist.append(i)
        new = Orthogroup(initlist)
        if purge:
            new._isSubset=True
            new._ignore=True
        #new._id = str(self._id)+"_"+str(other._id)
        return(new)
    def toString(self):
        s = "\t".join(self._list)
        return(s)
    


 
def checkConsistency(grp, outfile, rename=False):
    knownGenes= set()
    kglist = []
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
            newOrtho= Orthogroup(tmp)
            for t in tmp[1:]:                
                if t in knownGenes:
                    tmplist = []
                    for o in godict[t]:
                        new = newOrtho.expand(o)
                        new._isSubset=False
                        o._isSubset = True
                        newOrtho._isSubset = True
                        tmplist.append(new)
                    for l in tmplist:
                        godict[t].append(l)
                else:
                    knownGenes.add(t)
                    kglist.append(t)
                    godict[t] = [newOrtho]
    g.close()
    try:
        g = open(grp,'r')
        #o = open(outfile,'w')
    except IOError as e:
        print(e)
    processed =[]
    for l in g:
        l = l.rstrip()
        tmp = l.split("\t")
        for t in tmp[1:]:        
            nmm = [t for t in godict[t] if  t._isSubset==False and t._ignore==False]
            out=None
            if nmm !=[]:
                out = nmm[-1]
                if out._id not in proc:
                    processed.append(out._id+"\t"+out.toString())
                   # o.write(out._id+"\t"+out.toString())
                   # o.write("\n")
                proc.add(out._id)
    g.close()
    #o.close()
    try:
        o = open(outfile,'w')
    except IOError as e:
        print(e)
    
    seen = set()
    cnt = 0
    ok = 0
    redundant = 0
    isOK = True
    #print(processed)
    for p in processed:
        isOK = True
        tmp = p.split("\t")
        #print(tmp)
        for t in tmp[1:]:
            if t in seen:
                isOK = False
            seen.add(t)
        if isOK:
            ok+=1
            if not rename:
                o.write(str(p))
            if rename:
                o.write(str(cnt)+"\t"+"\t".join(p.split("\t")[1:]))
            o.write("\n")
            cnt+=1
        else:
            #print("RED",tmp)
            redundant+=1
    print("# non redundant ",ok)
    print("# redundant ",redundant)      
    o.close()
    o = open(outfile+"_bla",'w')
    for k in kglist:
        o.write(k)
        o.write("\n")
    o.close()
#===============================================================================
# def filterRedundancy(grp, out):
#     seen = set()
#     ok = 0
#     redundant = 0
#     isOK = True
#     try:
#         g = open(grp,'r')
#         o = open(out,'w')
#     except IOError as e:
#         print(e)
#     for l in g:
#             l = l.rstrip()
#             isOK = True
#             tmp = l.split("\t")
#             for t in tmp[1:]:
#                 #print(t)
#                 if t in seen:
#                     isOK = False
#                 
#                 seen.add(t)
#             if isOK:
#                 ok+=1
#                 o.write(l)
#                 o.write("\n")
#             else:
#                 #print("RED",tmp)
#                 redundant+=1
#     print("# non redundant ",ok)
#     print("# redundant ",redundant)  
#     g.close()
#     o.close()
#===============================================================================
def filterGroups(grp, outfile, numspec=None,rename=False):
    kglist = []
    
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
        out = open(outfile,'w')
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
            for m in ln.split("\t")[1:]:
                #print(m)
                kglist.append(m)
    out.close()
    #return(numspec)
    o = open(outfile+"_bla",'w')
    print("bla")
    for k in kglist:
        o.write(k)
        o.write("\n")
    o.close()
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
    
    numspec= filterGroups(grp,"tmp1", numspec, rename)
    checkConsistency(grp,"tmp2",rename)
                     #outNoSubsets, rename)
    filterGroups("tmp2", "tmp3", numspec, rename)
     #            outNoSubsets, out, numspec, rename)
    #filterRedundancy("tmp", outNoSubsets)
    #filterRedundancy(outNoSubsets,"tmp")

    #filterRedundancyPlusPlus(grp,outNoSubsets, rename)
    #unionOrthogroups(grp, outNoSubsets, numspec, rename)
    #filterGroups(outNoSubsets, out, numspec, rename=False)
    #filterRedundancy(outNoSubsets,"tmp")
    #filterRedundancy(out,"tmp")

    
if __name__=="__main__":
    main()