import httplib2, sys,json
from decimal import Decimal
from datetime import datetime, date, timedelta
import mysql.connector
import os
import time
###################REST
def getSequences(stableids,outdir, specname):
    res = ""
    http = httplib2.Http(".cache")
    server = "http://beta.rest.ensembl.org"
    out = open(outdir+os.sep+specname+".fa",'w')
    for g in stableids:           
            ext = "/sequence/id/"+g
            resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"text/x-fasta"})
            print(content.decode("utf8"))
            out.write(content.decode("utf8"))
            time.sleep(0.4)
    out.close()
    
    
def getHomology(targetspecies, queryspecies, querystableids, outdir):
    out = open(outdir+os.sep+queryspecies+"_"+targetspecies+".tsv",'w')
    print("target:",targetspecies)
    print("query",queryspecies)
    print("stableids",querystableids)
    res = ""
    http = httplib2.Http(".cache")
    server = "http://beta.rest.ensembl.org"
    for g in  querystableids:           
        ext = "/homology/id/"+g+"?content-type=application/json;format=condensed;type=orthologues;target_species="+targetspecies
        print(ext)
        resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
        print(content.decode("utf8"))
        data = json.loads(content.decode("utf8"))
        print("data:",data)
        
        for f in data.values():
            try:
                #print("f:",f)
                #print("f0",f[0])
                #print(dict(f))
                tmp1 = f[0]["id"]
                tmp2 = f[0]["homologies"][0]
                tmp2 = tmp2["id"]
            #print(tmp1,tmp2)
        #print(g, dict(content.decode("utf8"))["id"] )
                out.write("\t".join([tmp1,tmp2]))
                out.write("\n")
            except IndexError as e:
                    pass
            time.sleep(0.4)
    
    
    
def getGeneProteinRelation( outdir, specname, release):
    stableids=[]
    http = httplib2.Http(".cache")
    server = "http://beta.rest.ensembl.org"
    cmd = 'show databases like "'+specname+'_core_'+str(release)+'_%";'
    cnx = mysql.connector.connect(user='anonymous', host='ensembldb.ensembl.org' )

    curA = cnx.cursor(buffered=True)
    curB = cnx.cursor(buffered=True)
    curA.execute(cmd)
    
    
    size =dict()
    gene = dict()
    parent = dict()
    
    
    for a in curA:
        cmd = "use "+a[0]+";"
        print(cmd)
        curB.execute(cmd)
        cmd = 'select stable_id from gene where biotype="protein_coding" and source ="ensembl" limit 3;'
        curB.execute(cmd)
        out = open(outdir+os.sep+specname+"_"+str(release)+".tsv",'w')
        seen = set()
       
        for b in curB:#tuples
            print(b)
            for bb in b:
                time.sleep(0.4)
                ext = "/feature/id/"+bb+"?feature=cds"
                resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
                data = json.loads(content.decode("utf8"))
                for f in data:
                    try:
                        size[f["ID"]]+=f["end"]-f["start"]+1
                        parent[bb].add(f["ID"])
                    except KeyError as e:
                        size[f["ID"]]=f["end"]-f["start"]+1
                        gene[f["ID"]]=bb
                    try:
                        parent[bb].add(f["ID"])
                    except KeyError as e:
                        parent[bb]=set()
    for g in parent:
        for gg in parent[g]:
            stableids.append(gg)
            out.write("\t".join([g,gg,str(size[gg])]))
            out.write("\n")
            #print(g,gg,size[gg] )
    curA.close()
    curB.close()
    return(stableids, parent.keys())

def useEnsemblDB(specs,rel):#, targets):
    res = ""
    http = httplib2.Http(".cache")
    server = "http://beta.rest.ensembl.org"
    #cmd = 'show databases like "'+spec+'_core_'+rel+'_%";'
    #print(cmd)
    
    for i in range(0,len(specs)):
        print(i)
        stableids, geneids = getGeneProteinRelation( "./", specs[i], rel[i])
        getSequences(stableids,  ".",specs[i] )
        if(i+1<len(specs)):
                getHomology(specs[i+1],specs[i],geneids, "./")
def specInfo():
    pass
    http = httplib2.Http(".cache") 
    server = "http://beta.rest.ensembl.org"
    ext = "/info/species"
    resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
    data = json.loads(content.decode("utf8"))
    return(data)
def selectSpecies(data, namelist):
    backupdct={}
    specs = []
    rels =[]
    for s in data["species"]:
        print(s["release"], s["name"], s["aliases"], s["groups"])
        if (s["name"] in namelist ):
            specs.append(s["name"])
            rels.append(s["release"])
        else:
            aliases = s["aliases"]
            print(aliases)
            #aliases = aliases.split(",")
            found = [a for a in aliases if a in namelist]
            if found:
                specs.append(s["name"])
                rels.append(s["release"])
    return(specs, rels)
def main():
    namelist = ["bos taurus", "human", "chimp"]
    data = specInfo()
    specs, rels = selectSpecies(data, namelist)
    #for i in ind:
    #   tmp = ""
    useEnsemblDB(specs, rels)
        #print(tmp)
  
if __name__ == "__main__":
    main()