import httplib2, sys,json
from decimal import Decimal
from datetime import datetime, date, timedelta
import mysql.connector
import os
import time

def getSequences(stableids,outdir, specname):
    res = ""
    http = httplib2.Http(".cache")
    server = "http://beta.rest.ensembl.org"
    path=outdir+os.sep+"fa"
    if not os.path.isdir(path):
            os.makedirs(path)
    print(path)
    out = open(outdir+os.sep+"fa"+os.sep+specname+".fa",'w')
    for g in stableids:           
            ext = "/sequence/id/"+g
            resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"text/x-fasta"})
            print(content.decode("utf8"))
            out.write(content.decode("utf8"))
            time.sleep(0.4)
    out.close()
    
    
def getHomology(targetspecies, queryspecies, querystableids, outdir):
    path=outdir+os.sep+"ensembl_ortho_tsv"
    if not os.path.isdir(path):
            os.makedirs(path)
    out = open(outdir+os.sep+"ensembl_ortho_tsv"+os.sep+queryspecies+"__"+targetspecies+".tsv",'w')
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
                tmp1 = f[0]["id"]
                tmp2 = f[0]["homologies"][0]
                print(tmp2)
                if (tmp2["type"] == "ortholog_one2one" ):
                    tmp2 = tmp2["id"]
                    out.write("\t".join([tmp1,tmp2]))
                    out.write("\n")
                
            except IndexError as e:
                    pass
            time.sleep(0.2)
    
    
    
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
        cmd = 'select stable_id from gene where biotype="protein_coding" and source ="ensembl" limit 10;'
        curB.execute(cmd)
        path=outdir+os.sep+"ensembl_gene_tsv"
        if not os.path.isdir(path):
            os.makedirs(path)
        out = open(outdir+os.sep+"ensembl_gene_tsv"+os.sep+specname+"_"+str(release)+".tsv",'w')
        seen = set()
       
        for b in curB:#tuples
            print(b)
            for bb in b:
                time.sleep(0.2)
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

def useEnsemblDB(specs,rel, outdir):#, targets):
    res = ""
    http = httplib2.Http(".cache")
    server = "http://beta.rest.ensembl.org"
    
    for i in range(0,len(specs)):
        print(i)
        stableids, geneids = getGeneProteinRelation( outdir, specs[i], rel[i])
        getSequences(stableids,  outdir,specs[i] )
        if(i+1<len(specs)):
                getHomology(specs[i+1],specs[i],geneids, outdir)
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
            found = [a for a in aliases if a in namelist]
            if found:
                specs.append(s["name"])
                rels.append(s["release"])
    return(specs, rels)
def main():
    namelist = ["bos taurus", "human", "chimp"]
    data = specInfo()
    outdir=None
    if not outdir:
        outdir = "./"
    specs, rels = selectSpecies(data, namelist)
    useEnsemblDB(specs, rels, outdir)
  
if __name__ == "__main__":
    main()