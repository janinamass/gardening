import httplib2, sys,json
from decimal import Decimal
from datetime import datetime, date, timedelta
import mysql.connector
import os
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
    out.close()
    
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
        curB.execute(cmd)
        cmd = 'select stable_id from gene where biotype="protein_coding" and source ="ensembl" and stable_id like "ENSG%" limit 20;'
        curB.execute(cmd)
        out = open(outdir+os.sep+specname+"_"+str(release)+".tsv",'w')
        seen = set()
       
        for b in curB:#tuples
            for bb in b:
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
    return(stableids)

def useEnsemblDB(specs,rel):#, targets):
    res = ""
    http = httplib2.Http(".cache")
    server = "http://beta.rest.ensembl.org"
    #cmd = 'show databases like "'+spec+'_core_'+rel+'_%";'
    #print(cmd)
    
    for i in range(0,len(specs)):
        print(i)
        stableids = getGeneProteinRelation( "./", specs[i], rel[i])
        getSequences(stableids,  ".",specs[i] )
    
    #curA = cnx.cursor(buffered=True)
    #curB = cnx.cursor(buffered=True)
    #curA.execute(cmd)
    #for a in curA:
        #print(a)
     #   cmd = "use "+a[0]+";"
     #   curB.execute(cmd)
        #print("exec curb")
        #print(curB)
     #   cmd = 'select stable_id from gene where biotype="protein_coding" and source ="ensembl" and stable_id like "ENSG%" limit 2;'
     #   curB.execute(cmd)
     #   seen = set()
     #   size =dict()
     #   gene = dict()
     #   parent = dict()
     #   for b in curB:#tuples
     #       for bb in b:
     #           ext = "/feature/id/"+bb+"?feature=cds"
     #           resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
     #           data = json.loads(content.decode("utf8"))
     #           for f in data:
     #               try:
     #                   size[f["ID"]]+=f["end"]-f["start"]+1
     #                   parent[bb].add(f["ID"])
     #               except KeyError as e:
     #                   size[f["ID"]]=f["end"]-f["start"]+1
     #                   gene[f["ID"]]=bb
     #               try:
     #                   parent[bb].add(f["ID"])
     #               except KeyError as e:
     #                   parent[bb]=set()
    #for g in parent:
    #    for gg in parent[g]:
    #        print(g,gg,size[gg] )
    #        
    #        #server = "http://beta.rest.ensembl.org"
    #        ext = "/sequence/id/"+gg
    #        resp1, content1 = http.request(server+ext, method="GET", headers={"Content-Type":"text/x-fasta"})
    #        print(content1.decode("utf8"))
    #        
    #    ext = "/homology/id/"+g+"?format=condensed;type=orthologues;target_species="+targets[0]
    #    resp2, content2 = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
    #        #if not resp.status == 200:
    #        #print "Invalid response: ", resp.status
    #        #sys.exit()
    #        #import json
# 
#        decoded = json.loads(content2.decode("utf8"))
#        
#        for d in decoded["data"]:
#            print(g,d["id"])

  #ext = "/homology/id/ENSG00000157764?target_taxon=10090;sequence=cdna;target_species=cow;type=orthologues"
            
            
            
    #for g in gene:
    #    print(g,gene[g],size[g] )
   # curA.close()
   # curB.close()
  #  return (res)
def main():
    http = httplib2.Http(".cache") 
    server = "http://beta.rest.ensembl.org"
    ext = "/info/species"
    resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
    data = json.loads(content.decode("utf8"))

    specs =[]
    rels = []
    targets = []
    for s in data["species"]:
        print(s["release"], s["name"], s["aliases"], s["groups"])
        specs.append(s["name"])
        rels.append(s["release"])
    
    
    ind =[34]
    specs = [specs[2],specs[34]]
    #for i in ind:
    #   tmp = ""
    useEnsemblDB(specs, rels)
        #print(tmp)
  
if __name__ == "__main__":
    main()