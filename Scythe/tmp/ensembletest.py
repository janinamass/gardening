import httplib2, sys,json
from decimal import Decimal
from datetime import datetime, date, timedelta
import mysql.connector

###################REST
def useEnsemblDB(spec,rel):
    res = ""
    http = httplib2.Http(".cache")
    server = "http://beta.rest.ensembl.org"
    cmd = 'show databases like "'+spec+'_core_'+rel+'_%";'
    print(cmd)
    cnx = mysql.connector.connect(user='anonymous', host='ensembldb.ensembl.org' )
    
    curA = cnx.cursor(buffered=True)
    curB = cnx.cursor(buffered=True)
    curA.execute(cmd)
    for a in curA:
        print(a)
        cmd = "use "+a[0]+";"
        curB.execute(cmd)
        print("exec curb")
        print(curB)
        cmd = 'select stable_id from gene where biotype="protein_coding" and source ="ensembl" and stable_id like "ENSG%" limit 5;'
        curB.execute(cmd)
        seen = set()
        res=""
        prev=""
        cnt =0
        for b in curB:
            #print(b)
            
            for bb in b:
                #seen = set()
                size =dict()
                tmpstr = ""
                print(bb)
                ext = "/feature/id/"+bb+"?feature=cds"
                resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
                data = json.loads(content.decode("utf8"))
                #print(data)
                #size needs to be added up and spit out at the end
                #print(data)
                for f in data:
                    if prev=="": # very first
                        prev=f["ID"]
                        #print("first",prev )
                        #size[f["ID"]]=f["end"]-f["start"]+1
                        #seen.add(f["ID"]) 
                    if not f["ID"] in seen: # of protein
                        #print("not yet seen",f["ID"],bb )
                        size[f["ID"]]=f["end"]-f["start"]+1
                        seen.add(f["ID"])        #print out last
                        if prev in size: #previous existed
                            #print(bb,prev+"\t"+str(size[prev]))
                            #res +=bb,prev+"\t"+str(size[prev])+"\n"
                            res +="\t" .join([bb,prev,str(size[prev]),"\n"])
                        else:
                            size[f["ID"]]+=f["end"]-f["start"]+1
                            #res +="\n"
                    else: #already seen
                        try:
                            size[f["ID"]]+=f["end"]-f["start"]+1
                        except KeyError as k:
                            print(k,bb,f)#,data)
                        #print(f )
                    prev=f["ID"] #set new prev
            #print(bb,prev+"\t"+str(size[prev])+"\t"+str(size[prev]/3))
            try:
                res +="\t".join([bb,prev,str(size[prev]),"\n"])
            except KeyError as k:
                print(k,bb)#,data)
                    #dont forget the last gene
    curA.close()
    curB.close()
    return (res)
def main():
    http = httplib2.Http(".cache") 
    server = "http://beta.rest.ensembl.org"
    ext = "/sequence/id/ENSG00000157764?"
    ext = "/feature/id/ENSMUSG00000020707?feature=cds"    
    ext = "/feature/id/ENSMUSG00000021712?feature=cds"
    ext = "/info/species"
    
    resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
    # 
    #if not resp.status == 200:
    #    print ("Invalid response: ", resp.status)
    #    #sys.exit()
    data = json.loads(content.decode("utf8"))
    #print(content)
    #print(data)
    #data = data.split(",")
    for d,e in data.items():
        try:
            #dt= dict(d)
            print(d,e)#["release"], d["name"])
            print("\n")
        except Error as e:
            print(e)
    print(data["species"])
    specs =[]
    rels = []
    for s in data["species"]:
        print(s["release"], s["name"], s["aliases"], s["groups"])
        specs.append(s["name"])
        rels.append(s["release"])
    #data2 = [dict(d) for d in data]
    #for f in data2:
    #    print(f)
    #    print(f["ID"])
    ##########################################
    #cnx = mysql.connector.connect(user='anonymous', host='ensembldb.ensembl.org' )
    
    #curA = cnx.cursor(buffered=True)
    #curB = cnx.cursor(buffered=True)
    
    ind =[34]
    #spec = "pan_troglodytes"
    #spec = specs[ind]
    #rel = str(rels[ind])
    for i in ind:
        tmp = useEnsemblDB(specs[i], str(rels[i]))
        print(tmp)
    #cmd = 'show databases like "'+spec+'_core_'+rel+'_%";'
    #print(cmd)
    #curA.execute(cmd)
    #for a in curA:
    #    print(a)
    #    cmd = "use "+a[0]+";"
    #    curB.execute(cmd)
    #    cmd = "select stable_id from gene limit 10;"
    #    curB.execute(cmd)
    #    for b in curB:
    #        print(b)
    #        ext = "/feature/id/"+b[0]+"?feature=cds"
    
    #resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
     
    #if not resp.status == 200:
    #    print ("Invalid response: ", resp.status)
        #sys.exit()
    #data = json.loads(content.decode("utf8"))
    #print(data)
    #data = data[0]
    #print(data["Parent"], data["ID"], data["rank"])
    #for d,e in data.items():
    #    try:
            #dt= dict(d)
    #        print(d,e)#["release"], d["name"])
    #        print("\n")
    #    except Error as e:
    #        print(e)
if __name__ == "__main__":
    main()