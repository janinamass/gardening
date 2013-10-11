import httplib2, sys,json
 
http = httplib2.Http(".cache")
 
server = "http://beta.rest.ensembl.org"
ext = "/sequence/id/ENSG00000157764?"

ext = "/feature/id/ENSMUSG00000020707?feature=cds"


ext = "/feature/id/ENSMUSG00000021712?feature=cds"


resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
 
if not resp.status == 200:
    print ("Invalid response: ", resp.status)
    #sys.exit()
data = json.loads(content.decode("utf8"))
print(content)
print(data)
data2 = [dict(d) for d in data]
for f in data2:
    print(f["ID"])