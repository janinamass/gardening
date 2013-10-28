import mysql.connector

def q(s):
    return ('"'+s+'"')
def getGenome_Db_Ids(specnames, release=71):
    #stableids=[]
    #http = httplib2.Http(".cache")
    #server = "http://beta.rest.ensembl.org"
    #cmd = 'show databases like "'+specname+'_core_'+str(release)+'_%";'
    cmd = 'use ensembl_compara_'+str(release)+";" 
    cnx = mysql.connector.connect(user='anonymous', host='ensembldb.ensembl.org' )
    curA = cnx.cursor(buffered=True)
    print (cmd)
    #curB = cnx.cursor(buffered=True)
    curA.execute(cmd)
    specnames = [q(s) for s in specnames]
    s = " OR name=".join(specnames)
    cmd = 'SELECT genome_db_id FROM genome_db WHERE name='+s+";"
    print(cmd)
    curA.execute(cmd)
    res = curA.fetchall()
    print(res)
    return(res)


res = getGenome_Db_Ids(["homo_sapiens", "gorilla_gorilla"])
res = [r[0] for r in res]
print(res)

def getSpecies_Set_Ids(specids, release=71):
    pass