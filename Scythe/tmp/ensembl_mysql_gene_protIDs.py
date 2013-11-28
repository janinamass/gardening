import mysql.connector

def q(s):
    return ('"'+s+'"')

def ensembl_gene_protein(database):
    cnx = mysql.connector.connect(user='anonymous', host='ensembldb.ensembl.org' )
    curA = cnx.cursor(buffered=True)
    cmd = 'use '+database+';' 
    curA.execute(cmd)
    selectcmd = 'SELECT DISTINCT g.stable_id, trn.stable_id, trs.seq_region_start,trs.seq_region_end \
    FROM gene g, transcript trs, translation trn \
    WHERE g.gene_id = trs.gene_id AND trs.transcript_id = trn.transcript_id LIMIT 50;'
    curA.execute(selectcmd)
    res = curA.fetchall()
    curA.close()
    return(res)

def ensembl_gene_protein_canonical(database):
    cnx = mysql.connector.connect(user='anonymous', host='ensembldb.ensembl.org' )
    curA = cnx.cursor(buffered=True)
    cmd = 'use '+database+';' 
    curA.execute(cmd)
    selectcmd = 'SELECT DISTINCT g.stable_id, trs.stable_id, trn.stable_id \
    FROM gene g, transcript trs, translation trn \
    WHERE g.gene_id = trs.gene_id AND trs.transcript_id = trn.transcript_id LIMIT 5;'
    curA.execute(selectcmd)
    res = curA.fetchall()
    curA.close()
    return(res)


def getDatabaseId(specname, release=73):
    cnx = mysql.connector.connect(user='anonymous', host='ensembldb.ensembl.org' )
    curA = cnx.cursor(buffered=True)
    cmd = 'SHOW DATABASES like "'+specname+'_core_'+str(release)+'_%";'
    curA.execute(cmd)
    res = curA.fetchone()[0]
    curA.close()
    return(res)


species = ["homo_sapiens","mus_musculus","gorilla_gorilla"]
for s in species:
    db = getDatabaseId(s)
    print(db)
    tmpfile = open(db+"gene_tr.tmp.tsv",'w')
    #tmpfile = open(db+"gene_tr.tsv",'w')
    res = ensembl_gene_protein(db)
    for r in res:
        r = [str(t) for t in r]
        s = "\t".join(r)
        s += "\t"+str(int(r[-1])-int(r[-2])+1)
        tmpfile.write(s)
        tmpfile.write("\n")
    
    tmpfile.close()
    
    