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
    #s = " OR name=".join(specnames)
    res = dict()
    for s in specnames:
        cmd = 'SELECT genome_db_id FROM genome_db WHERE name='+s+";"
        print(cmd)
        curA.execute(cmd)
        res[s] = curA.fetchall()[0][0]
    print(res)
    return(res)




def getSpecies_Set_Ids(genomedb_ids, release=71):
    cmd = 'use ensembl_compara_'+str(release)+";" 
    cnx = mysql.connector.connect(user='anonymous', host='ensembldb.ensembl.org' )
    curA = cnx.cursor(buffered=True)
    print (cmd)
    curA.execute(cmd)
    genomedb_ids=[str(s) for s in genomedb_ids]
    s = "genome_db_id ="+" OR genome_db_id=".join(genomedb_ids)
    cmd = "SELECT species_set_id,genome_db_id, count(species_set_id )"
    cmd += "from species_set  where ("
    cmd += s
    cmd += ") group by species_set_id having count(species_set_id) > 1  ;"
    print (cmd)
    curA.execute(cmd)
    res = curA.fetchall()
    return (res)

def getMethodLinkSpecies_Set_Ids(speciesSet_ids, release=71):
    cmd = 'use ensembl_compara_'+str(release)+";" 
    cnx = mysql.connector.connect(user='anonymous', host='ensembldb.ensembl.org' )
    curA = cnx.cursor(buffered=True)
    print (cmd)
    curA.execute(cmd)
    speciesSet_ids=[str(s) for s in speciesSet_ids]
    s = "species_set_id ="+" OR species_set_id=".join(speciesSet_ids)
    cmd = "SELECT method_link_species_set_id, species_set_id"
    cmd += " FROM method_link_species_set  WHERE "
    cmd += s
    cmd += " ;"
    print (cmd)
    curA.execute(cmd)
    res = curA.fetchall()
    curA.close()
    return (res)

def getHomologyId(method_link_species_set_ids, release=71):
    cmd = 'use ensembl_compara_'+str(release)+";" 
    cnx = mysql.connector.connect(user='anonymous', host='ensembldb.ensembl.org' )
    curA = cnx.cursor(buffered=True)
    print (cmd)
    curA.execute(cmd)
    method_link_species_set_ids=[str(s) for s in method_link_species_set_ids]
    s = "method_link_species_set_id ="+" OR method_link_species_set_id=".join(method_link_species_set_ids)
    cmd = "SELECT homology_id, method_link_species_set_id"
    cmd += " FROM homology WHERE ("
    cmd += s
    cmd += ");"
    print (cmd)
    curA.execute(cmd)
    res = curA.fetchall()
    print(res)
    curA.close()
    return (res)

def getHomologyMemberId(homology_ids, release=71):
    cmd = 'use ensembl_compara_'+str(release)+";" 
    cnx = mysql.connector.connect(user='anonymous', host='ensembldb.ensembl.org' )
    curA = cnx.cursor(buffered=True)
    print (cmd)
    curA.execute(cmd)
    homology_ids=[str(s) for s in homology_ids]
    s = "homology_id ="+" OR homology_id=".join(homology_ids)
    cmd = "SELECT homology_member.member_id, member.member_id,member.stable_id, homology_id "
    cmd += " FROM homology_member LEFT JOIN  member ON "
    cmd+= "(homology_member.member_id = member.member_id) WHERE ("
    cmd += s
    cmd += ");"
    print (cmd)
    curA.execute(cmd)
    res = curA.fetchall()
    print(res)
    curA.close()
    return (res)
 #SELECT * FROM genome_db LEFT JOIN species_set on species_set_id LEFT JOIN
 # method_link_species_set ON species_set.species_set_id LEFT JOIN homology 
 # ON method_link_species_set.method_link_species_set_id LEFT JOIN 
 # homology_member on homology.homology_id LEFT JOIN member on 
 # homology_member.member_id WHERE ( genome_db.genome_db_id=90 
#AND homology.description="ortholog_one2one") limit 5;

# Select homology.homology_id,member.member_id,member.stable_id from 
# homology left join homology_member on homology_member.homology_id=homology.homology_id left join member on (homology_member.member_id=member.member_id ) where  method_link_species_set_id=28216  limit 6;









def test(method_link_species_set_ids,release = 71):
    method_link_species_set_ids=[str(s) for s in method_link_species_set_ids]

    cmd = 'use ensembl_compara_'+str(release)+";" 
    cnx = mysql.connector.connect(user='anonymous', host='ensembldb.ensembl.org' )
    curA = cnx.cursor(buffered=True)
    print (cmd)
    curA.execute(cmd)


    s = "method_link_species_set.method_link_species_set_id ="+" OR method_link_species_set.method_link_species_set_id=".join(method_link_species_set_ids)

    
    cmd = " SELECT gene_member_id,stable_id"
    cmd +=" FROM genome_db "
    cmd+=" LEFT JOIN species_set ON species_set_id"
    cmd +=" LEFT JOIN method_link_species_set ON "
    cmd +=" (species_set.species_set_id=method_link_species_set.species_set_id) "
    cmd +=" LEFT JOIN homology ON "
    cmd +=" (method_link_species_set.method_link_species_set_id = homology.method_link_species_set_id) "
    cmd +=" LEFT JOIN homology_member ON (homology.homology_id = homology_member.homology_id)"
    cmd +=" LEFT JOIN member on (homology_member.member_id = member.member_id)"
    cmd +=" WHERE ("
    #genome_db.genome_db_id=90 AND 
    cmd +='homology.description="ortholog_one2one")'
    cmd +=" AND ("+s
    cmd +=") limit 6;"

    print (cmd)
    curA.execute(cmd)
    res = curA.fetchall()
    print(res)
    curA.close()
    return (res)

def test2(method_link_species_set_ids,release = 71):
    method_link_species_set_ids=[str(s) for s in method_link_species_set_ids]

    cmd = 'use ensembl_compara_'+str(release)+";" 
    cnx = mysql.connector.connect(user='anonymous', host='ensembldb.ensembl.org' )
    curA = cnx.cursor(buffered=True)
    print (cmd)
    curA.execute(cmd)


    s = "method_link_species_set_id ="+" OR method_link_species_set_id=".join(method_link_species_set_ids)

    
    cmd = " Select homology.homology_id,member.member_id,member.stable_id from"
    cmd +=" homology left join homology_member on "
    cmd+=" homology_member.homology_id=homology.homology_id "
    cmd +=" left join member on (homology_member.member_id=member.member_id ) "
    #cmd +=" where  method_link_species_set_id=28216  "
    cmd +=" WHERE ("
    #genome_db.genome_db_id=90 AND 
    cmd +='homology.description="ortholog_one2one")'
    cmd +=" AND (method_link_species_set_id= "+s
    cmd +=") limit 60;"

    print (cmd)
    curA.execute(cmd)
    res = curA.fetchall()
    print(res)
    curA.close()
    return (res)

# Select homology.homology_id,member.member_id,member.stable_id from 
# homology left join homology_member on 
#homology_member.homology_id=homology.homology_id 
#left join member on (homology_member.member_id=member.member_id ) 
#where  method_link_species_set_id=28216  limit 6;
def test3(method_link_species_set_ids,release = 71):
    method_link_species_set_ids=[str(s) for s in method_link_species_set_ids]

    cmd = 'use ensembl_compara_'+str(release)+";" 
    cnx = mysql.connector.connect(user='anonymous', host='ensembldb.ensembl.org' )
    curA = cnx.cursor(buffered=True)
    print (cmd)
    curA.execute(cmd)


    s = "method_link_species_set.method_link_species_set_id ="+" OR method_link_species_set.method_link_species_set_id=".join(method_link_species_set_ids)

    
    cmd = " Select species_set.genome_db_id, homology.homology_id,member.member_id,member.stable_id from"
    cmd +=" homology left join homology_member on "
    cmd+=" homology_member.homology_id=homology.homology_id "
    cmd +=" left join member on (homology_member.member_id=member.member_id ) "
    
    cmd+=" left join method_link_species_set ON "
    cmd += "(homology.method_link_species_set_id = method_link_species_set.method_link_species_set_id) "
    cmd+=" LEFT JOIN species_set ON "
    cmd += "(species_set.species_set_id = method_link_species_set.species_set_id) "

    #cmd +=" where  method_link_species_set_id=28216  "
    cmd +=" WHERE ("
    #genome_db.genome_db_id=90 AND 
    cmd +='homology.description="ortholog_one2one")'
    cmd +=" AND (method_link_species_set.method_link_species_set_id= "+s
    cmd +=");"

    print (cmd)
    curA.execute(cmd)
    res = curA.fetchall()
    print(res)
    curA.close()
    return (res)

def test4(method_link_species_set_ids,release = 71):
    method_link_species_set_ids=[str(s) for s in method_link_species_set_ids]

    cmd = 'use ensembl_compara_'+str(release)+";" 
    cnx = mysql.connector.connect(user='anonymous', host='ensembldb.ensembl.org' )
    curA = cnx.cursor(buffered=True)
    print (cmd)
    curA.execute(cmd)


    s = "method_link_species_set.method_link_species_set_id ="+" OR method_link_species_set.method_link_species_set_id=".join(method_link_species_set_ids)

    
    cmd = " Select homology.homology_id,  member.genome_db_id,  member.stable_id from"
    cmd +=" homology left join homology_member on "
    cmd+=" homology_member.homology_id=homology.homology_id "
    cmd +=" left join member on (homology_member.member_id=member.member_id ) "
    
    cmd+=" left join method_link_species_set ON "
    cmd += "(homology.method_link_species_set_id = method_link_species_set.method_link_species_set_id) "
    cmd+=" LEFT JOIN species_set ON "
    cmd += "(species_set.species_set_id = method_link_species_set.species_set_id) "

    #cmd +=" where  method_link_species_set_id=28216  "
    cmd +=" WHERE ("
    #genome_db.genome_db_id=90 AND 
    cmd +='homology.description="ortholog_one2one")'
    cmd +=" AND (method_link_species_set.method_link_species_set_id= "+s
    cmd +=");"

    print (cmd)
    curA.execute(cmd)
    res = curA.fetchall()
    print(res)
    curA.close()
    return (res)

def makeTable(ids, orthologs2, orthoids):
    print(orthologs2)
    res = dict()
    specset = dict()
    #current = None
    #print(ids)
    #for o in orthoids:
    #    for i in ids:
    #        current = i
    #        if (o,i) in orthologs2:
    #            
    #            if orthologs2[(o,i)][0] not in mapping:
    #                mapping[orthologs2[(o,i)][0]]= [(o)]
    #            else:
    #                mapping[orthologs2[(o,i)][0]].append((o))
    #print(mapping)
    #seen = []
    #for m in mapping:
    #    for homology_id in mapping[m]:
    #        #print(homology_id)
    #        for i in ids:
    #            if (homology_id,i) in orthologs2:
    #                if orthologs2[(homology_id,i)][0] not in seen:
    #                    print(orthologs2[(homology_id,i)][0])
    #                    seen.append(orthologs2[(homology_id,i)][0])
                #else:
                #    print(" ")
    #    print("#\n")
    #seen = dict()
    for o in orthoids:
        for i in ids:
            if (o,i) in orthologs2:
                if not o in res:
                    res[o]=[i]
                else:
                    res[o].append(i)
    #            res+=orthologs2[(o,i)][0]+"\t"
    #            seen[]
    #        else:
    #            res+="-"+"\t"
    #for i in ids:
    #    res+="\n"
    #for o in orthologs2:
    #    print(o, orthologs2[o])
    #spec[orthologs[o][0][0]]=o
    
    print(res)
    for r in res:
        if not ((res[r][0],res[r][1])) in specset:
                 specset[(res[r][0],res[r][1])]=[(orthologs2[(r,res[r][0])],orthologs2[(r,res[r][1] )]) ]
    
        else:
             specset[(res[r][0],res[r][1])].append((orthologs2[(r,res[r][0])],orthologs2[(r,res[r][1])]))
             #specset[(res[r][0],res[r][1])].append(orthologs2[(r,res[r][1])])
    
    print(specset) 
    for s in specset:
        st = str(s[0])+"_"+str(s[1])
        print(st)
        out = open(st,"w")
        tuplist = specset[s]
        for tup in tuplist:
            tup = str(tup[0])+"\t"+str(tup[1])
            print(tup)
            out.write(tup)
            out.write("\n")

orthoids = [] # key: homology_id, value:list of tuples (species_id,stable_id) 
orthologs2 = dict()
members = dict()
res = getGenome_Db_Ids(["homo_sapiens", "gorilla_gorilla", "pan_troglodytes"])
#res = [r[0] for r in res]
print(res)
genomedbIDs = res
res = res.values()
#res = [r[0] for r in res]
###

res2 = getSpecies_Set_Ids(res,71)
print(res2)
res2 = [r[0] for r in res2]
print(res2)
res3 = getMethodLinkSpecies_Set_Ids(res2, 71)
print(res3)
res3 = [r[0] for r in res3]
#print(res3
restest= test4(res3,release = 71)
print(restest)
for r in restest:
    if not (r[0],r[1]) in orthologs2:
        if not r[0] in orthoids:
            orthoids.append(r[0])
       # orthologs[r[0]] = [r[1:]]
        orthologs2[(r[0],r[1])]=r[2:][0]
   # else:
   #     if r[1:] not in orthologs[r[0]]:
   #         orthologs[r[0]].append(r[1:])
   #         orthologs2[(r[0],r[1])]=r[2:][0]
   # if not r[1] in members:
   #     members[r[1]] = [(r[0],r[2:])]
   # else:
   #     members[r[1]].append((r[0],r[2:]))
#print(orthologs)
print(orthologs2)

makeTable(res, orthologs2,orthoids)

                
#print (members)
#somehow format this right...

#res3 = [r[0] for r in res3]
#print(res3)
#res4 =getHomologyId(res3,71)
#print(res4)
#res4 = [r[0] for r in res4]
#res5 =getHomologyMemberId(res4,71)
#print(res5)
#res4 = test(res3, 71)
#print(res4)


