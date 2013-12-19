from ftplib import FTP
import os
import tarfile
import zipfile
import gzip
#ftp://ftp.ensembl.org/pub/release-73/fasta/homo_sapiens/pep/README
def ensembl_ftp_fasta(release=74, specieslist=[]):
    dirlist = []
    if (specieslist==[]):
        print("Warning: No set of species selected")
        exit(0)
    ftp = FTP('ftp.ensembl.org')
    ftp.login()                    
    ftp.cwd('pub') #ftp://ftp.ensembl.org/pub/
    ftp.retrlines('LIST', callback=dirlist.append)           # list directory contents
    #print(dirlist)
    dirlist = [r for r in dirlist if "release-"+str(release)  in r and "fasta" in r]
    print(dirlist)
    if(len(dirlist) ==0):
        print("Nothing available for release-"+str(release))
        exit(1)
    ftprelhome = dirlist[0].split(" ")[-1]
    ftp.cwd(ftprelhome)
    ftp.retrlines('LIST', callback=dirlist.append)
    for s in specieslist:
        tmp = [d for d in dirlist if s in d][0]
        spec = tmp.split(" ")[-1]
        ftp.cwd(spec)
        ftp.cwd('pep')        
        falist = []
        ftp.retrlines('LIST', callback=falist.append)
        falist = [f for f in falist if "all.fa" in f][0]
        fafile = falist.split(" ")[-1]
        outfaname = spec+".fa.gz"
        outfa = open(outfaname,'wb')
        ftp.retrbinary("RETR "+fafile,outfa.write)
        outfa.close()
        xtract(outfaname)
        ftp.cwd("/pub/"+ftprelhome)
   
    ftp.quit()



def xtract(cfile, outpath = "."):
    
    if cfile.endswith('.gz'):
        dzf = ".".join(cfile.split(".")[:-1])
        #gzip 
        gzf = gzip.open(cfile,'rb')
        content = gzf.read()
        out = open(outpath+os.sep+dzf,'wb')
        out.write(content)
        out.close
        gzf.close()
    else: 
        print("Can't extract "+cfile)
    

def main():
    ensembl_ftp_fasta(release=74, specieslist=["homo_sapiens","mus_musculus","pan_troglodytes"])

main()