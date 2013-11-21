from ftplib import FTP
import os
import tarfile
import zipfile
import gzip
#ftp://ftp.ensembl.org/pub/release-73/fasta/homo_sapiens/pep/README
def ensembl_ftp_fasta(release=73, specieslist=["homo_sapiens"]):
    dirlist = []
    ftp = FTP('ftp.ensembl.org')
    ftp.login()                    
    ftp.cwd('pub') #ftp://ftp.ensembl.org/pub/
    ftp.retrlines('LIST', callback=dirlist.append)           # list directory contents
    dirlist = [r for r in dirlist if "release-73"  in r and "fasta" in r]
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
        #print(tmp.split(" ")[-1])
        
        
    #ftp.retrbinary('RETR README', open('README', 'wb').write)
    ftp.quit()



def xtract(cfile, outpath = "."):
    #xtr = None
    #mode = None
    #if cfile.endswith('.tar.gz') or cfile.endswith('.tgz'):
    #    xtr, mode = tarfile.open, 'r:gz'
    #elif cfile.endswith('.tar.bz2') or cfile.endswith('.tbz'):
    #    xtr, mode = tarfile.open, 'r:bz2'
    #elif cfile.endswith(".zip") or cfile.endswith('.gz'):
        #xtr, mode = zipfile.ZipFile, 'r'
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
    
    #cwd = os.getcwd()
    #os.chdir(outpath)
    
    #try:
    #    xfile = xtr(cfile, mode)
    #    try: 
    ##        xfile.extractall()
     #   finally:
     #       xfile.close()
    #finally:
    #    os.chdir(cwd)

def main():
    ensembl_ftp_fasta(release=73, specieslist=["homo_sapiens","mus_musculus"])

main()