import tkinter as tk
from tkinter import filedialog
import Scythe_gff2loc 
import Scythe_ensembl2loc 

import os
import random
import string
class ScytheConvertDialogProteinortho():
    def __init__(self):
        root = tk.Toplevel()
        root.title("Convert Proteinortho output to .grp file")
        

    def callScythe_proteinOrtho2grp(self, infile, locfile, outfile):
        pass
        #read_proteinortho2grp(infile, locfile, outfile)

class ScytheConvertDialogLoc():
    def __init__(self):
        self.root = tk.Toplevel()
        self.root.title("Convert gff/tab-separated files to .loc files")
        
        #convert whole folder
        self.lab_dir = tk.Label(self.root,text="gff/tsv directory")
        self.st_dir = tk.StringVar()
        self.st_dir.set("")
        self.ent_dir = tk.Entry(self.root, width = 30, 
                                     textvariable = self.st_dir, state = tk.NORMAL)
        
        #convert one file
        self.lab_file = tk.Label(self.root,text="gff/tsv file")
        self.st_file = tk.StringVar()
        self.st_file.set("")
        self.ent_file = tk.Entry(self.root, width = 30, 
                                     textvariable = self.st_file, state = tk.NORMAL)
     
        
        #Buttons
        self.b_dir = tk.Button(self.root,text="open...",command=self.askDir, state=tk.NORMAL)#self.askopenfilename)
        self.b_file = tk.Button(self.root,text="open...",command=self.askFile, state=tk.NORMAL)#self.askopenfilename)
        self.b_ok = tk.Button(self.root,text="OK",command=self.onOK, state=tk.NORMAL)#self.askopenfilename)
        self.b_cancel = tk.Button(self.root,text="Cancel",command=self.onCancel, state=tk.NORMAL)#self.askopenfilename)


        self.lab_dir.grid(row=0, column=0)
        self.ent_dir.grid(row=0, column=1)
        self.b_dir.grid(row=0, column=2)
        self.lab_file.grid(row=1, column=0)
        self.ent_file.grid(row=1, column=1)
        self.b_file.grid(row=1, column=2)
        self.b_ok.grid(row=2, column=1,sticky=tk.E, padx=75)
        self.b_cancel.grid(row=2, column=1, sticky=tk.E)

        #self.b_convertOrtho.grid(row=3,column=3)
        
    def onOK(self):
        #dir
        if  (self.st_dir !=""):
            if self.st_dir.get():
                self.convertFolder(self.st_dir.get())
            if self.st_file.get():
                self.convertSingleFile(self.st_file.get())
        #or file
        
    def onCancel(self):
        pass
    def askDir(self):
        tmp= filedialog.askdirectory()
        print(tmp)
        self.ent_dir.config(state=tk.NORMAL)
        self.st_dir.set(tmp)


    def askFile(self):
        formats = [('GFF','*.gff'),('GFF','*.gff3'),('tab separated','*.tsv')]
        tmp= filedialog.askopenfilename(filetypes=formats)
        print(tmp)
        self.ent_file.config(state=tk.NORMAL)
        #self.st_file.set(tmp)
        
    def callConvertGff(self,file):  
        Scythe_gff2loc.read_gff2loc(file, file+".loc")

        
    def callConvertTsv(self,file):  
        #pass
        print("callConvertTsv",file)
        Scythe_ensembl2loc.readEnsemblLoc(file, file+".loc")
    def mkdir(self,path):
        
        
        rand = ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(2))
        #os.urandom(2).decode('utf-8')
        print(rand)
        if not os.path.isdir(path):
            os.makedirs(path)
        else:
            self.mkdir(path+rand)
            return(path+rand)
        return(path)
    def concatFiles(self,listoffiles,outputfilename):
        res =""
        for f in listoffiles:
            for l in open(f,'r'):
                res+=l
        outhandle = open(outputfilename,'w')
        outhandle.write(res)
    def convertFolder(self, folder):
        print("convertFolder")
        import os
        gfffiles = []
        tsvfiles = []
        for f in os.listdir(folder):
            if f.endswith(".gff"):
                print (f)
                gfffiles.append(f)
            elif f.endswith(".gff3"):
                print (f)
                gfffiles.append(f)
            elif f.endswith(".tsv"):
                print (f)
                tsvfiles.append(f)
        print(tsvfiles, gfffiles)
        for g in gfffiles:
            self.callConvertGff(folder+os.sep+g)
        for t in tsvfiles:
            print("tsvfiles", t)
            if t:
                self.callConvertTsv(folder+os.sep+t)
        #mkdir 
        locpath = self.mkdir(folder+os.sep+"loc")
        locfiles = []
        for f in os.listdir(folder):
            if f.endswith(".loc"):
                os.rename(folder+os.sep+f, locpath+os.sep+f)
                locfiles.append(locpath+os.sep+f)
        self.concatFiles(locfiles, folder+os.sep+"all.loc")
    def convertSingleFile(self, f):
        print("convertFile")
        import os
        gfffiles = []
        tsvfiles = []
        #for f in os.listdir(folder):
        if f.endswith(".gff"):
            print (f)
            self.callConvertGff(f)
                #gfffiles.append(f)
        elif f.endswith(".gff3"):
            print (f)
            self.callConvertGff(f)
            #    gfffiles.append(f)
        elif f.endswith(".tsv"):
            print (f)
            self.callConvertTsv(f)
            #    tsvfiles.append(f)
        #print(tsvfiles, gfffiles)
        #for f in locfiles:    
        #    os.rename(folder+os.sep+f, locpath+os.sep+f)
        