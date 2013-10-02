import tkinter as tk
import configparser
from tkinter import filedialog
from tkinter import messagebox

root=tk.Tk()
root.title("Scythe GUI alpha")
root.iconbitmap('@scy.xbm')
import os
#todo if ../loc has fa: fill in
# vice versa

global LOGSTR
LOGSTR = ""
global CURRENTCONFIG
CURRENTCONFIG = configparser.ConfigParser()
def logged(f):
    global LOGSTR
    def wrapped(*args, **kargs):
        global LOGSTR
        print ("%s called..." % f.__name__)
        try:
            LOGSTR=LOGSTR+f.__name__+str(args)+str(kargs)
            return f(*args, **kargs)
        finally:
            print ("..Done.")
            print(LOGSTR)
    return wrapped

class ConfigHandler( ):
    def __init__(self):
        self._currentconfig = configparser.ConfigParser()
    @property
    def currentconfig(self):
        return self._currentconfig
    @currentconfig.setter
    def currentconfig(self,config):
        self._currentconfig=config
    
        
    #    initConfigHandler()
    #def initConfigHandler(self):
    #    global CURRENTCONFIG
    
class ScytheConvertDialog():
    pass

class Infobox():
    @logged
    def todo(self):
        message="Soon (tm)."
        messagebox.showinfo(title="Todo...",message= message )
    @logged
    def about(self):
        message="Scythe GUI 0.4 \nOctober 2013\n"
        messagebox.showinfo(title="About Scythe",message= message )
        tmp="blabl"
        return(tmp)
    
    def wanttosave(self):
        pass
    def bepatient(self):
        message="Scythe is running\n This may take some time.\n"
        messagebox.showinfo(title="Running",message= message )
    def saveConfig(self):
        formats = [('Scythe configuration','*.scy')]
        tmp= tk.filedialog.asksaveasfilename(parent=self.parent,filetypes=formats ,title="Save configuration as...")
        
class ScytheMenu(tk.Frame):
  
    def __init__(self, parent):
        tk.Frame.__init__(self, parent)   
        self.parent = parent        
        self.initGUI()
        self.confighandler = ConfigHandler()
        self.scythewizard= ScytheWizard(self.parent)
    def initGUI(self):
        menubar = tk.Menu(self.parent)
        self.parent.config(menu=menubar)
        
        fileMenu = tk.Menu(menubar)
        fileMenu.add_command(label="New run...", command=self.onNewRun)
        #fileMenu.add_command(label="Convert files...", command=self.onConvertFiles)
        fileMenu.add_command(label="Load configuration...", command=self.onLoadConfig)
        fileMenu.add_command(label="Save configuration...", command=self.onSaveConfig)
        
        convertMenu = tk.Menu(fileMenu)
        convertMenu.add_command(label="convert orthology information to .grp", command=self.onConvertToGrp)
        convertMenu.add_command(label="convert loci/transcript information to .loc", command=self.onConvertToLoc)

        fileMenu.add_cascade(label='Convert files...', menu=convertMenu, underline=0)

        fileMenu.add_command(label="Exit", command=self.onExit)

        infoMenu = tk.Menu(menubar)
        infoMenu.add_command(label="Show log...", command=self.onShowLog)
        
        helpMenu = tk.Menu(menubar)
        helpMenu.add_command(label="About...", command=self.onAbout)
        
        menubar.add_cascade(label="File", menu=fileMenu)
        menubar.add_cascade(label="Info", menu=infoMenu)
        menubar.add_cascade(label="Help", menu=helpMenu)
        #self.onNewRun()
    def onConvertToGrp(self):
        pass
    def onConvertToLoc(self):
        pass
    def onExit(self):
        self.quit()
    def onNewRun(self):
        self.scythewizard=ScytheWizard(self.parent)
    def onLoadConfig(self):
        formats = [('Scythe configuration','*.scy')]  
        tmp = tk.filedialog.askopenfilename(parent=self.parent,filetypes=[('Scythe configuration','*.scy')],title="Load configuration...")
        cfg = self.confighandler.currentconfig.read(tmp)
        if self.confighandler.currentconfig.get('Mode','use_local_files')=='yes':
            self.scythewizard.st_fastaDir.set(self.confighandler.currentconfig.get('Local_directories','fasta_directory') )
            self.scythewizard.st_locDir.set(self.confighandler.currentconfig.get('Local_directories','loc_directory') )
            self.scythewizard.st_grpFile.set(self.confighandler.currentconfig.get('Local_directories','grp_file') )

            self.scythewizard.cb_use_local.select()
            self.scythewizard.cb_use_local.configure(state=tk.NORMAL)
            self.scythewizard.cb_use_ensembl.configure(state=tk.DISABLED)
            self.scythewizard.ent_fastaDir.configure(state=tk.ENABLED)
            self.scythewizard.ent_locDir.configure(state=tk.ENABLED)
            self.scythewizard.ent_grpFile.configure(state=tk.ENABLED)

        #for section in cfg.sections():
        #    print(section)
        #    for option in cfg.options(section):
        #        print (" ", option, "=", cfg.get(section, option))
        
       
        
        #self.scythewizard.st_fastaDir.set(cfg["Local_directories"]['fasta_directory'])
        #print(self.confighandler.cuurentconfig())
        return tmp
    def onSaveConfig(self):
        formats = [('Scythe configuration','*.scy')]  
        
        tmp= tk.filedialog.asksaveasfilename(parent=self.parent,filetypes=formats ,title="Save configuration as...")
        #self.ent_grpFile.config(state=tk.NORMAL)
        #self.st_grpFile.set(tmp)
        return tmp
    def onShowLog(self):
        pass
    def onEnsembl(self):
        pass
    
    def onLocal(self):
        pass
    def onAbout(self):
        Infobox().about()

    def onConvertFiles(self):
        pass







class ScytheWizard(tk.Tk):
    def __init__(self, parent):
        self.parent = parent     
        self.initWizard()
        
        
    def initWizard(self):
       
        #Labels
        self.lab_fastaDir = tk.Label(text="Fasta Directory")
        self.lab_locDir = tk.Label(text=".loc Directory")
        self.lab_grpFile = tk.Label(text=".grp File")
        #Ints
        self.int_ensembl = tk.IntVar()
        self.int_local = tk.IntVar()
        
        #Strings
        self.st_fastaDir = tk.StringVar()
        self.st_fastaDir.set("")
        self.st_locDir= tk.StringVar()
        self.st_locDir.set("")
        self.st_grpFile= tk.StringVar()
        self.st_grpFile.set("")
        #Entries
        self.ent_fastaDir = tk.Entry(root, width = 30, 
                                     textvariable = self.st_fastaDir, state = tk.DISABLED)
        self.ent_locDir = tk.Entry(root, width = 30, 
                                     textvariable = self.st_locDir, state = tk.DISABLED)
        self.ent_grpFile = tk.Entry(root, width = 30, 
                                     textvariable = self.st_grpFile, state = tk.DISABLED)
        #Buttons
        self.b_loadConfig = tk.Button()
        self.b_saveConfig = tk.Button()
        self.b_fastaDir = tk.Button(text="open...",command=self.askFastaDir, state=tk.DISABLED)#self.askopenfilename)
        self.b_locDir = tk.Button(text="open...",command=self.askLocDir, state=tk.DISABLED)#self.askopenfilename)
        self.b_grpFile = tk.Button(text="open...",command=self.askGrpFile, state=tk.DISABLED)#self.askopenfilename)
        
        #self.b_convertOrtho = tk.Button(text="convert file to .grp")
        #todo convert orthomcl, proteinortho, tabsep
        #self.b_convertLoci = tk.Button(text="convert files to .loc")
        #todo convert tabsep gff
        #self.b_saveConfig = tk.Button(text='Save Config...')
        #self.b_loadConfig = tk.Button(text='Load Config...')
        #Button(master, text='Load Config...', command=None).grid(row=3, column=4, sticky=W, pady=0)
        
        #ChechButtons
        self.cb_use_ensembl = tk.Checkbutton(root, text='use ENSEMBL',variable=self.int_ensembl,
                                         command=self.useEnsembl)
        self.cb_use_local = tk.Checkbutton(root, text='use local files',variable=self.int_local,
                                         command=self.useLocal, state=tk.NORMAL)

       
        
   #     self.ent_fastaDir = tk.Entry(root, width = 100, 
   #                                  textvariable = self.st_fastaDir, state = tk.DISABLED).grid(row=3, column=0, sticky=W, pady=4)
        #add to grid
        self.cb_use_ensembl.grid(row=0, column=0 )
        self.cb_use_local.grid(row=0, column=1)
        
        self.lab_fastaDir.grid(row=1, column=0)
        self.ent_fastaDir.grid(row=1, column=1)
        self.b_fastaDir.grid(row=1, column=2)
        
        self.lab_locDir.grid(row=2, column=0)
        self.ent_locDir.grid(row=2, column=1)
        self.b_locDir.grid(row=2, column=2)
        
        
        self.lab_grpFile.grid(row=3, column=0)
        self.ent_grpFile.grid(row=3, column=1)
        self.b_grpFile.grid(row=3, column=2)
        #self.b_convertOrtho.grid(row=3,column=3)
        
    def useLocal(self):
        if self.int_local.get() ==1:
            #self.st_fastaDir
            self.ent_fastaDir.config(state=tk.NORMAL)
            self.ent_locDir.config(state=tk.NORMAL)
            self.ent_grpFile.config(state=tk.NORMAL)
            self.b_fastaDir.config(state=tk.NORMAL)
            self.b_locDir.config(state=tk.NORMAL)
            self.b_grpFile.config(state=tk.NORMAL)
            
            self.cb_use_ensembl.config(state=tk.DISABLED)
            
        else:
            #self.st_fastaDir="some"
            self.ent_fastaDir.config(state=tk.DISABLED)
            self.ent_locDir.config(state=tk.DISABLED)
            self.ent_grpFile.config(state=tk.DISABLED)
            self.b_fastaDir.config(state=tk.DISABLED)
            self.b_locDir.config(state=tk.DISABLED)
            self.b_grpFile.config(state=tk.DISABLED)
            self.cb_use_ensembl.config(state=tk.NORMAL)

    def askFastaDir(self):
        tmp= filedialog.askdirectory()
        print(tmp)
        self.ent_fastaDir.config(state=tk.NORMAL)
        self.st_fastaDir.set(tmp)
        if self.ent_locDir.get()=="":
            if os.path.isdir(os.path.split(tmp)[0]+os.sep+"loc"):
                self.st_locDir.set(os.path.split(tmp)[0]+os.sep+"loc")
        return tmp
    def askLocDir(self):
        tmp= filedialog.askdirectory()
        print(tmp)
        self.ent_locDir.config(state=tk.NORMAL)
        self.st_locDir.set(tmp)
        if self.ent_locDir.get()=="":
            if os.path.isdir(os.path.split(tmp)[0]+os.sep+"fa"):
                self.st_locDir.set(os.path.split(tmp)[0]+os.sep+"fa")
        return filedialog.tmp()
    def askGrpFile(self):
        tmp= filedialog.askopenfilename()
        self.ent_grpFile.config(state=tk.NORMAL)
        self.st_grpFile.set(tmp)
        return tmp
    
    def useEnsembl(self):
        print(self.int_ensembl)
        print(self.ent_fastaDir)
        if self.int_ensembl.get() ==1:
            print("True")
            #self.st_fastaDir="none"
            print("True", self.st_fastaDir)
            self.ent_fastaDir.config(state=tk.DISABLED)
            self.ent_locDir.config(state=tk.DISABLED)
            self.b_fastaDir.config(state=tk.DISABLED)
            self.b_locDir.config(state=tk.DISABLED)
            self.cb_use_local.config(state=tk.DISABLED)
        else:
            #self.st_fastaDir="some"
            print("False", self.st_fastaDir,self.ent_fastaDir )
            #self.ent_fastaDir.config(state=tk.NORMAL)
            #self.ent_locDir.config(state=tk.NORMAL)

            self.cb_use_local.config(state=tk.NORMAL)
   # def useLocal(self):
   #     if self.b_ensembl == 1:
   #         self.ent.configure(state='disabled')
   #     else:
   #         self.ent.configure(state='normal')       

#app=ScytheWizard()
app=ScytheMenu(root)
#ap=App()
root.mainloop()

#config = configparser.ConfigParser()
#config.sections()
#config.read('sampleconfig.scy')
#print(config.sections())
#global logger
#logger = ""

#master = Tk()
#master.iconbitmap('@scy.xbm')
#master.title("Scythe Wizard alpha")
#master.geometry("300x300")


#cb_use_ensembl= Checkbutton(master, text="Use Ensembl", variable=None).grid(row=0, column=1, sticky=W, pady=0)
#cb_use_local = Checkbutton(master, text="Use local files", variable=None).grid(row=0, column=2, sticky=W, pady=0)

#Label(master, text="Fasta Directory").grid(row=1,sticky=W, pady=0)
#Label(master, text="Loc Directory").grid(row=2,sticky=W, pady=0)

#e1 = Entry(master)
#e2 = Entry(master)

#e1.grid(row=1, column=1)
#e2.grid(row=2, column=1)
#Button(master, text='Quit', command=master.quit).grid(row=3, column=0, sticky=W, pady=4)
#Button(master, text='Show', command=None).grid(row=3, column=1, sticky=W, pady=4)
#Button(master, text='Save Config...', command=None).grid(row=3, column=3, sticky=W, pady=0)
#Button(master, text='Load Config...', command=None).grid(row=3, column=4, sticky=W, pady=0)
#Button(master, text='Dummy2', command=None).grid(row=3, column=5, sticky=W, pady=0)

#master.mainloop()
