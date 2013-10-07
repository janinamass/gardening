import tkinter as tk
import configparser
from tkinter import filedialog
from tkinter import messagebox
from tkinter import ttk
from tkinter import OptionMenu
from tkinter import Scale
root=tk.Tk()
root.title("Scythe GUI alpha")
root.iconbitmap('@scy.xbm')
import os


global LOGSTR
LOGSTR = ""
global CURRENTCONFIG
CURRENTCONFIG = configparser.ConfigParser()
BACKUPCONFIG = configparser.ConfigParser()
print(CURRENTCONFIG,BACKUPCONFIG)
############ labels ##################
CF_MODE = "Mode"
CF_MODE_use_ensembl = "use_ensembl_api"
CF_MODE_use_local_files = "use_local_files"
CF_PATHS = "Paths"
CF_PATHS_fasta_directory = "fasta_directory"
CF_PATHS_loc_directory = "loc_directory"
CF_PATHS_grp_file = "grp_file"
CF_PATHS_output_directory = "output_directory"
CF_OUTPUT ="Output"
CF_OUTPUT_attach_output_prefix="attach_output_prefix"
CF_OUTPUT_output_prefix="output_prefix"
CF_OUTPUT_output_orthogroups = "output_orthogroups"
CF_OUTPUT_output_species_fasta = "output_species_fasta"
CF_OUTPUT_merge_species_fasta_with_defaults = "merge_species_fasta_with_defaults"
CF_CLEANUP = "Cleanup"
CF_CLEANUP_clean_up_directories = "clean_up_directories"
CF_RUN="Run_options"
CF_RUN_max_threads ="max_threads"
CF_RUN_split_input="split_input"
CF_RUN_split_each="split_each"
CF_RUN_use_seqan="use_seqan"
CF_PENALTIES = "Penalties"
CF_PENALTIES_gap_open_cost = "gap_open_cost"
CF_PENALTIES_gap_extend_cost="gap_extend_cost"
CF_PENALTIES_substitution_matrix="substitution_matrix"
CF_ALGORITHM = "Algorithm"
CF_ALGORITHM_use_global_max ="use_global_max"
CF_ALGORITHM_use_default="use_default"
CF_ALGORITHM_use_global_sum="use_global_sum"
CF_PARALOGS="Paralogs"
CF_PARALOGS_include_paralogs = "include_paralogs"
##################options###########################
OPTIONS = {}
#dropdown menus
yn =["yes","no"]
for o in [CF_PARALOGS_include_paralogs, CF_ALGORITHM_use_global_max,CF_ALGORITHM_use_default,
          CF_ALGORITHM_use_global_sum, CF_RUN_use_seqan, CF_RUN_split_input, CF_CLEANUP_clean_up_directories,
          CF_OUTPUT_merge_species_fasta_with_defaults, CF_OUTPUT_output_species_fasta, CF_OUTPUT_output_orthogroups,
          CF_OUTPUT_attach_output_prefix,CF_MODE_use_local_files,CF_MODE_use_ensembl]:
    OPTIONS[o]=yn

#PENALTIES:int entry
#PENALTIES:str entry



###################################################
SECTIONS = [CF_MODE,CF_PATHS, CF_OUTPUT, CF_CLEANUP, CF_RUN, CF_PENALTIES,CF_ALGORITHM,CF_PARALOGS]
MAXCONFIG = configparser.ConfigParser()
for i in SECTIONS:
    MAXCONFIG.add_section(i)
for i in [CF_MODE_use_ensembl,CF_MODE_use_local_files]:
    MAXCONFIG.set(CF_MODE,i,"unset")
for i in [CF_PATHS_fasta_directory,CF_PATHS_loc_directory,CF_PATHS_grp_file,CF_PATHS_output_directory ]:
    MAXCONFIG.set(CF_PATHS,i,"unset")
for i in [CF_OUTPUT_attach_output_prefix,CF_OUTPUT_output_prefix,CF_OUTPUT_output_orthogroups,CF_OUTPUT_output_species_fasta,CF_OUTPUT_merge_species_fasta_with_defaults ]:
    MAXCONFIG.set(CF_OUTPUT,i,"unset")
for i in [CF_CLEANUP_clean_up_directories]:
    MAXCONFIG.set(CF_CLEANUP,i,"unset")
for i in [CF_RUN_max_threads,CF_RUN_split_input,CF_RUN_split_each,CF_RUN_use_seqan]:
    MAXCONFIG.set(CF_RUN,i,"unset")
for i in [CF_PENALTIES_gap_open_cost,CF_PENALTIES_gap_extend_cost,CF_PENALTIES_substitution_matrix]:
    MAXCONFIG.set(CF_PENALTIES,i,"unset")
for i in [CF_ALGORITHM_use_global_max,CF_ALGORITHM_use_default,CF_ALGORITHM_use_global_sum   ]:
    MAXCONFIG.set(CF_ALGORITHM,i,"unset")
for i in [CF_PARALOGS_include_paralogs  ]:
    MAXCONFIG.set(CF_PARALOGS,i,"unset")
for i in MAXCONFIG.items():
    print(i)
######################################
#q'n'd
def initConfCurrent():
    global CURRENTCONFIG
    global MAXCONFIG
    for i in MAXCONFIG.sections():
        CURRENTCONFIG.add_section(i)
        for j in MAXCONFIG.options(i):
            print(i,j)
            CURRENTCONFIG.set(i,j,MAXCONFIG.get(i,j))
def backupConf():
    global CURRENTCONFIG
    global BACKUPCONFIG
    for i in CURRENTCONFIG.sections():
        try:
            BACKUPCONFIG.add_section(i)
        except configparser.DuplicateSectionError as e:
             pass
        for j in CURRENTCONFIG.options(i):
            BACKUPCONFIG.set(i,j,CURRENTCONFIG.get(i,j) )
def backupConfTo(newconf):
    global CURRENTCONFIG
    for i in CURRENTCONFIG.sections():
        try:
            newconf.add_section(i)
        except configparser.DuplicateSectionError as e:
             pass
        for j in CURRENTCONFIG.options(i):
            newconf.set(i,j,CURRENTCONFIG.get(i,j) )
def restoreConf():
     global BACKUPCONFIG
     global CURRENTCONFIG
     for i in BACKUPCONFIG.sections():
         try:
             CURRENTCONFIG.add_section(i)
         except configparser.DuplicateSectionError as e:
             pass
         for j in  BACKUPCONFIG.options(i):
             CURRENTCONFIG.set(i,j,BACKUPCONFIG.get(i,j))
######################################
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

class ConfigHandler():
    def __init__(self):
        self._currentconfig = configparser.ConfigParser()
    @property
    def currentconfig(self):
        return self._currentconfig
    @currentconfig.setter
    def currentconfig(self,config):
        self._currentconfig=config
        
    def reset(self):
        initConfCurrent()
        print("full reset")
    
class ScytheConvertDialog():
    pass

class Infobox():
    @logged
    def todo(self):
        message="Soon (tm)."
        messagebox.showinfo(title="Todo...", message = message )
    @logged
    def about(self):
        message="Scythe GUI 0.4 \nOctober 2013\nJ. Mass\n"
        messagebox.showinfo(title="About Scythe", message = message )
    
    def wanttosave(self):
        pass
    def bepatient(self):
        message="Scythe is running\n This may take some time.\n"
        messagebox.showinfo(title="Running", message = message )
    def saveConfig(self):
        formats = [('Scythe configuration','*.scy')]
        tmp= tk.filedialog.asksaveasfilename(parent=self.parent,filetypes=formats ,title="Save configuration as...")
    def showConfig(self):
        global CURRENTCONFIG
        tmp = ConfigHandler()
        print(CURRENTCONFIG)
        tmp.currentconfig=CURRENTCONFIG
        #read configuration
        print("tmp", tmp)
        print("tmp.cc", tmp.currentconfig)
        message = ""
        for section in tmp.currentconfig.sections():
            message += "["+section +"]\n"
            for option in tmp.currentconfig.options(section):
                message += " "+ option+ "="+ tmp.currentconfig.get(section, option)+"\n"
        #messagebox.showinfo(title="Config", message = "" )
        top = tk.Toplevel(root)
        top.title("Configuration")
        txt = tk.Text(top)
        scrollv = tk.Scrollbar(top, command=txt.yview)
        txt.insert(tk.INSERT,message)
        txt.configure(yscrollcommand=scrollv.set, state=tk.DISABLED, background="black", foreground="green" )
        
        txt.grid(row=0, column=0)
        scrollv.grid(row=0, column=1)
       

        
        
        
        
        
class ScytheMenu(tk.Frame):
  
    def __init__(self, parent):
        tk.Frame.__init__(self, parent)   
        self.parent = parent        
        self.initGUI()
        self.confighandler = ConfigHandler()
        initConfCurrent()
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

        optionsMenu = tk.Menu(menubar)
        optionsMenu.add_command(label="Show configuration...", command=self.onShowOptions)
        optionsMenu.add_command(label="Set configuration...", command=self.onSetOptions)



        infoMenu = tk.Menu(menubar)
        infoMenu.add_command(label="Show log...", command=self.onShowLog)
        
        helpMenu = tk.Menu(menubar)
        helpMenu.add_command(label="About...", command=self.onAbout)
        
        menubar.add_cascade(label="File", menu=fileMenu)
        menubar.add_cascade(label="Options", menu=optionsMenu)
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
        self.confighandler.reset()
    def onLoadConfig(self):
        formats = [('Scythe configuration','*.scy')]  
        tmp = tk.filedialog.askopenfilename(parent=self.parent,filetypes=[('Scythe configuration','*.scy')],title="Load configuration...")
        cfg = self.confighandler.currentconfig.read(tmp)
        global CURRENTCONFIG
        CURRENTCONFIG=self.confighandler.currentconfig
        if self.confighandler.currentconfig.get(CF_MODE,'use_local_files')=='yes':
            self.scythewizard.st_fastaDir.set(self.confighandler.currentconfig.get(CF_PATHS,'fasta_directory') )
            self.scythewizard.st_locDir.set(self.confighandler.currentconfig.get(CF_PATHS,'loc_directory') )
            self.scythewizard.st_grpFile.set(self.confighandler.currentconfig.get(CF_PATHS,'grp_file') )
            self.scythewizard.st_outDir.set(self.confighandler.currentconfig.get(CF_PATHS,'output_directory') )
            self.scythewizard.cb_use_local.select()
            self.scythewizard.cb_use_local.configure(state=tk.NORMAL)
            self.scythewizard.cb_use_ensembl.configure(state=tk.DISABLED)
            self.scythewizard.ent_fastaDir.configure(state=tk.NORMAL)
            self.scythewizard.ent_locDir.configure(state=tk.NORMAL)
            self.scythewizard.ent_grpFile.configure(state=tk.NORMAL)
            self.scythewizard.ent_outDir.configure(state=tk.NORMAL)
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
        print(tmp)
        global CURRENTCONFIG
        print(CURRENTCONFIG)
        try:
            out = open(tmp,"w")
            CURRENTCONFIG.write(out)
            out.close()
        except Error as e:
            print (e)
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
    def onShowOptions(self):
        Infobox().showConfig()
    def onSetConfigApply(self):
        print("configapply")
        Infobox().todo()
    def onSetConfigOK(self):
        #Infobox().todo()
        backupConf()
        print("BACKED UP")
        restoreConf()
        print("RESTORED-->CURRENTCONF set")
    def onSetConfigCancel(self):
        #self.restoreOldConfig()
        print("Config CANCEL")
        
        #Infobox().todo()
    #def restoreOldConfig(self):
    #    CURRENTCONFIG = configparser.ConfigParser(BACKUPCONFIG)
    #def saveNewConfig(self):
    #    CURRENTCONFIG = configparser.ConfigParser(self.tempconfig)
    #    BACKUPCONFIG = configparser.ConfigParser(CURRENTCONFIG)
    
    def applyPenalties(self):
        pass
    def applyOutput(self):
        pass
    def applyCleanup(self):
        pass
    def applyRun(self):
        pass
    def applyAlgorithm(self):
        pass
    def applyParalogs(self):
        pass
    def onSetOptions(self):
        global CURRENTCONFIG
        global MAXCONFIG
        global CF_MODE
        tmpconfig= configparser.ConfigParser()
        backupConfTo(tmpconfig)
        
        top = tk.Toplevel()
        top.title("Set configuration")
        nb = ttk.Notebook(top)
        b_config_ok = tk.Button(top, text="OK", command=top.destroy)
        b_config_ok.bind('<ButtonRelease-1>',self.onSetConfigOK())
        b_config_apply = tk.Button(top, text="Apply", command=self.onSetConfigApply)
        b_config_cancel = tk.Button(top, text="Cancel", command=top.destroy)
        b_config_cancel.bind('<ButtonRelease-1>',self.onSetConfigCancel())
        
        
        fr_paths = tk.Frame(nb,width=200, height=100)
        fr_penalties = tk.Frame(nb,width=200, height=100)
        fr_mode = ttk.Frame(nb,width=200, height=100)
        fr_output = ttk.Frame(nb,width=200, height=100)
        fr_cleanup = ttk.Frame(nb,width=200, height=100)
        fr_run = ttk.Frame(nb,width=200, height=100)
        fr_algorithm = ttk.Frame(nb,width=200, height=100)
        fr_paralogs = ttk.Frame(nb,width=200, height=100)
        
        #######labels########################
        txt_sec=[]
        txt_subsec={}
        for section in MAXCONFIG.sections():
            print( "["+section +"]\n")
            txt_sec.append(section)
            for opt in MAXCONFIG.options(section):
                try:
                    txt_subsec[section].append(opt)
                except KeyError as e:
                    txt_subsec[section]=[opt]
        lab_sec=[]
        lab_subsec={}
        dd_subsec={}
        var_subsec={}
        for t in txt_sec:
            lab_sec.append(tk.Label(fr_paths,text = t)) 
        for t in txt_subsec:
            print(t,txt_subsec[t])
            for u in txt_subsec[t]:
                if t == CF_MODE:
                    fr = fr_mode
                elif t == CF_PATHS:
                    fr = fr_paths
                elif t == CF_OUTPUT:
                    fr = fr_output
                elif t == CF_CLEANUP:
                    fr = fr_cleanup
                elif t == CF_RUN:
                    fr = fr_run
                elif t == CF_PENALTIES:
                    fr = fr_penalties
                elif t == CF_ALGORITHM:
                    fr = fr_algorithm
                elif t == CF_PARALOGS:
                    fr = fr_paralogs
                else:
                    print("No such section:",t)
                try:
                    lab_subsec[t].append(tk.Label(fr,text = u))
                    var_subsec[t].append(tk.StringVar(fr))
                    if u in OPTIONS:
                        dd_subsec[t].append(OptionMenu(fr,var_subsec[t][-1],*OPTIONS[u]))
                    else:
                        dd_subsec[t].append("")
                except KeyError as e:
                    try:
                        lab_subsec[t]=[tk.Label(fr,text = u)]
                        var_subsec[t]=[tk.StringVar(fr)]
                        if u in OPTIONS:
                            dd_subsec[t] = [OptionMenu(fr,var_subsec[t][-1],*OPTIONS[u])]
                        else:
                            dd_subsec[t] = [""]
                    except KeyError as e:
                        print(e)
                        dd_subsec[t].append("")
                        
        for t in lab_subsec:
            r=0
            c=0
            for i in  lab_subsec[t]:
                print(i.cget("text"))
                i.grid(row=r,column=c, sticky=tk.E)
                r+=1
                print(r,i.cget("text"))
        for t in dd_subsec:
            c=1
            r=0
            for i in dd_subsec[t]:
                print(i)
                if i is not "":
                    i.grid(row=r,column=c,sticky=tk.N)
                r+=1
                print(r)
        ######################################
        st_submat = tk.StringVar()
        sc_config_numthreads = Scale(fr_run, from_=1, to=8, orient=tk.HORIZONTAL)
        sc_config_numthreads.grid(row=0, column=1, sticky=tk.E)
        en_config_gapopen=tk.Entry(fr_penalties)
        en_config_gapextend=tk.Entry(fr_penalties)
        om_config_submat=tk.OptionMenu(fr_penalties,st_submat, *["EBLOSUM62","EDNAFULL"])
        om_config_submat.grid(row=2,column=1 )
        en_config_outpref=tk.Entry(fr_output, width=6)
        en_config_gapopen.grid(row=0, column=1)
        en_config_gapextend.grid(row=1, column=1)
        #en_config_submat.grid(row=2, column=1)
        en_config_outpref.grid(row=1, column=1)
        #nb.add(fr_mode, text=CF_MODE)
        #nb.add(fr_paths, text=CF_PATHS)
        nb.add(fr_penalties, text=CF_PENALTIES)
        nb.add(fr_output, text=CF_OUTPUT)
        nb.add(fr_cleanup, text=CF_CLEANUP)
        nb.add(fr_run, text=CF_RUN)
        nb.add(fr_algorithm, text=CF_ALGORITHM)
        nb.add(fr_paralogs, text=CF_PARALOGS)

        nb.grid()
        b_config_cancel.grid(row=1, column=0, sticky=tk.E,padx=115)
        b_config_apply.grid(row=1, column=0, sticky=tk.E,padx=50)
        b_config_ok.grid(row=1, column=0, sticky=tk.E)
       
    
class ScytheWizard(tk.Tk):
    def __init__(self, parent):
        self.parent = parent     
        self.initWizard()
        
        
    def initWizard(self):
       
        #Labels
        self.lab_fastaDir = tk.Label(text="Fasta Directory")
        self.lab_locDir = tk.Label(text=".loc Directory")
        self.lab_grpFile = tk.Label(text=".grp File")
        self.lab_outDir = tk.Label(text="Output Directory")

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
        self.st_outDir= tk.StringVar()
        self.st_outDir.set("")
        #Entries
        self.ent_fastaDir = tk.Entry(root, width = 30, 
                                     textvariable = self.st_fastaDir, state = tk.DISABLED)
        self.ent_locDir = tk.Entry(root, width = 30, 
                                     textvariable = self.st_locDir, state = tk.DISABLED)
        self.ent_grpFile = tk.Entry(root, width = 30, 
                                     textvariable = self.st_grpFile, state = tk.DISABLED)
        self.ent_outDir = tk.Entry(root, width = 30, 
                                     textvariable = self.st_outDir, state = tk.DISABLED)
        #Buttons
        self.b_loadConfig = tk.Button()
        self.b_saveConfig = tk.Button()
        self.b_fastaDir = tk.Button(text="open...",command=self.askFastaDir, state=tk.DISABLED)#self.askopenfilename)
        self.b_locDir = tk.Button(text="open...",command=self.askLocDir, state=tk.DISABLED)#self.askopenfilename)
        self.b_grpFile = tk.Button(text="open...",command=self.askGrpFile, state=tk.DISABLED)#self.askopenfilename)
        self.b_outDir = tk.Button(text="open...",command=self.askOutDir, state=tk.DISABLED)#self.askopenfilename)

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
        
        self.lab_outDir.grid(row=4, column=0, sticky="E")
        self.ent_outDir.grid(row=4, column=1, sticky="E")
        self.b_outDir.grid(row=4, column=2, sticky="E")
        
        #self.b_convertOrtho.grid(row=3,column=3)
        
    def useLocal(self):
        global CURRENTCONFIG
        print("CURRENTCONFIG",CURRENTCONFIG)
        if self.int_local.get() ==1:
            CURRENTCONFIG.set(CF_MODE, "use_local_files", "yes")
            CURRENTCONFIG.set(CF_MODE, "use_ensembl_api", "no")
            #self.st_fastaDir
            self.ent_fastaDir.config(state=tk.NORMAL)
            self.ent_locDir.config(state=tk.NORMAL)
            self.ent_grpFile.config(state=tk.NORMAL)
            self.ent_outDir.config(state=tk.NORMAL)

            self.b_fastaDir.config(state=tk.NORMAL)
            self.b_locDir.config(state=tk.NORMAL)
            self.b_grpFile.config(state=tk.NORMAL)
            self.b_outDir.config(state=tk.NORMAL)

            self.cb_use_ensembl.config(state=tk.DISABLED)
            
        else:
            CURRENTCONFIG.set("Mode", "use_local_files", "no")
            CURRENTCONFIG.set("Mode", "use_ensembl_api", "no")
            #self.st_fastaDir="some"
            self.ent_fastaDir.config(state=tk.DISABLED)
            self.ent_locDir.config(state=tk.DISABLED)
            self.ent_grpFile.config(state=tk.DISABLED)
            self.ent_outDir.config(state=tk.DISABLED)

            self.b_fastaDir.config(state=tk.DISABLED)
            self.b_locDir.config(state=tk.DISABLED)
            self.b_grpFile.config(state=tk.DISABLED)
            self.b_outDir.config(state=tk.DISABLED)

            self.cb_use_ensembl.config(state=tk.NORMAL)

    def askFastaDir(self):
        tmp= filedialog.askdirectory()
        print(tmp)
        self.ent_fastaDir.config(state=tk.NORMAL)
        self.st_fastaDir.set(tmp)
        if self.ent_locDir.get()=="":
            if os.path.isdir(os.path.split(tmp)[0]+os.sep+"loc"):
                self.st_locDir.set(os.path.split(tmp)[0]+os.sep+"loc")
                CURRENTCONFIG.set("Paths", "loc_directory", os.path.split(tmp)[0]+os.sep+"loc")
        CURRENTCONFIG.set("Paths", "fasta_directory",tmp)
        return tmp
    def askLocDir(self):
        tmp= filedialog.askdirectory()
        print(tmp)
        self.ent_locDir.config(state=tk.NORMAL)
        self.st_locDir.set(tmp)
        if self.ent_locDir.get()=="":
            if os.path.isdir(os.path.split(tmp)[0]+os.sep+"fa"):
                self.st_locDir.set(os.path.split(tmp)[0]+os.sep+"fa")
                CURRENTCONFIG.set("Paths", "fasta_directory",os.path.split(tmp)[0]+os.sep+"fa")
        CURRENTCONFIG.set("Paths", "loc_directory",tmp)
        return filedialog.tmp()
    def askGrpFile(self):
        tmp= filedialog.askopenfilename()
        self.ent_grpFile.config(state=tk.NORMAL)
        self.st_grpFile.set(tmp)
        CURRENTCONFIG.set("Paths", "grp_file",tmp)
        return tmp
    def askOutDir(self):
        tmp= filedialog.askdirectory()
        print(tmp)
        self.ent_outDir.config(state=tk.NORMAL)
        self.st_outDir.set(tmp)
        CURRENTCONFIG.set("Paths", "out_dir",tmp)
        #if self.ent_locDir.get()=="":
        #    if os.path.isdir(os.path.split(tmp)[0]+os.sep+"fa"):
        #        self.st_Dir.set(os.path.split(tmp)[0]+os.sep+"fa")
        return filedialog.tmp()
    
    def useEnsembl(self):
        CURRENTCONFIG.set("Mode", "use_local_files", "no")
        CURRENTCONFIG.set("Mode", "use_ensembl_api", "yes")
        print(self.int_ensembl)
        print(self.ent_fastaDir)
        if self.int_ensembl.get() ==1:
            print("True")
            #self.st_fastaDir="none"
            print("True", self.st_fastaDir)
            self.ent_fastaDir.config(state=tk.DISABLED)
            self.ent_locDir.config(state=tk.DISABLED)
            self.ent_outDir.config(state=tk.NORMAL)

            self.b_fastaDir.config(state=tk.DISABLED)
            self.b_locDir.config(state=tk.DISABLED)
            self.b_outDir.config(state=tk.NORMAL)
            self.cb_use_local.config(state=tk.DISABLED)
        else:
            self.cb_use_local.config(state=tk.NORMAL)
            self.b_outDir.config(state=tk.DISABLED)
            self.ent_outDir.config(state=tk.DISABLED)
            CURRENTCONFIG.set("Mode", "use_local_files", "no")
            CURRENTCONFIG.set("Mode", "use_ensembl_api", "no")
            #self.cb_use_local.config(state=tk.NORMAL)

            #self.st_fastaDir="some"
            #print("False", self.st_fastaDir,self.ent_fastaDir )
            #self.ent_fastaDir.config(state=tk.NORMAL)
            #self.ent_locDir.config(state=tk.NORMAL)

            
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