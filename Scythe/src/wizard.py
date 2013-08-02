#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from tkinter import *
from tkinter.filedialog import askopenfilename, askdirectory
from datetime import datetime
import os
import threading
import time
import queue

#askopenfilename
clickedDct = {}
tk = Tk()
tk.title("Scythe GUI alpha")



class StdoutRedirector(Text):
    def flush(self):
        pass

    def write(self, str):
        self.insert(END,str)
class WizFrame(Frame):
    def setName(self, name):
        self.name = name
    def addPrev(self, prev=None):
        self.prev = prev
    def addNxt(self, nxt=None):
        self.nxt = nxt
    def printWiz(self):
        if self.prev:
            print(self.prev)
        if self.nxt:
            print(self.nxt)
        if self.name:
            print(self.name)
    def gotoPrev(self):
        self.pack_forget()
        self.prev.next = self
        self.prev.pack()
    def gotoNxt(self):
        if self.nxt is not None:
            self.pack_forget()
            self.nxt.prev = self
            self.nxt.pack()
        else:
            print("was none", self)

def goPrev(wiz):
    wiz.gotoPrev()
    print("goPrev",wiz,wiz.prev)
def goNext(wiz):
    wiz.gotoNxt()
    print("goNext",wiz,wiz.nxt)

def startFromScratch(wiz):
    global current
    page1 = WizFrame(tk)
    Button(master=page1,text='Quit', command = quitProg).pack(side=LEFT,fill=X)
    Button(master=page1,text='<- previous', command= startOver).pack(side=LEFT,fill=X)
    Button(master=page1,text='next ->').pack(fill=X)
    wiz.pack_forget()
    page1.pack()
    current = page1

def openDir(dir):
    global inDir
    inDir = askdirectory()
    print(inDir)
    dir.delete(0,END)
    dir.insert(0,inDir)
    return(inDir)
def openOutDir(dir):
    global outDir
    outDir = askdirectory()
    print(outDir)
    dir.delete(0,END)
    dir.insert(0,outDir)
    return(outDir)
def openFile(grp):
    global inGrp
    inGrp = askopenfilename()
    print(inGrp)
    grp.delete(0,END)
    grp.insert(0,inGrp)
    
    return(inGrp)

def chooseOutDir(wiz):
    global current
    global outDir
    outDir = os.getcwd()
    page1 = WizFrame(tk)
    page1.addPrev(wiz)
    wiz.addNxt(page1)
    odr=Entry(master=page1)
    if outDir is not None:
        odr.delete(0,END)
        odr.insert(0,outDir)
    
    odr.pack(side=TOP)
    Button(master=page1,text='Select Output Directory', command = lambda d=odr:  openOutDir(d)).pack(side=TOP)
    
    Button(master=page1,text='<- previous', command=lambda w = page1: goPrev(w)).pack(fill=X)
    Button(master=page1,text='next ->',command=lambda w =page1:chooseParam(w)).pack(fill=X)
    Button(master=page1,text='Quit', command = quitProg).pack(side=TOP,fill=X)
    wiz.pack_forget()
    page1.pack()
    current = page1



def startResume(wiz):
    global current
    global inDir
    global inGrp
    page1 = WizFrame(tk)
    page1.addPrev(wiz)
    wiz.addNxt(page1)
    
    dr=Entry(master=page1)
    if inDir is not None:
        dr.delete(0,END)
        dr.insert(0,inDir)
    
    dr.pack(side=TOP)
    Button(master=page1,text='Select Dir', command = lambda d=dr:  openDir(d)).pack(side=TOP)
        
    grp=Entry(master=page1)
    if inGrp is not None:
        grp.delete(0,END)
        grp.insert(0,inGrp)
    grp.pack(side=TOP)
    Button(master=page1,text='Select File', command = lambda g=grp: openFile(g)).pack(side=TOP)
    Button(master=page1,text='<- previous', command=lambda w = page1: goPrev(w)).pack(fill=X)
    Button(master=page1,text='next ->',command=lambda w =page1:chooseOutDir(w)).pack(fill=X )
    Button(master=page1,text='Quit', command = quitProg).pack(side=TOP,fill=X)
    wiz.pack_forget()
    page1.pack()
    current = page1
    #TODO: dirs can't be None

def chooseParam(wiz):
#TODO:check if folder contains fa and loc
#TODO: check if their names match
    global current
    global inDir
    global inGrp
    global outDir
    global gapOpen
    global gapExtend
    
    page1 = WizFrame(tk)
    page1.addPrev(wiz)
    wiz.addNxt(page1)
    goEntry=Entry(master=page1)
    goEntry.pack()
    geEntry=Entry(master=page1)
    geEntry.pack()
    if gapOpen is not None:
        goEntry.delete(0,END)
        goEntry.insert(0,gapOpen)
    if gapExtend is not None:
        geEntry.delete(0,END)
        geEntry.insert(0,gapExtend)
    Button(master=page1,text='<- previous', command=lambda w = page1: goPrev(w)).pack(fill=X)
    Button(master=page1,text='next ->',command=lambda w =page1,go=goEntry,ge=geEntry:runScythe(w,go,ge)).pack(fill=X )
    Button(master=page1,text='Quit', command = quitProg).pack(side=TOP,fill=X)
    wiz.pack_forget()
    page1.pack()
    current = page1
    pass
    #TODO: penalties can't be None, can't be letters





def runScythe(wiz,gapOpenEntry,gapExtendEntry):
    global current
    global inDir
    global inGrp
    global outDir
    global faFiles
    global locFiles
    global gapOpen
    global gapExtend
    global log_startTime
    gapOpen=gapOpenEntry.get()
    print(gapOpen)
    gapExtend=gapExtendEntry.get()
    print(gapExtend)
    
    page1 = WizFrame(tk)
    page1.addPrev(wiz)
    wiz.addNxt(page1)

    Label(master=page1,text='Here goes a summary').pack(fill=X)
    Button(master=page1,text='<- previous', command=lambda w = page1: goPrev(w)).pack(fill=X)
    Button(master=page1,text='run Scythe',command=lambda w =page1:beBusy(w)).pack(fill=X )
    Button(master=page1,text='Quit', command = quitProg).pack(side=TOP,fill=X)
    wiz.pack_forget()
    page1.pack()
    current = page1

def doSth(num):
    for i in range(1,num):
        print(i)
        time.sleep(1)
    
def beBusy(wiz):
    global current
    global inDir
    global inGrp
    global outDir
    global faFiles
    global locFiles
    global gapOpen
    global gapExtend
    global log_startTime
    print(gapOpen)
    print(gapExtend)
    textRedir = StdoutRedirector()
    sys.stdout= textRedir
    textRedir.pack()
    page1 = WizFrame(tk)
    page1.addPrev(wiz)
    wiz.addNxt(page1)
    Label(master=page1,text='Scythe running').pack(fill=X)
    Button(master=page1,text='Cancel', command = cancelProg).pack(side=TOP,fill=X)
    Button(master=page1,text='Start over', command = startOver).pack(side=TOP,fill=X)
    Button(master=page1,text='Quit', command = quitProg).pack(side=TOP,fill=X)

    wiz.pack_forget()
    page1.pack()
    current = page1
    log_startTime = str(datetime.now().time())
    t = threading.Thread(target=doSth, args=(6,))
    t.start()
    
    print(log_startTime)
    
def cancelProg():
    global inDir
    global inGrp
    global outDir
    global gapOpen
    global gapExtend
    log_cancelTime = str(datetime.now().time())
    s = ""
    s +="\n Run cancelled:\t"+log_cancelTime
    s +="\ninDir:\t"+inDir
    s +="\ninGrp:\t"+inGrp
    s +="\noutDir:\t"+outDir
    s +="\n gap Open:\t"+gapOpen
    s +="\n gap Extend:\t"+gapExtend
    print(s)

def startOver():
    sys.stdout=sys.__stdout__

    global current
    current.pack_forget()
    page0 = WizFrame(tk)
    page0.addPrev()
    page0.addNxt()
    page0.setName("pg0")
    page0.printWiz()
    Button(master=page0,text='Quit', command = quitProg).pack(side=LEFT,fill=X)
    Button(master=page0,text='Start from scratch', command = lambda w = page0: startFromScratch(w)).pack(side=LEFT,fill=X)
    Button(master=page0,text='Resume / start from prepared files', command = lambda w = page0: startResume(w)).pack(fill=X)

    page0.pack()
    current = page0

def quitProg():
    sys.exit(0)

inDir=None
inGrp=None
outDir=None
gapOpen=str(10.0)
gapExtend=str(5.0)
log_startTime=None

page0 = WizFrame(tk)
current = page0
startOver()
####################
# call runScythe() #
####################

#runScythe(groups=groups, delim=delim, 
#asID=asID, faFileList=faFileList, 
#namesList=namesList, cleanUp=cleanUp, 
#stopAfter=stopAfter, inDir=inDir, outDir=outDir, gapOpen=gapOpen, gapExtend=gapExtend)
#
#


tk.mainloop()