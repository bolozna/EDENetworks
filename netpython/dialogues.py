"""
EDENetworks, a genetic network analyzer
Copyright (C) 2011  Aalto University

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
"""
"""

Collection of dialogue windows


"""

import pynet,os,netio,netext,visuals,eden,transforms
import random
import heapq
import string
import percolator
import shutil
from math import ceil
from Tkinter import *
import tkMessageBox
#from pylab import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2TkAgg

# NEW DIALOGUE WINDOWS / JS / MAY-JUNE 09


class MySimpleDialog(Toplevel):
    '''Master class for a dialog popup window.
       Functions body() and apply() to be overridden
       with whatever the dialog should be doing.'''

    def __init__(self,parent,title=None):

        Toplevel.__init__(self,parent)
        self.transient(parent)


        self.title(title)

        self.parent=parent
        self.result=None

        body=Frame(self)
        self.initial_focus=self.body(self,body)
        body.pack(padx=5,pady=5)

        self.buttonbox()
        self.grab_set()

        if not self.initial_focus:
            self.initial_focus(self)

        self.protocol("WM_DELETE_WINDOW",self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.wait_window(self)

    def body(self,masterclass,masterwindow):
        pass

    def buttonbox(self):
        """OK and Cancel buttons"""

        box=Frame(self)
        w=Button(box,text="OK",width=10,command=self.ok,default=ACTIVE)
        w.pack(side=LEFT,padx=5,pady=5)
        w=Button(box,text="Cancel",width=10,command=self.cancel)
        w.pack(side=LEFT,padx=5,pady=5)

        self.bind("&lt;Return>",self.ok)
        self.bind("&lt;Escape",self.cancel)

        box.pack()

    def ok(self,event=None):

        if not self.validate():
            self.initial_focus.focus_set()
            return

        self.withdraw()
        self.update_idletasks()
        self.applyme()
        self.cancel()

    def cancel(self,event=None):

        self.parent.focus_set()
        self.destroy()

    def validate(self):
        return 1

    def applyme(self):
        pass

    def displayBusyCursor(self):

        self.parent.configure(cursor='watch')
        self.parent.update()
        self.parent.after_idle(self.removeBusyCursor)

    def removeBusyCursor(self):
        self.parent.configure(cursor='arrow')
        

class WLogbinDialog(MySimpleDialog):
    """Asks for the number of bins for log binning
       and allows linear bins for 1...10"""


    def __init__(self,parent,title=None):

        Toplevel.__init__(self,parent)
        self.configure(bg='Gray80')
        self.transient(parent)

        if title:
            self.title=title

        self.parent=parent
        self.result=None
        self.linfirst=IntVar()
        self.numbins=StringVar()

        body=Frame(self,bg='Gray80')
        self.initial_focus=self.body(self,body)
        body.pack(padx=5,pady=5)

        self.buttonbox()
        self.grab_set()

        if not self.initial_focus:
            self.initial_focus(self)

        self.protocol("WM_DELETE_WINDOW",self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.wait_window(self)

    def body(self,masterclass,masterwindow):

        self.b1=Checkbutton(masterwindow,text='Use linear bins for 1..10',variable=masterclass.linfirst,state=ACTIVE,bg='Gray80')
        self.b1.grid(row=0,column=0,columnspan=2)
        Label(masterwindow,text='Number of bins:',bg='Gray80').grid(row=1,column=0)
        self.c1=Entry(masterwindow,textvariable=masterclass.numbins,bg='Gray95')
        masterclass.numbins.set('30')
        self.c1.grid(row=1,column=1)
        return self.c1

    def applyme(self):
        self.result=[self.linfirst.get(),float(self.numbins.get())]

class LoadMatrixDialog(MySimpleDialog):
    """Asks for the number of bins for log binning
       and allows linear bins for 1...10"""


    def __init__(self,parent,title='Please provide information:'):

        Toplevel.__init__(self,parent)
       # self.configure(bg='Gray80')
       # self.transient(parent)

        
        self.title(title)

        self.parent=parent
        self.result=None
        self.clones=StringVar()
        self.measuretype=StringVar()

        body=Frame(self)
        self.initial_focus=self.body(self,body)
        body.pack(padx=5,pady=5)

        self.buttonbox()
        self.grab_set()

        if not self.initial_focus:
            self.initial_focus(self)

        self.protocol("WM_DELETE_WINDOW",self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.wait_window(self)

    def body(self,masterclass,masterwindow):

        self.c1=Label(masterwindow,text='What distance measure has been used?',bg='DarkOliveGreen2',anchor=W)
        self.c1.grid(row=0,column=0)

        r1=Radiobutton(masterwindow,text='Non-shared alleles',value='nsa',variable=masterclass.measuretype)
        r2=Radiobutton(masterwindow,text='Linear Manhattan',value='lm',variable=masterclass.measuretype)
        r3=Radiobutton(masterwindow,text='Allele parsimony',value='ap',variable=masterclass.measuretype)
        r4=Radiobutton(masterwindow,text='Hybrid',value="hybrid",variable=masterclass.measuretype)
        r5=Radiobutton(masterwindow,text='Other',value="other",variable=masterclass.measuretype)
        r1.grid(row=1,column=0,sticky=W)
        r2.grid(row=2,column=0,sticky=W)
        r3.grid(row=3,column=0,sticky=W)
        r4.grid(row=4,column=0,sticky=W)
        r5.grid(row=5,column=0,sticky=W)

        self.c2=Label(masterwindow,text='How have clones been handled?',bg='DarkOliveGreen2',anchor=W)
        self.c2.grid(row=6,column=0)
        r6=Radiobutton(masterwindow,text='Removed',value='collapsed',variable=masterclass.clones)
        r7=Radiobutton(masterwindow,text='Kept',value='included',variable=masterclass.clones)
        r8=Radiobutton(masterwindow,text='Unknown',value='unknown',variable=masterclass.clones)
        r6.grid(row=7,column=0,sticky=W)
        r7.grid(row=8,column=0,sticky=W)
        r8.grid(row=9,column=0,sticky=W)
        
        masterclass.measuretype.set('other')
        masterclass.clones.set('unknown')

        return self.c1

    def applyme(self):
        self.result=[self.measuretype.get(),self.clones.get()]

class MetaHelpWindow(MySimpleDialog):


    def __init__(self,parent,title=None,datatype='msat'):

        Toplevel.__init__(self,parent)
        self.configure(bg='Gray80')
        self.transient(parent)
        self.datatype=datatype

        if title:
            self.title=title

        self.parent=parent
        self.result=None
        self.linfirst=IntVar()
        self.numbins=StringVar()

        body=Frame(self,bg='Gray80')
        self.initial_focus=self.body(self,body)
        body.pack(padx=5,pady=5)

        self.buttonbox()
        self.grab_set()

        if not self.initial_focus:
            self.initial_focus(self)

        self.protocol("WM_DELETE_WINDOW",self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.wait_window(self)

    def body(self,masterclass,masterwindow):

        self.text=Text(self,bg='Gray90')
        self.text.pack(expand=YES,fill=BOTH)

        str1="Auxiliary data files are used for reading node properties, such as labels, classes, and sampling sites. "
        str2="File format is ASCII, such that each row lists properties for a node. \n"
        str3="The first row must contain HEADERS, i.e. labels for the properties. \n\n"
        str4="Example for the first row: \n node_label node_site node_latitude node_longitude node_geoclass \n"

        if self.datatype=='net':

            str4=str4+"\nWhen the input data is a network (.edg,.gml), THE FIRST HEADER COLUMN MUST BE node_label, "
            str4=str4+"and there must be a row for each node in the original network, using original node labels."
            str4=str4+"If you have saved the node properties in EDEN Analyzer, this has been taken care of already."

        if self.datatype=='msat':
            str4=str4+"\nThere must be one row for each row in the original microsatellite data file. "

        if self.datatype=="dmat":
            str4=str4+"\nThere must be one row for each row in the original distance matrix file. "

        self.text.insert(INSERT,str1+str2+str3+str4)
        

class AskNumberOfBins(MySimpleDialog):
    """Asks for number of bins for binning"""

    def __init__(self,parent,title=None):

        Toplevel.__init__(self,parent)
        #self.configure(bg='Gray80')
      #  self.transient(parent)

        if title:
            self.title=title

        self.parent=parent
        self.result=None
       # self.linfirst=IntVar()
        self.numbins=StringVar()

        body=Frame(self)
        self.initial_focus=self.body(self,body)
        body.pack(padx=5,pady=5)

        self.buttonbox()
        self.grab_set()

        if not self.initial_focus:
            self.initial_focus(self)

        self.protocol("WM_DELETE_WINDOW",self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.wait_window(self)

    def body(self,masterclass,masterwindow):

       # self.b1=Checkbutton(masterwindow,text='Use linear bins for 1..10',variable=masterclass.linfirst,state=ACTIVE,bg='Gray80')
       # self.b1.grid(row=0,column=0,columnspan=2)
        Label(masterwindow,text='Number of bins:').grid(row=1,column=0)
        self.c1=Entry(masterwindow,textvariable=masterclass.numbins,bg='Gray95')
        masterclass.numbins.set('30')
        self.c1.grid(row=1,column=1)
        return self.c1

    def applyme(self):
        self.result=float(self.numbins.get())

    def validate(self):
        userstr=self.numbins.get()
        try:
            nbins=int(userstr)
        except Exception:
            tkMessageBox.showerror(
                    "Error:",
                    "Number of bins must be an integer.")
            return 0
        if nbins<2:
            tkMessageBox.showerror(
                    "Error:",
                    "Number of bins must be larger than one.")
            return 0
        return 1    
        

class ProjectLaunchDialog(MySimpleDialog):
    """First window shown when launching a new analysis wizard.
        Inquires if the user wants to load microsatellite data,
        a distance matrix, or a network file."""
    

    def __init__(self,parent,title=None):

        Toplevel.__init__(self,parent)
        #self.configure(bg='Gray80')
      #  self.transient(parent)

        self.title("New analysis project")

        self.parent=parent
        self.result=None

        self.datatype=StringVar()
      
       # self.linfirst=IntVar()

        body=Frame(self)
        self.initial_focus=self.body(self,body)
        body.pack(padx=5,pady=5)

        self.buttonbox()
        self.grab_set()

        if not self.initial_focus:
            self.initial_focus(self)

        self.protocol("WM_DELETE_WINDOW",self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.wait_window(self)

    def body(self,masterclass,masterwindow):

       # self.b1=Checkbutton(masterwindow,text='Use linear bins for 1..10',variable=masterclass.linfirst,state=ACTIVE,bg='Gray80')
       # self.b1.grid(row=0,column=0,columnspan=2)

        self.wholeframe=Frame(masterwindow,relief='sunken',borderwidth=2)
        
        self.clabel=Label(self.wholeframe,text="Select input data type:",justify=LEFT,anchor=W,bg='gray90',relief='groove',borderwidth=1)
        self.clabel.pack(side=TOP,expand=YES,fill=X,ipadx=5,ipady=5)

        self.bottompart=Frame(self.wholeframe)

        r1=Radiobutton(self.bottompart,text='Genotype matrix, haploid, individual centred',value='ms_haploid',variable=masterclass.datatype)
        r125=Radiobutton(self.bottompart,text='Genotype matrix, diploid, individual centred',value='ms_diploid',variable=masterclass.datatype)
        r15=Radiobutton(self.bottompart,text='Genotype matrix, haploid, sampling site based',value='mpop_haploid',variable=masterclass.datatype)
        r175=Radiobutton(self.bottompart,text='Genotype matrix, diploid, sampling site based',value='mpop_diploid',variable=masterclass.datatype)
        r1875=Radiobutton(self.bottompart,text='Presence/absence matrix',value='presabs',variable=masterclass.datatype)
        r19=Radiobutton(self.bottompart,text='Presence/abundancy matrix',value='presabu',variable=masterclass.datatype)
        r2=Radiobutton(self.bottompart,text='Distance matrix',value='dmat',variable=masterclass.datatype)
        r3=Radiobutton(self.bottompart,text='Network data',value='net',variable=masterclass.datatype)
        r1.grid(row=1,column=0,sticky=W)
        r125.grid(row=3,column=0,sticky=W)
        r15.grid(row=2,column=0,sticky=W)
        r175.grid(row=4,column=0,sticky=W)
        r1875.grid(row=5,column=0,sticky=W)
        r19.grid(row=6,column=0,sticky=W)
        r2.grid(row=7,column=0,sticky=W)
        r3.grid(row=8,column=0,sticky=W)

        self.bottompart.pack(side=TOP,expand=YES,fill=BOTH,ipadx=7,ipady=7)

        masterclass.datatype.set('ms_haploid')

        self.wholeframe.pack(side=TOP,expand=YES,fill=BOTH)
       
        return self.wholeframe

    def applyme(self):
        self.result=(self.datatype.get())

class ChooseMatrixNodeNames(MySimpleDialog):
    """First window shown when launching a new analysis wizard.
        Inquires if the user wants to load microsatellite data,
        a distance matrix, or a network file."""
    

    def __init__(self,parent,title=None,titlemsg="How to set node labels?"):

        Toplevel.__init__(self,parent)
        #self.configure(bg='Gray80')
      #  self.transient(parent)

        self.title(title)

        self.titlemsg=titlemsg
        self.parent=parent
        self.result=None

        self.measuretype=StringVar()
      
       # self.linfirst=IntVar()

        body=Frame(self)
        self.initial_focus=self.body(self,body)
        body.pack(padx=5,pady=5)

        self.buttonbox()
        self.grab_set()

        if not self.initial_focus:
            self.initial_focus(self)

        self.protocol("WM_DELETE_WINDOW",self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.wait_window(self)

    def body(self,masterclass,masterwindow,titlemsg="How to set node labels?"):


        self.wholeframe=Frame(masterwindow,relief='sunken',borderwidth=2)
        
        self.clabel=Label(self.wholeframe,text=self.titlemsg,justify=LEFT,anchor=W,bg='gray90',relief='groove',borderwidth=1)
        self.clabel.pack(side=TOP,expand=YES,fill=X,ipadx=5,ipady=5)

        self.bottompart=Frame(self.wholeframe)


        r1=Radiobutton(self.bottompart,text='From file',value='file',variable=masterclass.measuretype)
        r2=Radiobutton(self.bottompart,text='1..N',value='numbers',variable=masterclass.measuretype)
        r1.grid(row=1,column=0,sticky=W)
        r2.grid(row=2,column=0,sticky=W)

        self.bottompart.pack(side=TOP,expand=YES,fill=BOTH,ipadx=7,ipady=7)
        
        masterclass.measuretype.set('nsa')

        self.wholeframe.pack(side=TOP,expand=YES,fill=BOTH)
       
        return self.wholeframe

    def applyme(self):
        self.result=(self.measuretype.get())

class ChooseDistanceMeasure(MySimpleDialog):
    """First window shown when launching a new analysis wizard.
        Inquires if the user wants to load microsatellite data,
        a distance matrix, or a network file."""
    

    def __init__(self,parent,title=None,titlemsg="Choose genetic distance measure"):

        Toplevel.__init__(self,parent)
        #self.configure(bg='Gray80')
      #  self.transient(parent)

        self.title(title)

        self.titlemsg=titlemsg
        self.parent=parent
        self.result=None

        self.measuretype=StringVar()
      
       # self.linfirst=IntVar()

        body=Frame(self)
        self.initial_focus=self.body(self,body)
        body.pack(padx=5,pady=5)

        self.buttonbox()
        self.grab_set()

        if not self.initial_focus:
            self.initial_focus(self)

        self.protocol("WM_DELETE_WINDOW",self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.wait_window(self)

    def body(self,masterclass,masterwindow,titlemsg="Choose genetic distance measure"):

       # self.b1=Checkbutton(masterwindow,text='Use linear bins for 1..10',variable=masterclass.linfirst,state=ACTIVE,bg='Gray80')
       # self.b1.grid(row=0,column=0,columnspan=2)

        self.wholeframe=Frame(masterwindow,relief='sunken',borderwidth=2)
        
        self.clabel=Label(self.wholeframe,text=self.titlemsg,justify=LEFT,anchor=W,bg='gray90',relief='groove',borderwidth=1)
        self.clabel.pack(side=TOP,expand=YES,fill=X,ipadx=5,ipady=5)

        self.bottompart=Frame(self.wholeframe)


        #r1=Radiobutton(self.bottompart,text='Non-shared alleles',value='nsa',variable=masterclass.measuretype)
        r1=Radiobutton(self.bottompart,text='Allele Sharing',value='ap',variable=masterclass.measuretype)
        r2=Radiobutton(self.bottompart,text='Linear Manhattan',value='lm',variable=masterclass.measuretype)
        #r3=Radiobutton(self.bottompart,text='Allele parsimony',value='ap',variable=masterclass.measuretype)
        #r4=Radiobutton(self.bottompart,text='Hybrid',value="hybrid",variable=masterclass.measuretype)
        r1.grid(row=1,column=0,sticky=W)
        r2.grid(row=2,column=0,sticky=W)
        #r3.grid(row=3,column=0,sticky=W)
        #r4.grid(row=4,column=0,sticky=W)

        self.bottompart.pack(side=TOP,expand=YES,fill=BOTH,ipadx=7,ipady=7)
        
        masterclass.measuretype.set('nsa')

        self.wholeframe.pack(side=TOP,expand=YES,fill=BOTH)
       
        return self.wholeframe

    def applyme(self):
        self.result=(self.measuretype.get())

class ImportMetadataYesNo(MySimpleDialog):
    """First window shown when launching a new analysis wizard.
        Inquires if the user wants to load microsatellite data,
        a distance matrix, or a network file."""
    

    def __init__(self,parent,title=None,titlemsg="Do you want to import auxiliary node data?",datatype='msat'):

        Toplevel.__init__(self,parent)
        #self.configure(bg='Gray80')
      #  self.transient(parent)

        self.title(title)

        self.datatype=datatype

        self.titlemsg=titlemsg
        self.parent=parent
        self.result=None

        self.metatype=IntVar()
      
       # self.linfirst=IntVar()

        body=Frame(self)
        self.initial_focus=self.body(self,body)
        body.pack(padx=5,pady=5)

        self.buttonbox(self.datatype)
        self.grab_set()

        if not self.initial_focus:
            self.initial_focus(self)

        self.protocol("WM_DELETE_WINDOW",self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.wait_window(self)

    def buttonbox(self,datatype='msat'):
        """OK, Cancel and Help buttons"""

        box=Frame(self)
        w=Button(box,text="OK",width=10,command=self.ok,default=ACTIVE)
        w.pack(side=LEFT,padx=5,pady=5)
        w=Button(box,text="Cancel",width=10,command=self.cancel)
        w.pack(side=LEFT,padx=5,pady=5)
        w=Button(box,text="Help",width=10,command=lambda s=self,t=datatype: s.displayhelp(t))
        w.pack(side=LEFT,padx=5,pady=5)

        self.bind("&lt;Return>",self.ok)
        self.bind("&lt;Escape",self.cancel)

        box.pack()

    def body(self,masterclass,masterwindow,titlemsg="Do you want to import auxiliary node data?"):

       # self.b1=Checkbutton(masterwindow,text='Use linear bins for 1..10',variable=masterclass.linfirst,state=ACTIVE,bg='Gray80')
       # self.b1.grid(row=0,column=0,columnspan=2)

        self.wholeframe=Frame(masterwindow,relief='sunken',borderwidth=2)
        
        self.clabel=Label(self.wholeframe,text=self.titlemsg,justify=LEFT,anchor=W,bg='gray90',relief='groove',borderwidth=1)
        self.clabel.pack(side=TOP,expand=YES,fill=X,ipadx=5,ipady=5)

        self.bottompart=Frame(self.wholeframe)

        r1=Radiobutton(self.bottompart,text='Yes',value=1,variable=masterclass.metatype)
        r2=Radiobutton(self.bottompart,text='No',value=0,variable=masterclass.metatype)
  
        r1.grid(row=1,column=0,sticky=W)
        r2.grid(row=2,column=0,sticky=W)
     
        masterclass.metatype.set(1)

        self.bottompart.pack(side=TOP,expand=YES,fill=BOTH,ipadx=7,ipady=7)

        self.wholeframe.pack(side=TOP,expand=YES,fill=BOTH)
       
        return self.wholeframe

    def applyme(self):
        self.result=(self.metatype.get())

    def displayhelp(self,datatype):
        MetaHelpWindow(self,datatype)




class MatrixDialog(MySimpleDialog):
    """Used when loading a matrix. Asks if the matrix contains weights or distances"""
    

    def __init__(self,parent,title=None):

        Toplevel.__init__(self,parent)
        #self.configure(bg='Gray80')
      #  self.transient(parent)

        if title:
            self.title=title

        self.parent=parent
        self.result=None

        self.mattype=IntVar()
      
       # self.linfirst=IntVar()

        body=Frame(self)
        self.initial_focus=self.body(self,body)
        body.pack(padx=5,pady=5)

        self.buttonbox()
        self.grab_set()

        if not self.initial_focus:
            self.initial_focus(self)

        self.protocol("WM_DELETE_WINDOW",self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.wait_window(self)

    def body(self,masterclass,masterwindow):

       # self.b1=Checkbutton(masterwindow,text='Use linear bins for 1..10',variable=masterclass.linfirst,state=ACTIVE,bg='Gray80')
       # self.b1.grid(row=0,column=0,columnspan=2)
        
        self.c1=Label(masterwindow,text='Matrix type:')
        self.c1.grid(row=0,column=0)

        r1=Radiobutton(masterwindow,text='Weight matrix',value=1,variable=masterclass.mattype)
        r2=Radiobutton(masterwindow,text='Distance matrix',value=0,variable=masterclass.mattype)
        r1.grid(row=0,column=1,sticky=W)
        r2.grid(row=1,column=1,sticky=W)
        
        masterclass.mattype.set(0)
       
        return self.c1

    def applyme(self):
        self.result=(self.mattype.get())

class MsatDialog(MySimpleDialog):
    """Used when loading a matrix. Asks if the matrix contains weights or distances"""
    

    def __init__(self,parent,title=None,titlemsg="Handling clones"):

        Toplevel.__init__(self,parent)
        #self.configure(bg='Gray80')
      #  self.transient(parent)

        self.title(title)

        self.titlemsg=titlemsg

        self.parent=parent
        self.result=None

        self.mattype=IntVar()
      
       # self.linfirst=IntVar()

        body=Frame(self)
        self.initial_focus=self.body(self,body)
        body.pack(padx=5,pady=5)

        self.buttonbox()
        self.grab_set()

        if not self.initial_focus:
            self.initial_focus(self)

        self.protocol("WM_DELETE_WINDOW",self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.wait_window(self)

    def body(self,masterclass,masterwindow,titlemsg="Handling clones"):

       # self.b1=Checkbutton(masterwindow,text='Use linear bins for 1..10',variable=masterclass.linfirst,state=ACTIVE,bg='Gray80')
       # self.b1.grid(row=0,column=0,columnspan=2)

        self.wholeframe=Frame(masterwindow,relief='sunken',borderwidth=2)
        
        self.clabel=Label(self.wholeframe,text=self.titlemsg,justify=LEFT,anchor=W,bg='gray90',relief='groove',borderwidth=1)
        self.clabel.pack(side=TOP,expand=YES,fill=X,ipadx=5,ipady=5)

        self.bottompart=Frame(self.wholeframe)

      
        r1=Radiobutton(self.bottompart,text='Collapse clones',value=1,variable=masterclass.mattype)
        r2=Radiobutton(self.bottompart,text='Leave clones',value=0,variable=masterclass.mattype)
        r1.grid(row=1,column=0,sticky=W)
        r2.grid(row=2,column=0,sticky=W)

        self.bottompart.pack(side=TOP,expand=YES,fill=BOTH,ipadx=7,ipady=7)

        self.wholeframe.pack(side=TOP,expand=YES,fill=BOTH)
        
        masterclass.mattype.set(1)
       
        return self.wholeframe

    def applyme(self):
        self.result=(self.mattype.get())

        
class VisualizationDialog(MySimpleDialog):
    """Asks options for network visualization"""

    def __init__(self,parent,title=None):

        Toplevel.__init__(self,parent)
        #self.configure(bg='Gray80')
      #  self.transient(parent)

        if title:
            self.title=title

        self.parent=parent
        self.result=None

        self.winsize=IntVar()
        self.vtxsize=StringVar()
        self.vtxcolor=StringVar()
        self.bgcolor=StringVar()
        self.showlabels=StringVar()
       # self.linfirst=IntVar()

        body=Frame(self)
        self.initial_focus=self.body(self,body)
        body.pack(padx=5,pady=5)

        self.buttonbox()
        self.grab_set()

        if not self.initial_focus:
            self.initial_focus(self)

        self.protocol("WM_DELETE_WINDOW",self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.wait_window(self)

    def body(self,masterclass,masterwindow):

       # self.b1=Checkbutton(masterwindow,text='Use linear bins for 1..10',variable=masterclass.linfirst,state=ACTIVE,bg='Gray80')
       # self.b1.grid(row=0,column=0,columnspan=2)
       
        self.c1=Label(masterwindow,text='Vertex color:')
        self.c1.grid(row=0,column=0)
        rowcount=-1
        for text, value in [('Black','000000'),('White','999999'),('Red','990000'),('Green','009900'),('Blue','000099'),('By strength','-1')]:
            rowcount+=1
            Radiobutton(masterwindow,text=text,value=value,variable=masterclass.vtxcolor).grid(row=rowcount,column=1,sticky=W)
        masterclass.vtxcolor.set('-1')
        Label(masterwindow,text='Vertex size:').grid(row=rowcount+1,column=0)
        for text, value in [('Small','0.4'),('Medium','0.7'),('Large','0.99'),('By strength','-1.0')]:
            rowcount=rowcount+1
            Radiobutton(masterwindow,text=text,value=value,variable=masterclass.vtxsize).grid(row=rowcount,column=1,sticky=W)
        masterclass.vtxsize.set('-1.0')
        Label(masterwindow,text='Show with:').grid(row=rowcount+1,column=0)
        for text, value in [('White background','white'),('Black background','black')]:
            rowcount=rowcount+1
            Radiobutton(masterwindow,text=text,value=value,variable=masterclass.bgcolor).grid(row=rowcount,column=1,sticky=W)
        masterclass.bgcolor.set('black')
        Label(masterwindow,text="Vertex labels:").grid(row=rowcount+1,column=0)
        for text, value in [('None','none'),('All','all'),('Top 10','top10')]:
            rowcount=rowcount+1
            Radiobutton(masterwindow,text=text,value=value,variable=masterclass.showlabels).grid(row=rowcount,column=1,sticky=W)
            masterclass.showlabels.set('all')
        return self.c1

    def applyme(self):
        self.result=(self.vtxcolor.get(),self.vtxsize.get(),self.bgcolor.get(),self.showlabels.get())

class AskThreshold(MySimpleDialog):
    """Asks threshold for thresholding"""

    def __init__(self,parent,title=None):

        Toplevel.__init__(self,parent)
        #self.configure(bg='Gray80')
      #  self.transient(parent)

        if title:
            self.title=title

        self.parent=parent
        self.result=None
       # self.linfirst=IntVar()
        self.threshold=StringVar()

        body=Frame(self)
        self.initial_focus=self.body(self,body)
        body.pack(padx=5,pady=5)

        self.buttonbox()
        self.grab_set()

        if not self.initial_focus:
            self.initial_focus(self)

        self.protocol("WM_DELETE_WINDOW",self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.wait_window(self)

    def body(self,masterclass,masterwindow):

       # self.b1=Checkbutton(masterwindow,text='Use linear bins for 1..10',variable=masterclass.linfirst,state=ACTIVE,bg='Gray80')
       # self.b1.grid(row=0,column=0,columnspan=2)
        Label(masterwindow,text='Threshold:').grid(row=1,column=0)
        self.c1=Entry(masterwindow,textvariable=masterclass.threshold,bg='Gray95')
        masterclass.threshold.set('0')
        self.c1.grid(row=1,column=1)
        return self.c1

    def applyme(self):
        self.result=float(self.threshold.get())

class PercolationDialog(MySimpleDialog):
    """First window shown when launching a new analysis wizard.
        Inquires if the user wants to load microsatellite data,
        a distance matrix, or a network file."""
    

    def __init__(self,parent,title=None,titlemsg="Set threshold distance",pdata=[],suscmax_thresh=0.0):

        Toplevel.__init__(self,parent)
        #self.configure(bg='Gray80')
      #  self.transient(parent)

        self.title(title)

        self.titlemsg=titlemsg
        self.parent=parent
        self.result=None
        self.data=pdata
        self.default_thresh=suscmax_thresh

        self.threshold=StringVar()
      
       # self.linfirst=IntVar()

        body=Frame(self)
        self.initial_focus=self.body(self,body)
        body.pack(padx=5,pady=5)

        self.buttonbox()
        self.grab_set()

        if not self.initial_focus:
            self.initial_focus(self)

        self.protocol("WM_DELETE_WINDOW",self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.wait_window(self)

    def body(self,masterclass,masterwindow,titlemsg="Choose genetic distance measure"):

       # self.b1=Checkbutton(masterwindow,text='Use linear bins for 1..10',variable=masterclass.linfirst,state=ACTIVE,bg='Gray80')
       # self.b1.grid(row=0,column=0,columnspan=2)

        self.wholeframe=Frame(masterwindow,relief='groove',borderwidth=2,bg="Gray95")
        
       # self.clabel=Label(self.wholeframe,text=self.titlemsg,bg='DarkOliveGreen2',relief='groove',borderwidth=1)
       # self.clabel.pack(side=TOP,expand=YES,fill=X,ipadx=5,ipady=5)

        self.midpart=Frame(self.wholeframe)

        myplot=visuals.ReturnPlotObject(self.data,plotcommand="plot",addstr=",color=\"#9e0b0f\"",titlestring="Largest component size:Susceptibility",xstring="Threshold distance",ystring="GCC size:Susceptibility",fontsize=9)
        myplot.canvas=FigureCanvasTkAgg(myplot.thisFigure,master=self.midpart)
        myplot.canvas.show()
        myplot.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=YES, padx=10,pady=10) 

        self.midpart.pack(side=TOP,expand=YES,fill=BOTH)

        self.bottompart=Frame(self.wholeframe)
        
        Label(self.bottompart,text='Estimated percolation threshold = %2.2f' % self.default_thresh).grid(row=1,column=0,columnspan=2)
        Label(self.bottompart,text='Threshold: ').grid(row=2,column=0)
        self.c1=Entry(self.bottompart,textvariable=masterclass.threshold,bg='Gray95')
        masterclass.threshold.set(str(self.default_thresh))
        self.c1.grid(row=2,column=1)
        self.bottompart.pack(side=TOP,expand=YES,fill=BOTH,ipadx=7,ipady=7)

        self.wholeframe.pack(side=TOP,expand=YES,fill=BOTH)
       
        return self.wholeframe

    def applyme(self):
        self.result=(self.threshold.get())
        

class VisualizationOptions(MySimpleDialog):
    """First window shown when launching a new analysis wizard.
        Inquires if the user wants to load microsatellite data,
        a distance matrix, or a network file."""
    

    def __init__(self,parent,network,title=None,titlemsg="Choose visualization options"):

        Toplevel.__init__(self,parent)
        #self.configure(bg='Gray80')
      #  self.transient(parent)

        self.title(title)

        self.titlemsg=titlemsg
        self.parent=parent
        self.result=None

        self.vcolor=StringVar()
        self.vsize=StringVar()
        self.bgcolor=IntVar()
      
       # self.linfirst=IntVar()

        body=Frame(self)
        self.initial_focus=self.body(self,body)
        body.pack(padx=5,pady=5)

        self.buttonbox(self.datatype)
        self.grab_set()

        if not self.initial_focus:
            self.initial_focus(self)

        self.protocol("WM_DELETE_WINDOW",self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.wait_window(self)

    def buttonbox(self,datatype='msat'):
        """OK, Cancel and Help buttons"""

        box=Frame(self)
        w=Button(box,text="OK",width=10,command=self.ok,default=ACTIVE)
        w.pack(side=LEFT,padx=5,pady=5)
        w=Button(box,text="Cancel",width=10,command=self.cancel)
        w.pack(side=LEFT,padx=5,pady=5)

        self.bind("&lt;Return>",self.ok)
        self.bind("&lt;Escape",self.cancel)

        box.pack()

    def body(self,masterclass,masterwindow,titlemsg="Select visualization options"):

       # self.b1=Checkbutton(masterwindow,text='Use linear bins for 1..10',variable=masterclass.linfirst,state=ACTIVE,bg='Gray80')
       # self.b1.grid(row=0,column=0,columnspan=2)
       

        self.wholeframe=Frame(masterwindow,relief='sunken',borderwidth=2)
        
        self.clabel=Label(self.wholeframe,text=self.titlemsg,justify=LEFT,anchor=W,bg='gray90',relief='groove',borderwidth=1)
        self.clabel.pack(side=TOP,expand=YES,fill=X,ipadx=5,ipady=5)

        self.bottompart=Frame(self.wholeframe)

        scrollbar = Scrollbar(self.bottompart, orient=VERTICAL)
        listbox = Listbox(self.bottompart, yscrollcommand=scrollbar.set)
        scrollbar.config(command=listbox.yview)
        scrollbar.pack(side=RIGHT, fill=Y)
        listbox.pack(side=LEFT, fill=BOTH, expand=1)
        
        listbox.insert(END,'none')
        listbox.insert(END,'by degree')

        plist=netext.getNumericProperties(network)

        for prop in plist:
            listbox.insert(END,prop)
       
     
        self.bottompart.pack(side=TOP,expand=YES,fill=BOTH,ipadx=7,ipady=7)

        self.wholeframe.pack(side=TOP,expand=YES,fill=BOTH)
       
        return self.wholeframe

    def applyme(self):
        #self.result=(self.metatype.get())
        pass

    def displayhelp(self,datatype):
        MetaHelpWindow(self,datatype)


class SliderDialog(MySimpleDialog):
    """Dialog for inputting a value using a slider (e.g. for font sizes etc)"""
    

    def __init__(self,parent,title=None,titlemsg="Select label font size",minval=0,maxval=100,currval=7):

        Toplevel.__init__(self,parent)
        #self.configure(bg='Gray80')
      #  self.transient(parent)

        self.title(title)

        self.titlemsg=titlemsg

        self.parent=parent
        self.result=None

        self.slidertype=IntVar()
        self.minval=minval
        self.maxval=maxval
        self.currval=currval
      
       # self.linfirst=IntVar()

        body=Frame(self)
        self.initial_focus=self.body(self,body)
        body.pack(padx=5,pady=5)

        self.buttonbox()
        self.grab_set()

        if not self.initial_focus:
            self.initial_focus(self)

        self.protocol("WM_DELETE_WINDOW",self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.wait_window(self)

    def body(self,masterclass,masterwindow,titlemsg="Select label font size"):

       # self.b1=Checkbutton(masterwindow,text='Use linear bins for 1..10',variable=masterclass.linfirst,state=ACTIVE,bg='Gray80')
       # self.b1.grid(row=0,column=0,columnspan=2)

        self.wholeframe=Frame(masterwindow,relief='sunken',borderwidth=2)
        
        self.clabel=Label(self.wholeframe,text=self.titlemsg,justify=LEFT,anchor=W,bg='gray90',relief='groove',borderwidth=1)
        self.clabel.pack(side=TOP,expand=YES,fill=X,ipadx=5,ipady=5)

        self.bottompart=Frame(self.wholeframe)

      
        self.scale=Scale(self.bottompart, from_=self.minval, to=self.maxval, orient=HORIZONTAL)
        self.scale.set(self.currval)


        self.scale.pack()

        self.bottompart.pack(side=TOP,expand=YES,fill=BOTH,ipadx=7,ipady=7)

        self.wholeframe.pack(side=TOP,expand=YES,fill=BOTH)
       
        return self.wholeframe

    def applyme(self):
        self.result=(self.scale.get())

class DoubleSliderDialog(MySimpleDialog):
    """Dialog for inputting a value using a slider (e.g. for font sizes etc)"""
    

    def __init__(self,parent,title=None,titlemsg="Select label font size",resolution=1,minval_1=0,maxval_1=100,currval_1=7,minval_2=0,maxval_2=0,currval_2=7,slidertext_1='min',slidertext_2='max'):

        Toplevel.__init__(self,parent)
        #self.configure(bg='Gray80')
      #  self.transient(parent)

        self.title(title)

        self.titlemsg=titlemsg

        self.parent=parent
        self.result=None

        self.slidertype_1=IntVar()
        self.minval_1=minval_1
        self.maxval_1=maxval_1
        self.currval_1=currval_1
        self.slidertext_1=slidertext_1
        self.resolution=resolution

        self.slidertype_2=IntVar()
        self.minval_2=minval_2
        self.maxval_2=maxval_2
        self.currval_2=currval_2
        self.slidertext_2=slidertext_2
      
       # self.linfirst=IntVar()

        body=Frame(self)
        self.initial_focus=self.body(self,body)
        body.pack(padx=5,pady=5)

        self.buttonbox()
        self.grab_set()

        if not self.initial_focus:
            self.initial_focus(self)

        self.protocol("WM_DELETE_WINDOW",self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.wait_window(self)

    def body(self,masterclass,masterwindow,titlemsg="Select label font size"):

       # self.b1=Checkbutton(masterwindow,text='Use linear bins for 1..10',variable=masterclass.linfirst,state=ACTIVE,bg='Gray80')
       # self.b1.grid(row=0,column=0,columnspan=2)

        self.wholeframe=Frame(masterwindow,relief='sunken',borderwidth=2)
        
        self.clabel=Label(self.wholeframe,text=self.titlemsg,justify=LEFT,anchor=W,bg='gray90',relief='groove',borderwidth=1)
        self.clabel.pack(side=TOP,expand=YES,fill=X,ipadx=5,ipady=5)

        self.bottompart=Frame(self.wholeframe)

        Label(self.bottompart,text=self.slidertext_1,justify=LEFT,anchor=W,bg='gray90',relief='groove',borderwidth=1).grid(row=0,column=0)
          
        self.scale_1=Scale(self.bottompart, from_=self.minval_1, to=self.maxval_1, orient=HORIZONTAL,resolution=self.resolution)
        self.scale_1.set(self.currval_1)


        self.scale_1.grid(row=0,column=1)

        Label(self.bottompart,text=self.slidertext_2,justify=LEFT,anchor=W,bg='gray90',relief='groove',borderwidth=1).grid(row=1,column=0)
          
        self.scale_2=Scale(self.bottompart, from_=self.minval_2, to=self.maxval_2, orient=HORIZONTAL,resolution=self.resolution)
        self.scale_2.set(self.currval_2)


        self.scale_2.grid(row=1,column=1)

        

        self.bottompart.pack(side=TOP,expand=YES,fill=BOTH,ipadx=7,ipady=7)

        self.wholeframe.pack(side=TOP,expand=YES,fill=BOTH)
       
        return self.wholeframe

    def applyme(self):
        self.result=[self.scale_1.get(),self.scale_2.get()]

class ColorMapDialog(MySimpleDialog):
    """Lists color maps and asks to choose one"""


    def __init__(self,parent,title='Select color map:',current_map='orange'):

        Toplevel.__init__(self,parent)
       # self.configure(bg='Gray80')
       # self.transient(parent)

        
        self.title(title)

        self.parent=parent
        self.result=None
        self.current_map=current_map
        self.colormap=StringVar()
   

        body=Frame(self)
        self.initial_focus=self.body(self,body)
        body.pack(padx=5,pady=5)

        self.buttonbox()
        self.grab_set()

        if not self.initial_focus:
            self.initial_focus(self)

        self.protocol("WM_DELETE_WINDOW",self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.wait_window(self)

    def body(self,masterclass,masterwindow):

        self.c1=Label(masterwindow,text='Please choose color map',anchor=W)
        self.c1.grid(row=0,column=0)

        colormaps=["orange","hsv","jet","spectral","winter","Accent","Paired","binary","gray"]

        rowcounter=1
        curr_frame=[]
        curr_plot=[]

        for cmap in colormaps:

            Radiobutton(masterwindow,text=cmap,value=cmap,variable=masterclass.colormap).grid(row=rowcounter,column=0,sticky=W)
            

            curr_frame.append(Frame(masterwindow))

            curr_plot.append(visuals.ReturnColorMapPlot(cmap))
            curr_plot[-1].canvas=FigureCanvasTkAgg(curr_plot[-1],curr_frame[-1])
            curr_plot[-1].canvas.show()
            curr_plot[-1].canvas.get_tk_widget().grid(row=rowcounter,column=1)

            curr_frame[-1].grid(row=rowcounter,column=1,sticky=W)

            rowcounter=rowcounter+1
            
        if self.current_map in colormaps:
            masterclass.colormap.set(self.current_map)
        else:
            masterclass.colormap.set('orange')
       

        return self.c1

    def applyme(self):
        self.result=self.colormap.get()
        
class WaitWindow(MySimpleDialog):
    """Used when loading a matrix. Asks if the matrix contains weights or distances"""
    

    def __init__(self,parent,title=None,titlemsg="Please wait...",hasProgressbar=False):

        Toplevel.__init__(self,parent)
        #self.configure(bg='Gray80')
        #self.transient(parent)

        self.title(title)
        self.titlemsg=titlemsg
        self.parent=parent
        self.hasProgressbar=hasProgressbar

        self.result=None

        body=Frame(self)
        self.initial_focus=self.body(self,body)
        body.pack(padx=5,pady=5)

        #self.buttonbox()
        #self.grab_set()

        if not self.initial_focus:
            self.initial_focus(self)

        self.protocol("WM_DELETE_WINDOW",self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,parent.winfo_rooty()+50))

        self.initial_focus.focus_set()



    def go(self,lambda_operation):
        self.result=lambda_operation()
        self.ok()
        #self.wait_window(self)

    def body(self,masterclass,masterwindow,titlemsg="Please wait..."):
        self.wholeframe=Frame(masterwindow,relief='sunken',borderwidth=2)
        self.clabel=Label(self.wholeframe,text=self.titlemsg,
                          justify=LEFT,
                          anchor=W,
                          bg='darkolivegreen2',
                          relief='groove',borderwidth=1,padx=3,pady=3,width=40)
        self.clabel.pack(side=TOP,expand=YES,fill=X,ipadx=5,ipady=5)
        self.bottompart=Frame(self.wholeframe)
        self.l1=Label(self.bottompart,text=self.titlemsg,justify=LEFT,
                      anchor=W,bg='white',relief='flat',padx=3,pady=3)
        self.l1.grid(row=1,column=0,sticky=W)
        self.l1['fg']='black'

        if self.hasProgressbar:
            self.progressbar=Meter(self.bottompart,value=0.0)
            self.progressbar.grid(row=4,column=0,sticky=W)
            self.updatePointer=self.progressbar.set

        self.bottompart.pack(side=TOP,expand=YES,fill=BOTH,ipadx=7,ipady=7)
        self.wholeframe.pack(side=TOP,expand=YES,fill=BOTH)

        return self.wholeframe



    def applyme(self):

        return self.result

class MsLoadWaiter(MySimpleDialog):
    """Used when loading a matrix. Asks if the matrix contains weights or distances"""
    

    def __init__(self,parent,inputfile,removeclones,measuretype,title="Processing microsatellite data",titlemsg="Please wait",nodeNames=None,datatype=None):

        Toplevel.__init__(self,parent)
        
        self.title(title)
        self.titlemsg=titlemsg
        self.parent=parent

        self.datatype=datatype

        self.result=None

        body=Frame(self)
        self.initial_focus=self.body(self,body)
        body.pack(padx=5,pady=5)

        #self.buttonbox()
        self.grab_set()
        self.points=['.','..','...']
        self.pointcounter=0

        self.displayBusyCursor()
        
        if not self.initial_focus:
            self.initial_focus(self)

        self.protocol("WM_DELETE_WINDOW",self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.l1['fg']='#b0b0b0'
        self.l1.update()

        if removeclones:
            self.l2['text']='Removing clones...'
            self.l2['fg']='black'
            self.l2.update()

        if datatype=="ms_diploid":
            msdata_initial=eden.MicrosatelliteData(inputfile)
        else:
            msdata_initial=eden.MicrosatelliteDataHaploid(inputfile)

        if removeclones:

            [msdata,keeptheserows]=msdata_initial.getUniqueSubset(returnOldIndices=True)
            clones='collapsed'

            if nodeNames!=None:
                newNodeNames=[]
                keeptheserows=sorted(keeptheserows)
                for row in keeptheserows: 
                    newNodeNames.append(nodeNames[row])
                nodeNames=newNodeNames
            else:
                nodeNames=sorted(keeptheserows)
                # trick to accommodate for (possible) node properties - we will use only "keeptheserows" containing indices to non-clonal samples
        else:
            msdata=msdata_initial
            clones='included'
            keeptheserows=None

        self.clones=clones
        self.keeptheserows=keeptheserows


        
        self.l3['text']='Calculating distance matrix...'
        self.l3['fg']='black'
        self.l3.update()

        self.l2['fg']='#b0b0b0'
        self.l2.update()
    
        m=msdata.getDistanceMatrix(measuretype,nodeNames=nodeNames,progressUpdater=self.progressbar.set)        

        if removeclones:
            Nclones=msdata_initial.getNumberOfNodes()-len(keeptheserows)
        else:
            Nclones=None
        self.Nclones=Nclones
        self.m=m
        self.msdata=msdata

        self.ok()


    def animate(self):

        self.pointcounter=self.pointcounter+1
        if self.pointcounter==3:
            self.pointcounter=0

        self.clabel['text']=self.titlemsg+self.points[self.pointcounter]
        self.clabel.update()

        self.after(500,self.animate)

    def body(self,masterclass,masterwindow,titlemsg="Please wait!"):

       # self.b1=Checkbutton(masterwindow,text='Use linear bins for 1..10',variable=masterclass.linfirst,state=ACTIVE,bg='Gray80')
       # self.b1.grid(row=0,column=0,columnspan=2)

        self.wholeframe=Frame(masterwindow,relief='sunken',borderwidth=2)
        self.clabel=Label(self.wholeframe,text=self.titlemsg,
                          justify=LEFT,
                          anchor=W,
                          bg='darkolivegreen2',
                          relief='groove',borderwidth=1,padx=3,pady=3,width=40)
        self.clabel.pack(side=TOP,expand=YES,fill=X,ipadx=5,ipady=5)
        self.bottompart=Frame(self.wholeframe)
        self.l1=Label(self.bottompart,text='Reading file...',justify=LEFT,
                      anchor=W,bg='white',relief='flat',padx=3,pady=3)
        self.l1.grid(row=1,column=0,sticky=W)
        self.l1['fg']='black'
        self.l2=Label(self.bottompart,text='-',justify=LEFT,
                      anchor=W,bg='white',relief='flat',padx=3,pady=3)
        self.l2.grid(row=2,column=0,sticky=W)
        self.l3=Label(self.bottompart,text='-',justify=LEFT,
                      anchor=W,bg='white',relief='flat',padx=3,pady=3)
        self.l3.grid(row=3,column=0,sticky=W)

        self.progressbar=Meter(self.bottompart,value=0.0)
        self.progressbar.grid(row=4,column=0,sticky=W)

        self.bottompart.pack(side=TOP,expand=YES,fill=BOTH,ipadx=7,ipady=7)
        self.wholeframe.pack(side=TOP,expand=YES,fill=BOTH)
               
        return self.wholeframe

    def applyme(self):

        self.result=[self.Nclones,self.clones,self.keeptheserows,self.m,self.msdata]



class LoadWaiter(MySimpleDialog):
    """Display progress bar while doing something."""
    

    def __init__(self,parent,function,args=[],kwargs=None,callback_function_parameter="progressUpdater",bodymsg="Processing...",titlemsg="Please wait"):
        Toplevel.__init__(self,parent)

        title=titlemsg

        if kwargs==None:
            kwargs={}

        self.title(title)
        self.titlemsg=titlemsg
        self.parent=parent
        self.result=None

        body=Frame(self)
        self.initial_focus=self.body(self,body,titlemsg=titlemsg,bodymsg=bodymsg)
        body.pack(padx=5,pady=5)

        #self.buttonbox()
        self.grab_set()
        self.points=['.','..','...']
        self.pointcounter=0

        self.displayBusyCursor()
        
        if not self.initial_focus:
            self.initial_focus(self)

        self.protocol("WM_DELETE_WINDOW",self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.l1['fg']='#b0b0b0'
        self.l1.update()

        kwargs[callback_function_parameter]=self.progressbar.set
        self.output=function(*args,**kwargs)        

        self.ok()


    def animate(self):

        self.pointcounter=self.pointcounter+1
        if self.pointcounter==3:
            self.pointcounter=0

        self.clabel['text']=self.titlemsg+self.points[self.pointcounter]
        self.clabel.update()

        self.after(500,self.animate)

    def body(self,masterclass,masterwindow,titlemsg="",bodymsg=""):
        self.wholeframe=Frame(masterwindow,relief='sunken',borderwidth=2)
        self.clabel=Label(self.wholeframe,text=self.titlemsg,
                          justify=LEFT,
                          anchor=W,
                          bg='darkolivegreen2',
                          relief='groove',borderwidth=1,padx=3,pady=3,width=40)
        self.clabel.pack(side=TOP,expand=YES,fill=X,ipadx=5,ipady=5)
        self.bottompart=Frame(self.wholeframe)
        self.l1=Label(self.bottompart,text=bodymsg,justify=LEFT,
                      anchor=W,bg='white',relief='flat',padx=3,pady=3)
        self.l1.grid(row=1,column=0,sticky=W)
        self.l1['fg']='black'
        """
        self.l2=Label(self.bottompart,text='-',justify=LEFT,
                      anchor=W,bg='white',relief='flat',padx=3,pady=3)
        self.l2.grid(row=2,column=0,sticky=W)
        self.l3=Label(self.bottompart,text='-',justify=LEFT,
                      anchor=W,bg='white',relief='flat',padx=3,pady=3)
        self.l3.grid(row=3,column=0,sticky=W)
        """
        self.progressbar=Meter(self.bottompart,value=0.0)
        self.progressbar.grid(row=4,column=0,sticky=W)

        self.bottompart.pack(side=TOP,expand=YES,fill=BOTH,ipadx=7,ipady=7)
        self.wholeframe.pack(side=TOP,expand=YES,fill=BOTH)
               
        return self.wholeframe

    def applyme(self):
        self.result=self.output



class GenericLoadWaiter(MySimpleDialog):
    """Generic class for wait windows displaying one item"""
    

    def __init__(self,parent,filename,removeclones,measuretype,title=None,titlemsg="Please wait",specmsg="Processing"):

        Toplevel.__init__(self,parent)
        #self.configure(bg='Gray80')
      #  self.transient(parent)

        self.title(title)
            
        self.titlemsg=titlemsg
        self.specmsg=specmsg

        self.parent=parent
        self.result=None

        body=Frame(self)
        self.initial_focus=self.body(self,body)
        body.pack(padx=5,pady=5)

       #self.buttonbox()
        self.grab_set()

        self.points=['.','..','...']
        self.pointcounter=0

        self.displayBusyCursor()
        
        if not self.initial_focus:
            self.initial_focus(self)

        self.protocol("WM_DELETE_WINDOW",self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        # INSERT YOUR FUNCTIONALITY HERE

        self.ok()

    def body(self,masterclass,masterwindow,titlemsg="Please wait!"):

        self.wholeframe=Frame(masterwindow,relief='sunken',borderwidth=2)
        
        self.clabel=Label(self.wholeframe,text=self.titlemsg,justify=LEFT,anchor=W,bg='darkolivegreen2',relief='groove',borderwidth=1,padx=3,pady=3,width=40)
        self.clabel.pack(side=TOP,expand=YES,fill=X,ipadx=5,ipady=5)

        self.bottompart=Frame(self.wholeframe)

        self.l1=Label(self.bottompart,text=self.specmsg,justify=LEFT,anchor=W,bg='white',relief='flat',padx=3,pady=3)
        self.l1.grid(row=1,column=0,sticky=W)
        self.l1['fg']='black'

        self.bottompart.pack(side=TOP,expand=YES,fill=BOTH,ipadx=7,ipady=7)

        self.wholeframe.pack(side=TOP,expand=YES,fill=BOTH)
               
        return self.wholeframe


class MetaLoadWaiter(GenericLoadWaiter):
    """Generic class for wait windows displaying one item"""
    

    def __init__(self,parent,m,keeptheserows,title=None,titlemsg="Please wait",specmsg="Processing"):

        Toplevel.__init__(self,parent)
        #self.configure(bg='Gray80')
      #  self.transient(parent)

        self.title(title)

        self.titlemsg=titlemsg
        self.specmsg=specmsg

        self.parent=parent
        self.result=None

        body=Frame(self)
        self.initial_focus=self.body(self,body)
        body.pack(padx=5,pady=5)

       # self.buttonbox()
        self.grab_set()

        self.points=['.','..','...']
        self.pointcounter=0

        self.displayBusyCursor()
        
        if not self.initial_focus:
            self.initial_focus(self)

        self.protocol("WM_DELETE_WINDOW",self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.distancematrix=transforms.filterNet(m,keeptheserows)

        self.ok()

    def applyme(self):

        self.result=self.distancematrix


class HimmeliWaiter(GenericLoadWaiter):
    """Generic class for wait windows displaying one item"""
    

    def __init__(self,parent,net,time,usemst,coordinates=None,title=None,titlemsg="Please wait",specmsg="Processing..."):

        Toplevel.__init__(self,parent)
        #self.configure(bg='Gray80')
        self.transient(parent)

        self.title(title)

        self.titlemsg=titlemsg
        self.specmsg=specmsg

        self.parent=parent
        self.result=None

        body=Frame(self)
        self.initial_focus=self.body(self,body)
        body.pack(padx=5,pady=5)

       # self.buttonbox()
        self.grab_set()

        self.points=['.','..','...']
        self.pointcounter=0

        self.displayBusyCursor()
        
        if not self.initial_focus:
            self.initial_focus(self)

        self.protocol("WM_DELETE_WINDOW",self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        if usemst:
            h=visuals.Himmeli(net,time=time,coordinates=coordinates,useMST=usemst)
        else:
            h=visuals.Himmeli(net,time=time,useMST=usemst)

        
        self.coords=h.getCoordinates()
        del h

        self.ok()

    def applyme(self):

        self.result=self.coords


class Meter(Frame):
    '''A simple progress bar widget.
    Made by Michael'''
    def __init__(self, master, fillcolor='orchid1', text='',
                 value=0.0, **kw):
        Frame.__init__(self, master, bg='white', width=350,
                               height=20)
        self.configure(**kw)
        
        self._c = Canvas(self, bg=self['bg'],
                         width=self['width'], height=self['height'],\
                             highlightthickness=0, relief='flat',
                         bd=0)
        self._c.pack(fill='x', expand=1)
        self._r = self._c.create_rectangle(0, 0, 0,
                                           int(self['height']), fill=fillcolor, width=0)
        self._t = self._c.create_text(int(self['width'])/2,
                                      int(self['height'])/2, text='')

        self.set(value, text)

    def set(self, value=0.0, text=None):
        #make the value failsafe:
        if value < 0.0:
            value = 0.0
        elif value > 1.0:
            value = 1.0
        if text == None:
            #if no text is specified get the default percentage string:
            text = str(int(round(100 * value))) + ' %'
        self._c.coords(self._r, 0, 0, int(self['width']) * value,
                       int(self['height']))
        self._c.itemconfigure(self._t, text=text)
        self.update()



class BootstrapPopulationsDialog(MySimpleDialog):
    """Asks for the number of bins for log binning
       and allows linear bins for 1...10"""


    def __init__(self,parent,title='Please provide information:',initRuns=100):

        Toplevel.__init__(self,parent)
        
        self.title(title)

        self.parent=parent
        self.result=None
        self.initRuns=initRuns
        self.nRuns=StringVar()

        body=Frame(self)
        self.initial_focus=self.body(self,body)
        body.pack(padx=5,pady=5)

        self.buttonbox()
        self.grab_set()

        if not self.initial_focus:
            self.initial_focus(self)

        self.protocol("WM_DELETE_WINDOW",self.cancel)
        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,parent.winfo_rooty()+50))
        self.initial_focus.focus_set()
        self.wait_window(self)

    def body(self,masterclass,masterwindow):

        self.c1=Label(masterwindow,text='Number of repetitions?',
                      anchor=W)
        self.c1.grid(row=0,column=0)
        self.e1=Entry(masterwindow,textvariable=self.nRuns,bg='Gray95')
        self.nRuns.set(str(self.initRuns))
        self.e1.grid(row=1,column=0)

        self.c2=Label(masterwindow,text='Percentage of nodes in each location?',
                      anchor=W)
        self.c2.grid(row=2,column=0)        
        self.scale_1=Scale(masterwindow, from_=0.0, to=1.0, 
                           orient=HORIZONTAL,resolution=0.01)
        self.scale_1.set(0.5)
        self.scale_1.grid(row=3,column=0)

        return self.c1

    def applyme(self):
        self.result=[self.scale_1.get(),self.e1.get()]

class SliderDialog2(MySimpleDialog):
    """Asks for the number of bins for log binning
       and allows linear bins for 1...10"""


    def __init__(self,parent,sMin,sMax,sRes,sStart,title='Please provide information:',bodyText="How much?"):

        Toplevel.__init__(self,parent)
        
        self.title(title)

        self.sMin=sMin
        self.sMax=sMax
        self.sRes=sRes
        self.sStart=sStart
        self.bodyText=bodyText
        self.parent=parent
        self.result=None
        self.clones=StringVar()
        self.measuretype=StringVar()

        body=Frame(self)
        self.initial_focus=self.body(self,body)
        body.pack(padx=5,pady=5)

        self.buttonbox()
        self.grab_set()

        if not self.initial_focus:
            self.initial_focus(self)

        self.protocol("WM_DELETE_WINDOW",self.cancel)
        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,parent.winfo_rooty()+50))
        self.initial_focus.focus_set()
        self.wait_window(self)

    def body(self,masterclass,masterwindow):
        self.c=Label(masterwindow,text='Percentage of nodes in each location?',
                      anchor=W)
        self.c.grid(row=0,column=0)        
        self.scale=Scale(masterwindow, from_=self.sMin, to=self.sMax, 
                           orient=HORIZONTAL,resolution=self.sRes)
        self.scale.set(self.sStart)
        self.scale.grid(row=1,column=0)

        return self.c

    def applyme(self):
        self.result=self.scale.get()


class AskNumberDialog(MySimpleDialog):
    """Asks for a number 
    """


    def __init__(self,parent,title='Please provide information:',bodyText="Number: ",initNumber=100):

        Toplevel.__init__(self,parent)
        
        self.title(title)

        self.bodyText=bodyText
        self.parent=parent
        self.result=None
        self.initNumber=initNumber
        self.theNumber=StringVar()
        #self.measuretype=StringVar()

        body=Frame(self)
        self.initial_focus=self.body(self,body)
        body.pack(padx=5,pady=5)

        self.buttonbox()
        self.grab_set()

        if not self.initial_focus:
            self.initial_focus(self)

        self.protocol("WM_DELETE_WINDOW",self.cancel)
        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,parent.winfo_rooty()+50))
        self.initial_focus.focus_set()
        self.wait_window(self)

    def body(self,masterclass,masterwindow):

        self.c1=Label(masterwindow,text=self.bodyText,
                      anchor=W)
        self.c1.grid(row=0,column=0)
        self.e1=Entry(masterwindow,textvariable=self.theNumber,bg='Gray95')
        self.theNumber.set(str(self.initNumber))
        self.e1.grid(row=1,column=0)

        return self.c1

    def applyme(self):
        self.result=self.e1.get()


class EDENRadioDialog(MySimpleDialog):
    """A dialog with a radio selector."""
    
    def __init__(self,parent,title=None,titlemsg="Choose one:",options=["Default"]):
        Toplevel.__init__(self,parent)
        self.title(title)
        self.titlemsg=titlemsg
        self.parent=parent
        self.options=options

        self.result=StringVar()
      
        body=Frame(self)
        self.initial_focus=self.body(self,body)
        body.pack(padx=5,pady=5)

        self.buttonbox()
        self.grab_set()

        if not self.initial_focus:
            self.initial_focus(self)

        self.protocol("WM_DELETE_WINDOW",self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.wait_window(self)

    def body(self,masterclass,masterwindow,titlemsg="Choose genetic distance measure"):
        self.wholeframe=Frame(masterwindow,relief='sunken',borderwidth=2)
        
        self.clabel=Label(self.wholeframe,text=self.titlemsg,justify=LEFT,anchor=W,bg='gray90',relief='groove',borderwidth=1)
        self.clabel.pack(side=TOP,expand=YES,fill=X,ipadx=5,ipady=5)

        self.bottompart=Frame(self.wholeframe)

        buttons=[]
        for i,value in enumerate(masterclass.options):
            r=Radiobutton(self.bottompart,text=value,value=value,variable=masterclass.result)
            r.grid(row=i+1,column=0,sticky=W)
        masterclass.result.set(masterclass.options[0])

        self.bottompart.pack(side=TOP,expand=YES,fill=BOTH,ipadx=7,ipady=7)
        self.wholeframe.pack(side=TOP,expand=YES,fill=BOTH)
       
        return self.wholeframe

    def applyme(self):
        self.result=(self.result.get())
        
