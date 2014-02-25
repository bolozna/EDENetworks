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
from __future__ import with_statement
# -*- coding: utf-8 -*-
from Tkinter import *
from netpython import *
import tkFileDialog,tkMessageBox
from pylab import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2TkAgg
from PIL import Image
from PIL import ImageTk
import os
import shutil,math
import random

# ----------------------------------------------------------------------------------- DATA WINDOW ----------------
        
class DataWindow(Toplevel):
    '''Plots input data using matplotlib. Inputs: root window,
       data in the form [[xdata1,ydata1],[xdata2,ydata2],...],
       plotcommand='plot|loglog', titlestring,xstring=x-axis label,
       ystring=y-axis label. If several data vectors, use ystring='ylabel1:ylabel2:...'
       '''
    def __init__(self,parent,data,namestring=None,plotcommand='plot',titlestring='',xstring='',ystring='',fontsize=12,addstring=''):
        Toplevel.__init__(self,parent,bg='gray90')
        if namestring:
            self.title(namestring)

        self.data=data      

        plotbody=Frame(self,bg='gray90')

        # create plot object using matplotlib
        myplot=visuals.ReturnPlotObject(data,plotcommand,titlestring,xstring,ystring,fontsize=fontsize,addstr=addstring)
        self.thisfigure=myplot.thisFigure

        # create a canvas out of the plot object, into the frame "plotbody"
        myplot.canvas=FigureCanvasTkAgg(myplot.thisFigure,master=plotbody)
        myplot.canvas.show()
        myplot.canvas._tkcanvas.pack(side=TOP,fill=BOTH,expand=YES)

        # put a toolbar below the plot canvas
        myplot.toolbar=NavigationToolbar2TkAgg(myplot.canvas,plotbody)  # these display plot navigation & save menu bar
        myplot.toolbar.update()

        # create menubar        
        mBar=Frame(self,relief='raised',borderwidth=2,bg='steelblue')

        FileBtn=self.makeFileMenu(mBar,self.data,namestring)
        PlotBtn=self.makePlotMenu(mBar,myplot)
        mBar.tk_menuBar=(FileBtn,PlotBtn)

        # pack both things: menubar is on top, frame plotbody containing plot canvas is below
        mBar.pack(side=TOP,fill=X,expand=YES)           
        plotbody.pack(side=BOTTOM,fill=X,expand=YES)

    def makeFileMenu(self,mBar,data,namestring):
        '''Defines the File menu of the data window'''        
        FileBtn=Menubutton(mBar,text="File",underline=0,bg='steelblue')
        FileBtn.pack(side=LEFT,padx="2m")
        FileBtn.menu=Menu(FileBtn)
        FileBtn.menu.add_command(label="Save graph",underline=0,command=lambda s=self,n=namestring: s.save_plot(n))
        FileBtn.menu.add_command(label="Export data",underline=0,command=lambda s=self, d=data: s.export_data(d))
        FileBtn.menu.add('separator')
        FileBtn.menu.add_command(label="Close",underline=0,command=self.close_window)
        FileBtn['menu']=FileBtn.menu
        return FileBtn

    def makePlotMenu(self,mBar,myplot):
        '''Defines the Axes menu of the data window'''        
        PlotBtn=Menubutton(mBar,text="Axes",underline=0,bg='steelblue')
        PlotBtn.pack(side=LEFT,padx="2m")
        PlotBtn.menu=Menu(PlotBtn)
        PlotBtn.menu.add_radiobutton(label="Linear",underline=0,command=lambda s=self, a=myplot.axes, c=myplot.canvas: s.setAxes(a,c,'linear','linear'))
        PlotBtn.menu.add_radiobutton(label="Double logarithmic",underline=0,command=lambda s=self, a=myplot.axes, c=myplot.canvas: s.setAxes(a,c,'log','log'))
        PlotBtn.menu.add_radiobutton(label="Linear x, log y",underline=0,command=lambda s=self, a=myplot.axes, c=myplot.canvas: s.setAxes(a,c,'linear','log'))
        PlotBtn.menu.add_radiobutton(label="Log x, linear y",underline=0,command=lambda s=self, a=myplot.axes, c=myplot.canvas: s.setAxes(a,c,'log','linear'))

        PlotBtn['menu']=PlotBtn.menu
        return PlotBtn

    def setAxes(self,a,c,xstring,ystring):
        setp(a,'xscale',xstring,'yscale',ystring)
        c.show()  

    def save_plot(self,namestring):
        filetypes=[("PNG file","*.png"),("PDF file","*.pdf"),('EPS file','*.eps'),("SVG file","*.svg")]
        types=["png","pdf","eps","svg"]
        newname=tkFileDialog.asksaveasfilename(initialfile=namestring+".png",title='Save as',filetypes=filetypes,defaultextension=".png")
        if len(newname)>0: #if user did not press cancel button
            ending=newname.split(".")[-1]
            if ending in types:
                self.thisfigure.savefig(newname)
            else:
                tkMessageBox.showerror(
                    "Error saving the figure",
                    "Unknown file format. Use png, pdf, eps or svg."
                    )

    def export_data(self,data):
        '''Saves displayed data to text file'''
        filename=tkFileDialog.asksaveasfilename()
        if (len(data)>1):
            tkMessageBox.showinfo('Multiple data vectors will be saved as %s_1.%s etc.' % (filename.split('.')[0],filename.split('.')[1]))
            for i in range(len(data)):

                fname=filename.split('.')[0]+'_'+str(i)+'.'+filename.split('.')[1]
                try:
                    fobj=open(fname,'w')
                    for j in range(len(data[i][0])):
                        fobj.write('%f %f\n' % (data[i][0][j],data[i][1][j]))
                    fobj=close()
                except IOerror,e:
                    tkMessageBox.showerror(message='Error saving file') 

        else:
            try:
                fobj=open(filename,'w')    
                for i in range(len(data[0][0])):
                    fobj.write('%f %f\n' % (data[0][0][i],data[0][1][i]))
                fobj=close()
                tkMessageBox.showinfo(message='Save succesful!')

            except IOError,e:
                tkMessageBox.showerror(message='Error saving file') 

    def close_window(self):
        self.destroy()
