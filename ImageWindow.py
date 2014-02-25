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
import matplotlib,numpy

from DataWindow import DataWindow

class Myplot(object):
    '''empty container object'''
    pass


class PlotContainer:
    def __init__(self,figure):
        self.figure=figure


# -------  NETWORK VISUALIZATION WINDOW -------


class ImageWindow(Toplevel):
    '''Visualizes the network'''    
    def __init__(self,parent,net,namestring=None,sizeX=600,sizeY=600,useMST=True,coords=[],time=20,edgeColorMap='primary',parentnet=None):
        Toplevel.__init__(self,parent,bg='gray90')
        if namestring:
            self.title(namestring)
        mBar=Frame(self,relief='raised',borderwidth=2,bg='steelblue')

        self.net=net
        self.parentnet=parentnet
        self.parent=parent
        self.time=time
        self.useMST=useMST
        usematplotlib=True 

        #Calculate the betweenness centrality
        bcWaiter=dialogues.WaitWindow(self,title="Processsing",titlemsg="Calculating betweenness centrality",hasProgressbar=True)
        bcFunction=lambda:netext.getBetweennessCentrality(net)
        bcWaiter.go(bcFunction)
        bc=bcWaiter.result
        netext.addNodeProperty(net,"betweenness")
        for node in net:
            net.nodeProperty["betweenness"][node]=bc[node]

        # first generate colors for each node, BY STRENGTH
        viz=PlotParams()
        viz.nodeColors=visuals.getNodeColors(net,colorwith='betweenness')
        viz.equalsize=False
        viz.nodeSizes=visuals.getNodeSizes(net,size_by="betweenness",minsize=viz.nodeSize_min,maxsize=viz.nodeSize_max)
        self.bgimage=None

        #Save plot parameters to this object
        self.viz=viz

        #Coordinates for the nodes:
        #we should also check that x and y are numeric
        if hasattr(net,"nodeProperty") and "x" in net.nodeProperty and "y" in net.nodeProperty:
            allowedCoordProperties=["int","float","number"]
            propertyTypes=netext.getPropertyTypes(net)
            if propertyTypes["x"] in allowedCoordProperties and propertyTypes["y"] in allowedCoordProperties:
                self.hasPropertyCoords=True
                self.xy=self.propertiesToCoords("x","y")
            else:
                self.hasPropertyCoords=False
        else:
            self.hasPropertyCoords=False

        if not self.hasPropertyCoords:
            self.xy=self.getHimmeliCoords(parent,net,time,useMST,coords)
        xy=self.xy

        # top 10 strength nodes are shown if "top 10" ticked
        plotbody=Frame(self,bg='gray90')
        self.plotbody=plotbody

        myplot=Myplot()
        bgfig=matplotlib.pyplot.figure()
        bgc=PlotContainer(bgfig)
        myplot.canvas=FigureCanvasTkAgg(bgfig,master=plotbody)
          
        #-> drawing code here        
        myplot.thisfigure=visuals.VisualizeNet(net,self.xy,nodeColors=viz.nodeColors,nodeColorMap=viz.nodeColorMap,equalsize=viz.equalsize,minnode=viz.nodeSize_min,maxnode=viz.nodeSize_max,nodeSize=viz.nodeSize,uselabels=viz.uselabels,labels=viz.nodeLabels,nodeSizes=viz.nodeSizes,showAllNodes=viz.showAllNodes,minwidth=viz.minwidth,maxwidth=viz.maxwidth,fontsize=viz.fontsize,edgeColorMap=viz.edgeColorMap,bgcolor=viz.bgcolor,interactive=True,keepAspectRatio=True,baseFig=bgc)
        myplot.thisfigure.subplots_adjust(left=0,right=1,top=1,bottom=0)
        self.thisfigure=myplot.thisfigure
        #<- drawing code ends


        myplot.canvas.show()

        myplot.toolbar=NavigationToolbar2TkAgg(myplot.canvas,plotbody)
        myplot.toolbar.update()
        myplot.canvas._tkcanvas.pack(side=TOP,fill=BOTH,expand=YES)

        FileBtn=self.makeFileMenu(mBar,myplot,net,xy,viz,namestring)
        OptBtn=self.makeOptionsMenu(mBar,self.parent,myplot,net,xy,viz,namestring,self.parentnet)
        mBar.tk_menuBar=(FileBtn,OptBtn)
        mBar.pack(side=TOP,fill=X,expand=YES)
        self.protocol("WM_DELETE_WINDOW",self.cancel)
        plotbody.pack(side=BOTTOM,fill=BOTH,expand=YES)

        myplot.thisfigure.startInteraction()

    def makeFileMenu(self,mBar,myplot,net,xy,viz,namestring='temp'):
        '''Defines the File menu of the data window'''
        
        FileBtn=Menubutton(mBar,text="File",underline=0,bg='steelblue')
        FileBtn.pack(side=LEFT,padx="2m")
        FileBtn.menu=Menu(FileBtn)
        FileBtn.menu.add_command(label="Save figure",underline=0,command=lambda s=self, n=namestring: s.save_figure(n))
        FileBtn.menu.add_command(label="Save display coords",underline=0,command=lambda s=self, n=net, x=xy : s.save_coords(n,x))
        FileBtn.menu.add_command(label="Close",underline=0,command=lambda s=self: s.close_window())
        FileBtn['menu']=FileBtn.menu
        return FileBtn

    def makeOptionsMenu(self,mBar,root,myplot,net,xy,viz,namestring='temp',parentnet=[]):
        '''Defines the Options menu'''

        OptBtn=Menubutton(mBar,text="Options",underline=0,bg='steelblue')
        OptBtn.pack(side=LEFT,padx="2m")
        OptBtn.menu=Menu(OptBtn)

        OptBtn.menu.nodecolors=Menu(OptBtn.menu)
        OptBtn.menu.nodecolors.add_command(label='None',underline=0,command=lambda s=self,n=net,x=xy,m=myplot,v=viz,p=parentnet: s.change_nodeColors(n,x,p,m,v,'white'))

        bc_calculated=False

        # ----------------- menu for changing node colors ----------------
        if hasattr(net,'nodeProperty'):
            props=list(net.nodeProperty)
            for prop in props:

                OptBtn.menu.nodecolors.add_command(label='By '+prop,underline=0,command=lambda s=self,n=net,x=xy,m=myplot,v=viz,p=parentnet,c=prop: s.change_nodeColors(n,x,p,m,v,c))


        # ----------------- menu for changing node labels   
        OptBtn.menu.labeltext=Menu(OptBtn.menu)
        OptBtn.menu.labeltext.add_command(label='None',underline=0,command=lambda s=self,n=net,x=xy,m=myplot,v=viz,p=parentnet: s.change_nodeLabels(n,x,p,m,v,'none'))
        OptBtn.menu.labeltext.add_command(label='Node name',underline=0,command=lambda s=self,n=net,x=xy,m=myplot,v=viz,p=parentnet: s.change_nodeLabels(n,x,p,m,v,'index'))

        if hasattr(net,'nodeProperty'):

            props=list(net.nodeProperty)

            for prop in props:

                OptBtn.menu.labeltext.add_command(label=prop,underline=0,command=lambda s=self,n=net,x=xy,m=myplot,v=viz,p=parentnet,c=prop: s.change_nodeLabels(n,x,p,m,v,c))

  
        # ----------- menu for changing node sizes --------------
        OptBtn.menu.nodesize=Menu(OptBtn.menu)
        OptBtn.menu.nodesize.add_command(label='Fixed size',underline=0,command=lambda s=self,n=net,x=xy,m=myplot,v=viz,p=parentnet: s.change_nodeSize(n,x,p,m,v,'fixed'))

        if hasattr(net,'nodeProperty'):
            props=netext.getNumericProperties(net)

            for prop in props:
                OptBtn.menu.nodesize.add_command(label=prop,underline=0,command=lambda s=self,n=net,x=xy,m=myplot,v=viz,p=parentnet,c=prop: s.change_nodeSize(n,x,p,m,v,c))

        OptBtn.menu.nodesizemenu=Menu(OptBtn.menu)
        OptBtn.menu.nodesizemenu.add_cascade(label='Size by',menu=OptBtn.menu.nodesize)
        OptBtn.menu.nodesizemenu.add_command(label='Adjust scale',underline=0,command=lambda s=self,n=net,x=xy,m=myplot,v=viz,p=parentnet,r=root: s.change_nodeScale(n,x,p,m,v,r))
                                            
        OptBtn.menu.nodelabels=Menu(OptBtn.menu)
        OptBtn.menu.nodelabels.add_cascade(label='Label text',menu=OptBtn.menu.labeltext)
        OptBtn.menu.nodelabels.add_command(label='Font size',underline=0,command=lambda s=self,n=net,x=xy,m=myplot,v=viz,p=parentnet,r=root: s.change_nodeLabelFontSize(n,x,p,m,v,r))

        OptBtn.menu.nodecolor_submenu=Menu(OptBtn.menu)
        OptBtn.menu.nodecolor_submenu.add_cascade(label='Color by',menu=OptBtn.menu.nodecolors)
        OptBtn.menu.nodecolor_submenu.add_command(label='Color map',underline=0,command=lambda s=self,n=net,x=xy,m=myplot,v=viz,p=parentnet,r=root: s.change_nodeColorMap(n,x,p,m,v,r))

        OptBtn.menu.nodes=Menu(OptBtn.menu)
        OptBtn.menu.nodes.add_cascade(label='Color',menu=OptBtn.menu.nodecolor_submenu)
        OptBtn.menu.nodes.add_cascade(label='Size',menu=OptBtn.menu.nodesizemenu)
        OptBtn.menu.nodes.add_cascade(label='Labels',menu=OptBtn.menu.nodelabels)

        OptBtn.menu.edges=Menu(OptBtn.menu)
        OptBtn.menu.edges.add_command(label='Color map',underline=0,command=lambda s=self,n=net,x=xy,m=myplot,v=viz,p=parentnet,r=root: s.change_edgeColorMap(n,x,p,m,v,r))
        OptBtn.menu.edges.add_command(label='Edge width',underline=0,command=lambda s=self,n=net,x=xy,m=myplot,v=viz,p=parentnet,r=root: s.change_edgeWidth(n,x,p,m,v,r))

        OptBtn.menu.background=Menu(OptBtn.menu)
        OptBtn.menu.background.add_command(label='Black',underline=0,command=lambda s=self,n=net,x=xy,m=myplot,v=viz,p=parentnet: s.change_background(n,x,p,m,v,'black'))
        OptBtn.menu.background.add_command(label='White',underline=0,command=lambda s=self,n=net,x=xy,m=myplot,v=viz,p=parentnet: s.change_background(n,x,p,m,v,'white'))
        OptBtn.menu.background.add_command(label='Image...',underline=0,command=lambda s=self,n=net,x=xy,m=myplot,v=viz,p=parentnet: s.change_background_image(n,x,p,m,v))

        OptBtn.menu.layout=Menu(OptBtn.menu)
        OptBtn.menu.layout.add_command(label='Calculate new layout',underline=0,command=lambda s=self,n=net,x=xy,m=myplot,v=viz,p=parentnet: s.change_layout(n,x,p,m,v,'new'))
        OptBtn.menu.layout.add_command(label='Load coordinate file',underline=0,command=lambda s=self,n=net,x=xy,m=myplot,v=viz,p=parentnet: s.change_layout(n,x,p,m,v,'load'))
        OptBtn.menu.layout.add_command(label='Use auxiliary data (x and y)',underline=0,command=lambda s=self,n=net,x=xy,m=myplot,v=viz,p=parentnet: s.change_layout(n,x,p,m,v,'aux'))

        if not self.hasPropertyCoords:
            OptBtn.menu.layout.entryconfig(3, state=DISABLED)

        OptBtn.menu.add_cascade(label='Nodes',menu=OptBtn.menu.nodes)
        OptBtn.menu.add_cascade(label='Edges',menu=OptBtn.menu.edges)
        OptBtn.menu.add_cascade(label='Background',menu=OptBtn.menu.background)
        OptBtn.menu.add_cascade(label='Layout',menu=OptBtn.menu.layout)

        OptBtn['menu']=OptBtn.menu
        return OptBtn

    def change_layout(self,network,xy,parentnet,myplot,viz,command):
        if command=="new":
            self.xy=self.getHimmeliCoords(self.parent,network,self.time,False,None)
        elif command=="aux":
            self.xy=self.propertiesToCoords("x","y")
        else:
            tempxy=self.load_coords()
            if tempxy!=None:
                self.xy=tempxy
        self.redraw_Network(network,None,parentnet,myplot,viz)

    def change_background(self,network,xy,parentnet,myplot,viz,command):
        self.bgimage=None
        viz.bgcolor=command
        self.redraw_Network(network,xy,parentnet,myplot,viz)

    def propertiesToCoords(self,xProperty,yProperty):
        """
        Uses xProperty and yProperty named node properties to
        return a coordinate map to be used with the visualization functions.
        """
        coords={}
        for node in self.net:
            x=self.net.nodeProperty[xProperty][node]
            y=self.net.nodeProperty[yProperty][node]
            coords[node]=(x,y)
        return coords

    def getHimmeliCoords(self,parent,net,time,useMST,coords):
        if useMST==True:
            temp=dialogues.HimmeliWaiter(parent,net,time,useMST,title="Visualizing network",coordinates=coords)
        else:
            temp=dialogues.HimmeliWaiter(parent,net,time,useMST,title="Visualizing network")
        return temp.coords

    def get_Betweenness(self,network,xy,parentnet,myplot,viz,command="size"):
        bc=netext.getBetweennessCentrality(network)
        netext.addNodeProperty(network,"betweenness")

        for node in network:
            network.nodeProperty["betweenness"][node]=bc[node]

        if command=="color":
            self.change_nodeColors(network,xy,parentnet,myplot,viz,"betweenness")
        elif command=="size":
            self.change_nodeSize(network,xy,parentnet,myplot,viz,"betweenness")


    def change_edgeWidth(self,network,xy,parentnet,myplot,viz,root):

        newsizes_get=dialogues.DoubleSliderDialog(root,title='Select width limits',titlemsg='Select width limits',resolution=0.2,minval_1=0.2,maxval_1=12,currval_1=viz.minwidth,minval_2=0.2,maxval_2=12,currval_2=viz.maxwidth)
        newsizes=newsizes_get.result
        if newsizes!=None:
            viz.minwidth=newsizes[0]
            viz.maxwidth=newsizes[1]
        
            self.redraw_Network(network,xy,parentnet,myplot,viz)

    def change_edgeColorMap(self,network,xy,parentnet,myplot,viz,root):
        newmap=dialogues.ColorMapDialog(root,viz.nodeColorMap)
        if newmap.result!=None:
            viz.edgeColorMap=newmap.result
            self.redraw_Network(network,xy,parentnet,myplot,viz)    

    def change_nodeColorMap(self,network,xy,parentnet,myplot,viz,root):
        newmap=dialogues.ColorMapDialog(root,viz.nodeColorMap)
        if newmap.result!=None:
            viz.nodeColorMap=newmap.result
            if not(viz.nodeColorData=='white'):

                viz.nodeColors=visuals.getNodeColors(network,colorwith=viz.nodeColorData,parentnet=parentnet,useColorMap=viz.nodeColorMap)      

            self.redraw_Network(network,xy,parentnet,myplot,viz)

    def change_nodeSize(self,network,xy,parentnet,myplot,viz,command):
        if command=="fixed":
            viz.equalsize=True
            viz.nodeSizeData="fixed"
        else:
            viz.nodeSizeData=command
            viz.nodeSizes=visuals.getNodeSizes(network,size_by=command,minsize=viz.nodeSize_min,maxsize=viz.nodeSize_max)
            viz.equalsize=False
        
        self.redraw_Network(network,xy,parentnet,myplot,viz)
                                             
    def change_nodeScale(self,network,xy,parentnet,myplot,viz,root):
        if viz.equalsize==True:
            newsize=dialogues.SliderDialog(root,title='Select node size',titlemsg='Select node size',minval=1,maxval=35,currval=viz.nodeSize)
            if newsize.result!=None:
                viz.nodeSize=newsize.result
        else:
            newsizes_get=dialogues.DoubleSliderDialog(root,title='Select size limits',titlemsg='Select size limits',minval_1=1,maxval_1=24,currval_1=viz.nodeSize_min,minval_2=1,maxval_2=24,currval_2=viz.nodeSize_max)
            if newsizes_get.result!=None:
                newsizes=newsizes_get.result
                viz.nodeSize_min=newsizes[0]
                viz.nodeSize_max=newsizes[1]
                viz.nodeSizes=visuals.getNodeSizes(network,size_by=viz.nodeSizeData,minsize=viz.nodeSize_min,maxsize=viz.nodeSize_max)
            
        self.redraw_Network(network,xy,parentnet,myplot,viz)

    def change_nodeLabelFontSize(self,network,xy,parentnet,myplot,viz,root):
        newsize=dialogues.SliderDialog(root,title='Select label font size',titlemsg='Select label font size',minval=4,maxval=24,currval=viz.fontsize)

        if newsize.result!=None:
            viz.fontsize=newsize.result
            self.redraw_Network(network,xy,parentnet,myplot,viz)

    def change_nodeColors(self,network,xy,parentnet,myplot,viz,command):
        viz.nodeColorData=command
        if command=='white':
            for node in viz.nodeColors.keys():
                viz.nodeColors[node]=(1,1,1)

        elif command=='strength':
            viz.nodeColors=visuals.getNodeColors(network,colorwith='strength',parentnet=[],useColorMap=viz.nodeColorMap)            
        else:
            if command=='betweenness':
                parentnet=[]          
            viz.nodeColors=visuals.getNodeColors(network,colorwith=command,parentnet=parentnet,useColorMap=viz.nodeColorMap)

        self.redraw_Network(network,xy,parentnet,myplot,viz)

    def change_nodeLabels(self,network,xy,parentnet,myplot,viz,command):
        '''Changes node labels in the property object viz. Again in case parentnet is used,
        i.e. a network containing more nodes than the network itself, the label dict is
        created for all of these.'''

        if command=='none':
            viz.uselabels='none'
        else:
            viz.uselabels='custom'
            if command=="betweenness":
                parentnet=[]

            if len(parentnet)>0:
                thisnet=parentnet
            else:
                thisnet=network
        
            if command=='index':
                for node in thisnet:
                    viz.nodeLabels[node]=str(node)
            elif hasattr(thisnet,'nodeProperty'):
                    plist=list(thisnet.nodeProperty)
                    if command in plist:
                        for node in thisnet:
                            viz.nodeLabels[node]=thisnet.nodeProperty[command][node]

        self.redraw_Network(network,xy,parentnet,myplot,viz)

    def redraw_Network(self,network,xy,parent,myplot,viz):
        if self.bgimage!=None:
            figsize=self.bgimage.shape[1],self.bgimage.shape[0]
        else:
            figsize=(1,1)

        myplot.thisfigure.stopInteraction()
        myplot.thisfigure.clf()
        myplot.thisfigure=visuals.VisualizeNet(network,self.xy,nodeColors=viz.nodeColors,nodeColorMap=viz.nodeColorMap,equalsize=viz.equalsize,minnode=viz.nodeSize_min,maxnode=viz.nodeSize_max,nodeSize=viz.nodeSize,uselabels=viz.uselabels,labels=viz.nodeLabels,nodeSizes=viz.nodeSizes,showAllNodes=viz.showAllNodes,minwidth=viz.minwidth,maxwidth=viz.maxwidth,fontsize=viz.fontsize,edgeColorMap=viz.edgeColorMap,bgcolor=viz.bgcolor,baseFig=myplot.canvas,interactive=True,keepAspectRatio=True,figsize=figsize)

        if self.bgimage!=None:
            ax=myplot.thisfigure.get_axes()[0];x=ax.get_xlim();y=ax.get_ylim()
            matplotlib.pyplot.imshow(self.bgimage,figure=myplot.canvas.figure,origin=self.bgimage_origin,extent=[x[0],x[1],y[0],y[1]])

        myplot.thisfigure.subplots_adjust(left=0,right=1,top=1,bottom=0)

        myplot.canvas.show()

        myplot.thisfigure.startInteraction()

    def change_background_image(self,network,xy,parentnet,myplot,viz,):
        filetypes=[("PNG file","*.png"),("JPG file","*.jpg")]
        types=["png","jpg"]
        bgimage_cand=tkFileDialog.askopenfilename(title='Open...',filetypes=filetypes,defaultextension=".png")
        if len(bgimage_cand)>0:
            ending=bgimage_cand.split(".")[-1]
            if ending.lower() in types:  
                #there is a bug in matplotlib imread when using pil to read it,
                #i.e. when the image is something else than png.
                self.bgimage=numpy.asarray(matplotlib.pyplot.imread(bgimage_cand))
                if ending.lower()=="png":
                    self.bgimage_origin="upper"
                else:
                    self.bgimage_origin="lower"
                self.redraw_Network(network,xy,parentnet,myplot,viz)
            else:
                tkMessageBox.showerror(
                    "Error loading image.",
                    "Unknown file format."
                    )

    def close_window(self):
        self.destroy()
        
    def save_figure(self,namestring):
        namestring=namestring.replace(" ","_")
        filetypes=[("PNG file","*.png"),("PDF file","*.pdf"),('EPS file','*.eps'),("SVG file","*.svg")]
        types=["png","pdf","eps","svg"]
        newname=tkFileDialog.asksaveasfilename(initialfile=namestring+".png",title='Save as',filetypes=filetypes,defaultextension=".png")
        if len(newname)>0: #if user did not press cancel button
            ending=newname.split(".")[-1]
            if ending in types:
                self.thisfigure.savefig(newname,facecolor=self.viz.bgcolor)
            else:
                tkMessageBox.showerror(
                    "Error saving the figure",
                    "Unknown file format. Use png, pdf, eps or svg."
                    )

    def save_png_figure(self,namestring):
        newname=tkFileDialog.asksaveasfilename(initialfile=namestring,title='Save as')
        if len(newname)>0:
            self.thisfigure.savefig(newname)

    def save_jpeg_file(self,namestring):
        filetypes=[('JPG file','*.jpg')]
        newfilename=namestring+'.jpg'
        newname=tkFileDialog.asksaveasfilename(filetypes=filetypes,initialfile=newfilename,title='Save as')
        if len(newname)>0:
            self.img2.save(newname)

    def load_coords(self):
        filetypes=[('Text file','*.txt')]
        filename=tkFileDialog.askopenfilename(title='Load coordinates',filetypes=filetypes)
        coords=None
        if len(filename)>0: #if user did not press cancel button
            f=open(filename,'rU')
            firstFields=f.readline().split()
            try:
                if firstFields!=["node_label","x","y"]:
                    raise Exception("First line should be:\n node_label x y")
                coords={}
                for i,line in enumerate(f):
                    fields=line.split()
                    if len(fields)!=3:
                        raise Exception("Syntax error in line: "+str(i+1))

                    #we only save coords for nodes which are in the net
                    if fields[0] in self.net:
                        coords[fields[0]]=(float(fields[1]),float(fields[2]))                        
                    elif fields[0].isdigit() and int(fields[0]) in self.net: #node names can be numerical, we'll check for that also
                        coords[int(fields[0])]=(float(fields[1]),float(fields[2]))                        
                    
                for node in self.net:
                    if node not in coords:
                        raise Exception("Node "+str(node)+" not found in the file.")
            except Exception,e:
                tkMessageBox.showerror(
                    "Invalid file format",
                    str(e)
                    )
        return coords

    def save_coords(self,net,xy):
        newname=tkFileDialog.asksaveasfilename(initialfile='',title='Save as')
        if len(newname)>0:
            fobj=open(newname,'w')
            fobj.write("node_label x y\n")
            for node in net:
                fobj.write(str(node)+' %f %f\n' % (self.xy[node][0],self.xy[node][1]))
            fobj=close()

    def __del__(self):
        if os.path.isfile(self.filename):
            os.remove(self.filename)

    def cancel(self,event=None):
        self.destroy()



class PlotParams(object):
    '''Container object for visualization parameters'''
    def __init__(self):
        self.nodeColors={}
        self.nodeColorData="strength"
        self.nodeLabels={}
        self.nodeLabelData="labels"
        self.bgcolor='black'
        self.nodeSizes={}
        self.nodeSizeData="strength"
        self.nodeSize=3.0
        self.nodeSize_min=2.0
        self.nodeSize_max=10.0
        self.nodeColorMap='orange'
        self.edgeColorMap='winter'
        self.fontsize=7
        self.uselabels='all'
        self.frame='false'
        self.maxwidth=2.0
        self.minwidth=0.4
        self.equalsize=False
        self.showAllNodes=True
