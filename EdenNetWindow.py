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
import pylab
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2TkAgg
from PIL import Image
from PIL import ImageTk
import os
import shutil,math
import random
from ImageWindow import ImageWindow
from ImageWindow import Myplot
from DataWindow import DataWindow

# ------- EDEN_NETWINDOW: MAIN NETWORK ANALYSIS WINDOW-------------
# -                                                               -
# - EdenNetWindow: displays network data & contains menu commands -
# -                                                               -
# -----------------------------------------------------------------

class EdenNetWindow(Toplevel):
    '''Class for displaying network data and related menus'''
    def __init__(self,parent,network,namestring='another network',coords=None,nettype=1,mstnet=None,parentnet=None):
        Toplevel.__init__(self,parent,bg='gray80')
        if namestring:
            self.title(namestring)

        self.parent=parent
        self.coords=coords
        self.nettype=nettype
        self.netname=namestring
        self.mstnet=mstnet
        self.parentnet=parentnet
        self.net=network

        #--- Calculating basic statistics
        N=len(self.net)
        E=len(self.net.edges)
        avgk=2.0*E/N

        self.minw=min(self.net.weights) #saved for dermining if log bins are to be used
        self.maxw=max(self.net.weights)
        if (self.minw==self.maxw):
            self.weighted=FALSE
        else:
            self.weighted=TRUE
            
        degs=netext.deg(network)
        maxk=max(degs.itervalues())
        mink=min(degs.itervalues())

        averageClustering=netanalysis.globalClustering(network)

        components=percolator.getComponents(network)
        Ncomponents=len(components)
        maxcomponent=components.getGiantSize()
        

        #--- Menubar
        mBar=Frame(self,relief='raised',borderwidth=2,bg='steelblue')
        FileBtn=self.makeFileMenu(mBar,self.net,namestring,self.parent,parentnet)
        AnalyzeBtn=self.makeAnalyzeMenu(mBar,self.net,namestring,self.parent,self.weighted)
        VisualizeBtn=self.makeVisualizeMenu(mBar,self.net,namestring,self.parent,self.coords,weighted=self.weighted,parentnet=self.parentnet)
        mBar.tk_menuBar=(FileBtn,AnalyzeBtn)

        # infoframe
        lowerhalf=Frame(self,bg='gray90')

        # ----------------------- network visualization frame ---------------
        netview=Frame(lowerhalf,bg='gray80',relief='raised',borderwidth=2)
     
        #Calculating the coordinates (and saving the Himmeli run?)
        if coords==None:
            h=visuals.Himmeli(self.net,time=20)
        else:
            h=visuals.Himmeli(self.net,coordinates=self.coords,time=20,useMST=True)
        filename=h.netName+"_0001.eps"
        if os.path.isfile(filename):
            os.remove(filename)
        self.himmeli=h 
        xy=self.himmeli.getCoordinates()

        labels={}
        myplot=Myplot()
        myplot.thisfigure=visuals.VisualizeNet(self.net,xy,edgeColorMap="primary",coloredvertices=False,labels=labels,equalsize=True,vsize=3,bgcolor='black',uselabels=False,figsize=(3,3))
        self.thisfigure=myplot.thisfigure
        myplot.canvas=FigureCanvasTkAgg(myplot.thisfigure,master=netview)
        myplot.canvas.show()
        myplot.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=YES, padx=10,pady=10) 

        lwidth=30

        body=Frame(lowerhalf,bg='gray80',relief='raised',borderwidth=2)
        
        Label(body,text="Network:",relief='flat',bg='DarkOliveGreen2',justify=RIGHT,width=lwidth).grid(row=0)
        Label(body,text="%s" % namestring,relief='flat',bg='DarkOliveGreen3',width=lwidth).grid(row=0,column=1)

        Label(body,text="Size:",relief='flat',bg='gray70',width=lwidth).grid(row=1)
        Label(body,text="N=%d nodes" % N,relief='flat',bg='gray60',width=lwidth).grid(row=1,column=1)

        Label(body,text="Edges:",relief='flat',bg='gray70',width=lwidth).grid(row=2)
        Label(body,text="E=%d edges" % E,relief='flat',bg='gray60',width=lwidth).grid(row=2,column=1)

        Label(body,text="Average degree:",relief='flat',bg='gray70',width=lwidth).grid(row=3)
        Label(body,text="<k>=%3.2f" % avgk,relief='flat',bg='gray60',width=lwidth).grid(row=3,column=1)

        Label(body,text="Max degree:",relief='flat',bg='gray70',width=lwidth).grid(row=4)
        Label(body,text="kmax=%d" % maxk,relief='flat',bg='gray60',width=lwidth).grid(row=4,column=1)

        Label(body,text="Average clustering:",relief='flat',bg='gray70',width=lwidth).grid(row=5)
        Label(body,text="<c>=%3.2f" % averageClustering,relief='flat',bg='gray60',width=lwidth).grid(row=5,column=1)

        # additional commands for weighted networks:

        if self.weighted:
            avgw=sum(self.net.weights)/float(E)

            if self.nettype==0:
                wstring1='distance'
                wstring2='d'
            else:
                wstring1='weight'
                wstring2='w'

            color1='gray65'
            color2='gray55'

            Label(body,text="Average %s:" % wstring1,relief='flat',bg=color1,width=lwidth).grid(row=6)
            Label(body,text="<%s>=%3.2f" % (wstring2,avgw),relief='flat',bg=color2,width=lwidth).grid(row=6,column=1)

            Label(body,text="Min and max %s:" % wstring1,relief='flat',bg=color1,width=lwidth).grid(row=7)
            Label(body,text="%smin=%3.2f, %smax=%3.2f" % (wstring2,self.minw,wstring2,self.maxw), relief='flat',bg=color2,width=lwidth).grid(row=7,column=1) 
	    rowcounter=8              
        else:
            Label(body,text="This network is unweighted.",relief='flat',justify=LEFT,bg='gray80',width=2*lwidth).grid(row=6,columnspan=2)
	    rowcounter=7

	rowcounter=rowcounter+1
        Label(body,text="Connectivity:",relief='flat',bg='gray70',width=lwidth).grid(row=rowcounter)

        if Ncomponents==1:
            Label(body,text="connected",relief='flat',bg='gray60',width=lwidth).grid(row=rowcounter,column=1)
        else:
            Label(body,text="%d components, largest %d" % (Ncomponents,maxcomponent),relief='flat',bg='gray60',width=lwidth).grid(row=rowcounter,column=1)

	rowcounter=rowcounter+1

	if hasattr(network,'nodeProperty'):

		Label(body,text="Imported attributes",relief='flat',bg='DarkOliveGreen2',width=lwidth).grid(row=rowcounter)
		Label(body,text="type",relief='flat',bg='DarkOliveGreen3',width=lwidth).grid(row=rowcounter,column=1)

		rowcounter=rowcounter+1

		propertydict=netext.getPropertyTypes(network)

		for prop in propertydict.keys():

			Label(body,text=prop,relief='flat',bg='gray70',width=lwidth).grid(row=rowcounter)
			Label(body,text=propertydict[prop],relief='flat',bg='gray65',width=lwidth).grid(row=rowcounter,column=1)
			rowcounter=rowcounter+1


        netview.grid(row=0,column=0,sticky=NW)
        body.grid(row=0,column=1,sticky=NW)
        
        mBar.pack(side=TOP,fill=X,expand=YES,anchor=NW)
        lowerhalf.pack(side=TOP,fill=BOTH,expand=YES,anchor=NW)

    def dummy(self):
        pass

    def makeFileMenu(self,mBar,network,namestring,root,parentnet):
        FileBtn=Menubutton(mBar,text="File",underline=0,bg='steelblue')
        FileBtn.pack(side=LEFT,padx="2m")
        FileBtn.menu=Menu(FileBtn)

        FileBtn.menu.add_command(label="Save network",underline=0,command=lambda s=self,net=network,n=namestring: s.save_network(net))
        FileBtn.menu.add_command(label="Export component data",underline=0,command=lambda s=self,net=network,n=namestring,p=parentnet: s.export_components(net,n,p))
        FileBtn.menu.add_command(label="Export node data",underline=0,command=lambda s=self,net=network,n=namestring,p=parentnet: s.export_nodes(net,n,p))

        FileBtn.menu.add('separator')
        FileBtn.menu.add_command(label="Close",underline=0,command=self.exit_network)
        FileBtn['menu']=FileBtn.menu
        return FileBtn

    def makeVisualizeMenu(self,mBar,network,namestring,root,coords,parentnet=None,weighted=FALSE):
        VisualizeBtn=Menubutton(mBar,text="View",underline=0,bg='steelblue')
        VisualizeBtn.pack(side=LEFT,padx="2m")
        VisualizeBtn.menu=Menu(VisualizeBtn)
        VisualizeBtn.menu.add_command(label="View network",underline=0,command=lambda s=self,x=network,c=coords,n=namestring,r=root,p=parentnet: s.draw_MST_network(x,c,n,r,p),state=NORMAL)
        VisualizeBtn['menu']=VisualizeBtn.menu

        return VisualizeBtn

    def makeAnalyzeMenu(self,mBar,network,namestring,root,weighted=FALSE):
        AnalyzeBtn=Menubutton(mBar,text="Analyze",underline=0,bg='steelblue')
        AnalyzeBtn.pack(side=LEFT,padx="2m")
        AnalyzeBtn.menu=Menu(AnalyzeBtn)

        AnalyzeBtn.menu.degrees=Menu(AnalyzeBtn.menu)
        AnalyzeBtn.menu.degrees.add_command(label='No binning',command=lambda s=self,net=network,n=namestring,r=root: s.calculate_degreedistribution(net,n,r,'lin'))
        AnalyzeBtn.menu.degrees.add_command(label='Logarithmic bins',command=lambda s=self,net=network,n=namestring,r=root: s.calculate_degreedistribution(net,n,r,'log'))
        AnalyzeBtn.menu.degrees.add_command(label='Cumulative',command=lambda s=self,net=network,n=namestring,r=root: s.calculate_degreedistribution(net,n,r,'cum'))

        AnalyzeBtn.menu.weights=Menu(AnalyzeBtn.menu)
        AnalyzeBtn.menu.weights.add_command(label='Linear bins',command=lambda s=self,net=network,n=namestring,r=root: s.weights_bin(net,n,r,'lin'))
        AnalyzeBtn.menu.weights.add_command(label='Cumulative',command=lambda s=self,net=network,n=namestring,r=root: s.weights_cumulative(net,n,r))
        AnalyzeBtn.menu.weights.add_command(label='Logarithmic bins',command=lambda s=self,net=network,n=namestring,r=root: s.weights_bin(net,n,r,'log'))
        #Disable logarithmic bins if data has zero distances
        if self.minw==0:
            AnalyzeBtn.menu.weights.entryconfig(3, state=DISABLED)

        AnalyzeBtn.menu.distributions=Menu(AnalyzeBtn.menu)

        AnalyzeBtn.menu.distributions.add_command(label='Shortest path distribution',command=lambda s=self,net=network,n=namestring,r=root: s.calculate_shortestpathdistribution())

        AnalyzeBtn.menu.distributions.add_cascade(label='Degree distribution',menu=AnalyzeBtn.menu.degrees)
        if not(weighted):
            AnalyzeBtn.menu.distributions.add_cascade(label='Distance distribution',state=DISABLED)
        else:
            AnalyzeBtn.menu.distributions.add_cascade(label='Distance distribution',menu=AnalyzeBtn.menu.weights)

        AnalyzeBtn.menu.clustering=Menu(AnalyzeBtn.menu)

        AnalyzeBtn.menu.clustering.add_command(label='No binning',command=lambda s=self,net=network,n=namestring,r=root: s.calculate_clustering(net,n,r,'lin'))
        AnalyzeBtn.menu.clustering.add_command(label='Logarithmic bins',command=lambda s=self,net=network,n=namestring,r=root: s.calculate_clustering(net,n,r,'log'))

        AnalyzeBtn.menu.assortativity=Menu(AnalyzeBtn.menu)
        AnalyzeBtn.menu.assortativity.options=Menu(AnalyzeBtn.menu)
        AnalyzeBtn.menu.assortativity.options.add_command(label='No binning',command=lambda s=self,net=network,n=namestring,r=root: s.knn(net,n,r,'lin'))
        AnalyzeBtn.menu.assortativity.options.add_command(label='Logarithmic bins',command=lambda s=self,net=network,n=namestring,r=root: s.knn(net,n,r,'log'))

        if self.weighted:
            AnalyzeBtn.menu.assortativity.w_options=Menu(AnalyzeBtn.menu)
            AnalyzeBtn.menu.assortativity.w_options.add_command(label='No binning',command=lambda s=self,net=network,n=namestring,r=root: s.knn(net,n,r,'lin',True))
            AnalyzeBtn.menu.assortativity.w_options.add_command(label='Logarithmic bins',command=lambda s=self,net=network,n=namestring,r=root: s.knn(net,n,r,'log',True))
        
        AnalyzeBtn.menu.assortativity.add_cascade(label='Average neighbor degree',menu=AnalyzeBtn.menu.assortativity.options)

        AnalyzeBtn.menu.percolation=Menu(AnalyzeBtn.menu)

        AnalyzeBtn.menu.percolation.weight_options=Menu(AnalyzeBtn.menu)
        AnalyzeBtn.menu.percolation.weight_options.add_command(label='Low distances added first',command=lambda s=self,n=namestring,r=root: s.percolation(r,namestring=n,method='weight',reverse=False))
        AnalyzeBtn.menu.percolation.weight_options.add_command(label='High distances added first',command=lambda s=self,n=namestring,r=root: s.percolation(r,namestring=n,method='weight',reverse=True))

        AnalyzeBtn.menu.percolation.fract_options=Menu(AnalyzeBtn.menu)
        AnalyzeBtn.menu.percolation.fract_options.add_command(label='Low distances added first',command=lambda s=self,n=namestring,r=root: s.percolation(r,namestring=n,method='fraction',reverse=False))
        AnalyzeBtn.menu.percolation.fract_options.add_command(label='High distances added first',command=lambda s=self,n=namestring,r=root: s.percolation(r,namestring=n,method='fraction',reverse=True))

        AnalyzeBtn.menu.percolation.add_cascade(label='By edge distance',menu=AnalyzeBtn.menu.percolation.weight_options)
        AnalyzeBtn.menu.percolation.add_cascade(label='By % of edges',menu=AnalyzeBtn.menu.percolation.fract_options)


        AnalyzeBtn.menu.add_cascade(label='Distributions',menu=AnalyzeBtn.menu.distributions)
        AnalyzeBtn.menu.add_cascade(label='Clustering',menu=AnalyzeBtn.menu.clustering)
        AnalyzeBtn.menu.add_cascade(label='Assortativity',menu=AnalyzeBtn.menu.assortativity)
        AnalyzeBtn.menu.add_cascade(label='Percolation',menu=AnalyzeBtn.menu.percolation)
        AnalyzeBtn['menu']=AnalyzeBtn.menu

        
        return AnalyzeBtn

    # DEFINE ANALYSIS COMMANDS

    def percolation(self,parent,namestring='',method='fraction',reverse=False):
        edges=list(self.net.edges)
        random.shuffle(edges)
        Nedges=len(edges)
        Nnodes=len(self.net)
        edges.sort(lambda x,y:cmp(x[2],y[2]),reverse=reverse)
        ee=percolator.EvaluationList(edges)
        if method=='weight':
            ee.setStrengthEvaluations()
        else:
            ee.setLinearEvaluations(0,len(edges),min(100,len(edges)))
        data=[]
        thresholds=[]
        threshfract=[]
        gcs=[]
        clusterings=[]
        susc=[]
        suscmax=0.0
        suscmax_thresh=0.0
        suscmax_fract=0.0
        counter=0
        gcctitle=''
        for cs in percolator.Percolator(ee):

            thresholds.append(cs.threshold)
            threshfract.append(cs.addedEdges/float(Nedges))
            gcs.append(cs.getGiantSize())
            s=cs.getSusceptibility(Nnodes)
            if (s>suscmax):
                suscmax=s
                suscmax_thresh=cs.threshold
                suscmax_fract=cs.addedEdges/float(Nedges)
            susc.append(s)

        if method=='weight':
            data.append([thresholds,gcs])
            data.append([thresholds,susc])
            tstring='w'
            susctitle='S peaks at f=%2.2f' % float(suscmax_thresh)
            if not(reverse):
                gcctitle='Links with d<threshold added'
            else:
                gcctitle='Links with d>threshold added'
                
        else:
            data.append([threshfract,gcs])
            data.append([threshfract,susc])
            tstring='%'
            susctitle='S peaks at f=%2.2f, w=%2.2f' % (float(suscmax_fract),float(suscmax_thresh))
            if not(reverse):
                gcctitle='% weakest links added'
            else:
                gcctitle='% strongest links added'
      
        t=DataWindow(parent,data,namestring+': percolation','plot',gcctitle+':'+susctitle,'threshold '+tstring,'GCC size:susceptibility')      

    def knn(self,network,namestring,parent,linorlog='log',weighted=False):
        '''Calculates and plots average nn degree'''

        if weighted:
            titlestring='Weighted avg nearest neighbour degree'
        else:
            titlestring='Average nearest neighbour degree'

        if linorlog=='log':
            plotstyle='loglog'
            choices=dialogues.AskNumberOfBins(parent,'Distance distribution: options')
            if choices.result==None:
                return
            if choices.result<11:
                tkMessageBox.showerror("Error","Number of bins must be larger than 10.")
                return
            Nbins=choices.result-10
            temp=netanalysis.knn_spectrum(network,'log',weighted,Nbins)

        else:
            plotstyle='plot'
            temp=netanalysis.knn_spectrum(network,'lin',weighted)

        t=DataWindow(parent,[[temp[0],temp[1]]],namestring,plotstyle,titlestring,'k','<knn(k)>')   

    def weights_cumulative(self,network,namestring,parent,maxPoints=1000):
        nDists=len(network.weights)
        dists=pylab.zeros(nDists)
        for i,w in enumerate(network.weights):
            dists[i]=w
        dists.sort()
        x=pylab.array(pylab.linspace(0,nDists,min(nDists,maxPoints)),dtype='uint')
        x[-1]=x[-1]-1
        dists=dists[x]
        x=pylab.array(x,dtype='float')/nDists
        DataWindow(self,[[dists,x]],namestring,'plot','Cumulative distance distribution','d','P(<d)')


    def weights_bin(self,network,namestring,parent,linorlog='log'):
        '''Calculates and plots weight distribution'''

        choices=dialogues.AskNumberOfBins(parent,'Distance distribution: options')
        if choices.result==None:
            return

        if linorlog=='lin':

            plotstyle='plot'
            temp=netanalysis.weight_distribution(network,'linbin',choices.result)

        else:

            plotstyle='loglog'
            temp=netanalysis.weight_distribution(network,'logbin',choices.result)

        t=DataWindow(parent,[[temp[0],temp[1]]],namestring,plotstyle,'Distance distribution','d','P(d)')
        
    def strength_bin(self,network,namestring,parent,linorlog='log'):
        '''Calculates and plots strength distribution'''
                
        choices=dialogues.AskNumberOfBins(parent,'Strength distribution: options')
        if choices.result==None:
            return

        if linorlog=='lin':
            plotstyle='plot'
            temp=netanalysis.strength_distribution(network,'linbin',choices.result)

        else:
            plotstyle='loglog'
            temp=netanalysis.strength_distribution(network,'logbin',choices.result)

        t=DataWindow(parent,[[temp[0],temp[1]]],namestring,plotstyle,'Strength distribution','s','P(s)')
        

    def calculate_clustering(self,network,namestring,parent,logorlin='log'):
        '''Calculates and plots clustering coefficient as function of degree'''

        if logorlin=='log':
             plottype='loglog'
             choices=dialogues.AskNumberOfBins(parent,'Clustering spectrum: options')
             if choices.result==None:
                 return
             if choices.result<11:
                 tkMessageBox.showerror("Error","Number of bins must be larger than 10.")
                 return
             Nbins=choices.result-10
             temp=netanalysis.clustering_spectrum(network,'logbin',Nbins)

        else:
            plottype='plot'
            temp=netanalysis.clustering_spectrum(network)

        t=DataWindow(parent,[[temp[0],temp[1]]],namestring,plottype,'Clustering spectrum','k','<c(k)>')


    def calculate_degreedistribution(self,network,namestring,parent,style='log'):

        titlestring='Degree distribution'
        ystring='P(k)'

        if style=='log':
            plottype='loglog'
            choices=dialogues.AskNumberOfBins(parent,'Degree distribution: options')            
            Nbins=choices.result
            if Nbins==None:
                return
            if choices.result<11:
                tkMessageBox.showerror("Error","Number of bins must be larger than 10.")
                return
            Nbins=Nbins-10
            temp=netanalysis.degree_distribution(network,'logbin',Nbins)

        elif style=='cum':
            plottype='loglog'
            temp=netanalysis.degree_distribution(network,'cumulative')
            titlestring='Cumulative degree distribution'
            ystring='P>(k)'

        else:
            plottype='plot'
            temp=netanalysis.degree_distribution(network,'linbin')

        t=DataWindow(parent,[[temp[0],temp[1]]],namestring,plottype,titlestring,'k',ystring)

    def calculate_shortestpathdistribution(self):
        dd=netext.getPathLengthDistribution(self.net,-1)
        diam=max(dd.keys())
        distance=range(1,diam+1)
        p=map(lambda x:dd.get(x,0),distance)
        avgdist=sum(map(lambda x,y:x*y,distance,p))
        leftCorner=map(lambda x:x-0.5,distance)
        t=DataWindow(self,[[leftCorner,p]],self.netname,'bar',"Shortest paths, diameter="+str(diam)+", avg l=%3.2f"%(avgdist),'l',"P(l)")


    # DEFINE NETWORK DERIVATION COMMANDS

    def get_coords(self,network):

        h=visuals.Himmeli(network)
        self.coords=h.getCoordinates()

    def threshold_manual(self,network,namestring,parent,mode):

         th=dialogues.AskThreshold(parent,'Choose threshold:')
         threshold=th.result

         if mode:
             titlestring='%s_edges_above_%3.2f' % (namestring,threshold)
         else:
             titlestring='%s_edges_below_%3.2f' % (namestring,threshold)

         newnet=transforms.threshold_by_value(network,threshold,mode)

         if len(self.coords)==0:

             # first time thresholding; generate MST visualization coords
             if self.nettype==1:
                 maximum=TRUE
             else:
                 maximum=FALSE
             mstnet=transforms.mst_kruskal(network,randomize=TRUE,maximum=maximum)
             h=visuals.Himmeli(mstnet)
             c=h.getCoordinates()
             self.coords=c             

         t=NetWindow(parent,newnet,titlestring,coords=self.coords,nettype=self.nettype)

    def generate_mst(self,network,namestring,parent,max_spanningtree):

        if max_spanningtree:
            titlestring='Max_spanning_tree_%s' % namestring
        else:
            titlestring='Min_spanning_tree_%s' % namestring

        newnet=transforms.mst_kruskal(network,randomize=TRUE,maximum=max_spanningtree)

        t=NetWindow(parent,newnet,titlestring,nettype=self.nettype)

    def generate_LCC(self,network,namestring,giantmembers,parent):

        titlestring=namestring+'_LCC'

        newnet=pynet.SymmNet()

        edgelist=list(network.edges)

        for edge in edgelist:

            if (edge[0] in giantmembers) and (edge[1] in giantmembers):

                newnet[edge[0]][edge[1]]=edge[2]

        t=NetWindow(parent,newnet,titlestring)

    def dist_to_weight(self,network,namestring,parent):

        netfilename='%s_DtoW' % namestring

        m=transforms.dist_to_weights(network)
        t=NetWindow(parent,m,netfilename,nettype=1)

    # DEFINE FILE COMMANDS

    def draw_MST_network(self,net,coords,namestring,root,parentnet):
        titlestring='Figure_%s' % namestring
        x=ImageWindow(root,net,titlestring,coords=coords,time=40,useMST=True,edgeColorMap='primary',parentnet=parentnet)


    def draw_network(self,net,namestring,root):
        titlestring='Figure_%s' % namestring

        N=len(net._nodes)

        if N>1000:
            time=80
        elif N>700:
            time=60
        elif N>400:
            time=40
        else:
            time=20
            
        x=ImageWindow(root,net,titlestring,sizeX=c[0],sizeY=c[0],time=time)
        
    def open_network(self):
        pass


    def export_nodes(self,network,namestring,parentnet=[]):

        t=self.netname.split('.')[0].replace(" ","_")+".txt"
        export_filename=tkFileDialog.asksaveasfilename(filetypes=[("Txt files",".txt"),("ASC files",".asc")],initialfile=t,defaultextension='.txt',title='Export as')
        if len(export_filename)==0:
            return

        f=open(export_filename,'w')
        hasproperties=hasattr(network,"nodeProperty")

        degs=netext.deg(network)
        clustering=netanalysis.clustering(network)

        f.write("Node\tDegree\tClustering\t")

        if hasproperties:
            p=list(network.nodeProperty)
            for prop in p:
                f.write(prop)
                f.write("\t")            
            f.write("\n")

        for node in network:
            f.write(str(node)+"\t"+str(degs[node])+"\t"+str(clustering[node])+"\t")
            if hasproperties:
                p=list(network.nodeProperty)
                for prop in p:
                    f.write(str(network.nodeProperty[prop][node]))
                    f.write("\t")
            f.write("\n")

        for node in parentnet:
            if not(node in network):
                f.write(str(node)+"\t0\t0\t")
                if hasattr(parentnet,"nodeProperty"):
                    p=list(parentnet.nodeProperty)
                    for prop in p:
                        f.write(str(parentnet.nodeProperty[prop][node]))
                        f.write("\t")
                f.write("\n")
        f.close()

        
    def export_components(self,network,namestring,parentnet=[]):
        """Writes contents of all network components into a text file"""
        t=self.netname.split('.')[0].replace(" ","_")+".txt"
        export_filename=tkFileDialog.asksaveasfilename(filetypes=[("Txt files",".txt"),("ASC files",".asc")],initialfile=t,defaultextension='.txt',title='Export as')
        if len(export_filename)==0:
            return
        f=open(export_filename,'w')
        components_temp=percolator.getComponents(network)        
        components=[]
        for comp in components_temp:
            components.append(list(comp))

        components.sort(key=len)
        components.reverse()

        hasproperties=hasattr(network,"nodeProperty")

        degs=netext.deg(network)
        clustering=netanalysis.clustering(network)

        compindex=0

        for index,comp in enumerate(components):
            # first write headers for this component
            f.write("Component %d: " % index)
            f.write("Size = %d nodes" % len(comp))
            f.write("\n")
            f.write("Node\tDegree\tClustering\t")
            if hasproperties:
                p=list(network.nodeProperty)
                for prop in p:
                    f.write(prop)
                    f.write("\t")            
            f.write("\n")

            for node in comp:
                f.write(str(node)+"\t"+str(degs[node])+"\t"+str(clustering[node])+"\t")

                if hasproperties:
                    p=list(network.nodeProperty)
                    for prop in p:
                        f.write(str(network.nodeProperty[prop][node]))
                        f.write("\t")
                f.write("\n")

            f.write("---\n")

            compindex=compindex+1

        for node in parentnet:
            if not(node in network):
                f.write("Component %d:" % compindex)
                f.write("Size = 1 nodes")
                f.write("\n")

                f.write("Node\tDegree\tClustering\t")
                if hasattr(parentnet,"nodeProperty"):
                    p=list(parentnet.nodeProperty)
                    for prop in p:
                        f.write(prop)
                        f.write("\t")            
                f.write("\n")
               
                f.write(str(node)+"\t0\t0\t")

                if hasattr(parentnet,"nodeProperty"):
                    p=list(parentnet.nodeProperty)
                    for prop in p:
                        f.write(str(parentnet.nodeProperty[prop][node]))
                        f.write("\t")

                f.write("\n")

                f.write("---\n")

                compindex=compindex+1                
        f.close()



    def save_network(self,network):
        namestring=self.netname.split('.')[0].replace(" ","_")

        filetypes=[("Edge list","*.edg"),("GML file","*.gml"),('Pajek file','*.net'),("Matrix format","*.mat")]
        types=["edg","gml","net","mat"]
        defaultType=types[0]

        newname=tkFileDialog.asksaveasfilename(initialfile=namestring+"."+defaultType,title='Save as',filetypes=filetypes,defaultextension="."+defaultType)
        if len(newname)>0: #if user did not press cancel button
            ending=newname.split(".")[-1]
            if ending in types:
                netio.writeNet(network,newname,headers=False)
            else:
                tkMessageBox.showerror(
                    "Error saving the figure",
                    "Unknown file format. Use edg, gml, net or mat."
                    )

    def exit_network(self):
        self.destroy()


        
