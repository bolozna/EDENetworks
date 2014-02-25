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
import pylab
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2TkAgg
from PIL import Image
from PIL import ImageTk
import os
import shutil,math
import random
from EdenNetWindow import EdenNetWindow
from DataWindow import DataWindow



class RawDataWindow(Toplevel):
    def __init__(self,parent,distancematrix,namestring=None,sizeX=800,sizeY=600,matrixtype='distance',distancemeasure='unknown',msdata=None,poplist=None,nodeTypes='individual'):
        """
        nodeTypes : individual or population
        """
        
        Toplevel.__init__(self,parent,bg='gray90')  # creates the window itself

        self.distancematrix=distancematrix
        self.parent=parent
        self.weighted=True
        self.matrixtype=0  # meaning this is a distance matrix instead of a weight matrix (type=1)
        self.coords=[]
        self.namestring=""
        self.msdata=msdata
        self.poplist=poplist
        self.nodeTypes=nodeTypes
        self.distancemeasure=distancemeasure
        
        if namestring:
            self.title("Raw data: "+namestring)
            self.namestring=namestring

        mBar=Frame(self,relief='raised',borderwidth=2,bg='steelblue') # this will be the menu bar

        lowerhalf=Frame(self,bg='gray90')

        # --- calculate main statistics etc
        N=len(self.distancematrix)
        [minw,maxw,avgw]=netanalysis.weightStats(self.distancematrix)
        self.minw=minw #saved for disabling log binning from menu bars 

        # --- menu bar ------------------------------
        mBar=Frame(self,relief='raised',borderwidth=2,bg='steelblue')
        FileBtn=self.makeFileMenu(mBar,self.distancematrix,namestring,self.parent)
        AnalyzeBtn=self.makeAnalyzeMenu(mBar,self.distancematrix,namestring,self.parent,self.weighted)
        DeriveBtn=self.makeDeriveMenu(mBar,self.distancematrix,namestring,self.parent,self.weighted)
        if self.nodeTypes=="population":
            PopBtn=self.makePopulationMenu(mBar)
        
        if self.nodeTypes=="population":
            mBar.tk_menuBar=(FileBtn,DeriveBtn,AnalyzeBtn,PopBtn)
        else:
            mBar.tk_menuBar=(FileBtn,DeriveBtn,AnalyzeBtn)


        # --- lower left: plot panel ---------------

        # generate plot frame
        plotframe=Frame(lowerhalf,relief='raised',borderwidth=2,bg='gray80')

        # first generate plot
        temp=netanalysis.weight_distribution(self.distancematrix,'linbin',14)
        datavector=[[temp[0],temp[1]]]

        width=0.75*(temp[0][-1]-temp[0][0])/float(len(temp[0]))

        myplot=visuals.ReturnPlotObject(datavector,'bar',"sample distances","d","P(d)",figsize=(3,2),fontsize=9,addstr=',width='+str(width)+',color=\"#9e0b0f\"')
        myplot.canvas=FigureCanvasTkAgg(myplot.thisFigure,master=plotframe)
        myplot.canvas.show()
        myplot.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=YES, padx=10,pady=10) 

        plotframe.grid(row=0,column=0,sticky=NW)

        # ----- lower right: info panel -------
        infoframe=Frame(lowerhalf,relief='raised',borderwidth=2,bg='gray80')
        lwidth=20 # sets width for the following labels
        rowcounter=0
        Label(infoframe,
              text="Data:",
              relief='flat',
              bg='DarkOliveGreen2',
              justify=RIGHT,width=lwidth).grid(row=rowcounter)
        Label(infoframe,
              text="%s" % namestring,
              relief='flat',
              bg='DarkOliveGreen3',
              width=lwidth).grid(row=rowcounter,column=1)
        rowcounter=rowcounter+1
        Label(infoframe,text="Size:",relief='flat',bg='gray70',width=lwidth).grid(row=rowcounter)
	sizeunit="samples"
	if self.nodeTypes=="population":
		sizeunit="populations"

        Label(infoframe,text="N=%d " % N + sizeunit,relief='flat',bg='gray60',width=lwidth).grid(row=rowcounter,column=1)
	rowcounter=rowcounter+1
       
	if hasattr(self.distancematrix,'N_origsamples'):
		Label(infoframe,text="Based on:",relief='flat',bg='gray70',width=lwidth).grid(row=rowcounter)
		Label(infoframe,text="Ns=%d samples" % self.distancematrix.N_origsamples,relief='flat',bg='gray60',width=lwidth).grid(row=rowcounter,column=1)
		rowcounter=rowcounter+1


	if not(self.nodeTypes=="population"):
        # label informing about how clones have been handled
	        Label(infoframe,text="Clones:",relief='flat',bg='gray70',width=lwidth).grid(row=rowcounter)

	        if hasattr(self.distancematrix,'clones'):      
        	    if self.distancematrix.clones=='collapsed':
                	if hasattr(self.distancematrix,'Nclones'):
                    	   clonetext='removed %d samples' % self.distancematrix.Nclones
                	else:
                  	   clonetext='removed'
           	    elif self.distancematrix.clones=='included':
                	clonetext='kept'
	            else:
               	        clonetext='unknown'
        	else:
              	    clonetext='unknown'
        
        
        	Label(infoframe,text="%s" % clonetext,relief='flat',bg='gray60',width=lwidth).grid(row=rowcounter,column=1)
        	rowcounter=rowcounter+1


        # label about distance measure
        if hasattr(self.distancematrix,'distancemeasure'):
            distancetext=self.distancematrix.distancemeasure
        else:
            distancetext='unknown'

        Label(infoframe,text="Distance measure:",relief='flat',bg='gray65',width=lwidth).grid(row=rowcounter)
        Label(infoframe,text="%s" % distancetext,relief='flat',bg='gray55',width=lwidth).grid(row=rowcounter,column=1)
        rowcounter=rowcounter+1

        Label(infoframe,text="Average distance:",relief='flat',bg='gray65',width=lwidth).grid(row=rowcounter)
        Label(infoframe,text="<d>=%3.2f" % avgw,relief='flat',bg='gray55',width=lwidth).grid(row=rowcounter,column=1)
        rowcounter=rowcounter+1
        
        Label(infoframe,text="Min distance:",relief='flat',bg='gray65',width=lwidth).grid(row=rowcounter)
        Label(infoframe,text="dmin=%3.2f" % minw,relief='flat',bg='gray55',width=lwidth).grid(row=rowcounter,column=1)
        rowcounter=rowcounter+1

        Label(infoframe,text="Max distance:",relief='flat',bg='gray65',width=lwidth).grid(row=rowcounter)
        Label(infoframe,text="dmax=%3.2f" % maxw,relief='flat',bg='gray55',width=lwidth).grid(row=rowcounter,column=1)
        rowcounter=rowcounter+1

        if hasattr(self.distancematrix,'nodeProperty'):
            Label(infoframe,text="Imported attributes", relief='flat',bg='DarkOliveGreen2',width=lwidth).grid(row=rowcounter)
            Label(infoframe,text="type",relief='flat',bg='DarkOliveGreen3',width=lwidth).grid(row=rowcounter,column=1)
            rowcounter=rowcounter+1
            propertydict=netext.getPropertyTypes(self.distancematrix)
            for prop in propertydict.keys():
                Label(infoframe,text=prop,relief='flat',bg='gray70',width=lwidth).grid(row=rowcounter)
                Label(infoframe,text=propertydict[prop],relief='flat',bg='gray65',width=lwidth).grid(row=rowcounter,column=1)
                rowcounter=rowcounter+1


        # -- SHOW IT ALL ---
        mBar.pack(side=TOP,fill=X,expand=YES,anchor=NW)
        infoframe.grid(row=0,column=1,sticky=NW)
        lowerhalf.pack(side=TOP,fill=BOTH,expand=YES)

    # ------- generators for menus ----------

    def makeFileMenu(self,mBar,network,namestring,root):
        FileBtn=Menubutton(mBar,text="File",underline=0,bg='steelblue')
        FileBtn.pack(side=LEFT,padx="2m")
        FileBtn.menu=Menu(FileBtn)

        FileBtn.menu.save=Menu(FileBtn.menu)
        FileBtn.menu.save.add_command(label="Square matrix",underline=0,command=lambda s=self,type="square":s.save_network(type))
        FileBtn.menu.save.add_command(label="Upper triangular",underline=0,command=lambda s=self,type="upperdiag":s.save_network(type))
        FileBtn.menu.save.add_command(label="Strictly upper triangluar",underline=0,command=lambda s=self,type="supperdiag":s.save_network(type))
        FileBtn.menu.save.add_command(label="Lower triangular",underline=0,command=lambda s=self,type="lowerdiag":s.save_network(type))
        FileBtn.menu.save.add_command(label="Strictly lower triangluar",underline=0,command=lambda s=self,type="slowerdiag":s.save_network(type))
        FileBtn.menu.add_cascade(label='Save distance matrix',menu=FileBtn.menu.save)

        FileBtn.menu.add('separator')
        FileBtn.menu.add_command(label="Close",underline=0,command=self.exit_network)
        FileBtn['menu']=FileBtn.menu
        return FileBtn

    def makeDeriveMenu(self,mBar,network,namestring,root,weighted=FALSE):
        DeriveBtn=Menubutton(mBar,text="Derive",underline=0,bg='steelblue')
        DeriveBtn.pack(side=LEFT,padx="2m")
        DeriveBtn.menu=Menu(DeriveBtn)

        DeriveBtn.menu.add_command(label='Minimum spanning tree',command=lambda s=self,net=self.distancematrix,n=namestring,r=root: s.generate_mst(net,n,r,False))
        DeriveBtn.menu.add_command(label='Manual thresholding',command=lambda s=self,n=namestring,r=root: s.percolation(r,namestring=n,method='weight',reverse=False))
        DeriveBtn.menu.add_command(label='Automatic thresholding',command=self.autothreshold)

        DeriveBtn["menu"]=DeriveBtn.menu
         
    def makeAnalyzeMenu(self,mBar,network,namestring,root,weighted=FALSE):
        AnalyzeBtn=Menubutton(mBar,text="Analyze",underline=0,bg='steelblue')
        AnalyzeBtn.pack(side=LEFT,padx="2m")
        AnalyzeBtn.menu=Menu(AnalyzeBtn)

        AnalyzeBtn.menu.weights=Menu(AnalyzeBtn.menu)
        AnalyzeBtn.menu.weights.add_command(label='Linear bins',command=lambda s=self,net=self.distancematrix,n=namestring,r=root: s.weights_bin(net,n,r,'lin'))
        AnalyzeBtn.menu.weights.add_command(label='Cumulative',command=lambda s=self,net=self.distancematrix,n=namestring,r=root: s.weights_cumulative(net,n,r))
        AnalyzeBtn.menu.weights.add_command(label='Logarithmic bins',command=lambda s=self,net=self.distancematrix,n=namestring,r=root: s.weights_bin(net,n,r,'log'))

        #Disable logarithmic bins if data has zero distances
        if self.minw==0:
            AnalyzeBtn.menu.weights.entryconfig(3, state=DISABLED)

        AnalyzeBtn.menu.add_cascade(label='Distance distribution',menu=AnalyzeBtn.menu.weights)
        AnalyzeBtn['menu']=AnalyzeBtn.menu
        
        return AnalyzeBtn

    def makePopulationMenu(self,mBar):
        PopBtn=Menubutton(mBar,text="Randomizations",underline=0,bg='steelblue')
        PopBtn.pack(side=LEFT,padx="2m")
        PopBtn.menu=Menu(PopBtn)
        PopBtn.menu.stats=Menu(PopBtn.menu)
        PopBtn.menu.single=Menu(PopBtn.menu)

        PopBtn.menu.single.add_command(label='Bootstrapping',command=lambda: self.singleBootstrapping())
        PopBtn.menu.stats.add_command(label='Betweenness: Bootstrapping',command=lambda: self.bootstrapping())
        PopBtn.menu.single.add_command(label='Shuffle samples',command=lambda: self.singleShuffledNodes())
        PopBtn.menu.stats.add_command(label='Clustering: Shuffle samples',command=lambda: self.statsShuffled("nodes"))
        PopBtn.menu.single.add_command(label='Shuffle alleles',command=lambda: self.singleShuffledAlleles())
        PopBtn.menu.stats.add_command(label='Clustering: Shuffle alleles',command=lambda:  self.statsShuffled("alleles"))
        
        PopBtn.menu.add_cascade(label="Single realizations",menu=PopBtn.menu.single)        
        PopBtn.menu.add_cascade(label="Statistics",menu=PopBtn.menu.stats)        
        PopBtn['menu']=PopBtn.menu
        
        return PopBtn

    # DEFINE ANALYSIS COMMANDS
    def bootstrapping(self):
        bsdialog=dialogues.BootstrapPopulationsDialog(self)
        bsP,rounds=bsdialog.result
        bsP,rounds=float(bsP),int(rounds)
        
        bsWaiter=dialogues.WaitWindow(self,title="Processsing...",titlemsg="Calculating statistics...",hasProgressbar=True)
        bsFunction=lambda:bootstrap(self.msdata,self.poplist,self.distancematrix,bsP,rounds,bsWaiter.progressbar.set)
        bsWaiter.go(bsFunction)
        nodeBc,nodeBcRank,nodeDegree,nodeDegreeRank=bsWaiter.result

        BootstrapResultsWindow(self.parent,nodeBc,nodeBcRank,"betweenness centrality","BC",namestring=self.namestring,bsP=bsP)

    def singleBootstrapping(self):
        bsdialog=dialogues.SliderDialog2(self,0.0,1.0,0.01,0.5,bodyText="Percentage of nodes in each location?")
        bsP=bsdialog.result
        bsNet,bsMsdata,bsPoplist=bootstrap_single(self.msdata,self.poplist,self.distancematrix,bsP)
	bsNet.matrixtype=0
        RawDataWindow(self.parent,bsNet,namestring=self.namestring+"_bootstrap_P="+str(bsP),sizeX=800,sizeY=600,matrixtype=self.matrixtype,distancemeasure=self.distancemeasure,msdata=bsMsdata,poplist=bsPoplist,nodeTypes=self.nodeTypes)

    def autothreshold(self):
        newnet,threshold=autoThreshold(self.distancematrix,outputThreshold=True)
        netext.copyNodeProperties(self.distancematrix,newnet)
        titlestring=self.namestring+" thresholded at %2.2f" % threshold
        self.makeMstAndCoords()
        EdenNetWindow(self.parent,newnet,namestring=titlestring,nettype=self.matrixtype,parentnet=self.distancematrix,coords=self.coords,mstnet=self.mstnet)

    def singleShuffledNodes(self):
        newmsdata=self.msdata.copy() #make a copy of the data
        newmsdata.shuffleNodes()
        goldstein_list,unique_poplist=eden.getGoldsteinLists(self.poplist)
        newdm=newmsdata.getGroupwiseDistanceMatrix(goldstein_list,self.distancematrix.distancemeasure,groupNames=unique_poplist)
        newdm.distancemeasure=self.distancematrix.distancemeasure
	newdm.matrixtype=0
        RawDataWindow(self.parent,newdm,namestring=self.namestring+"_shuffledSamples",sizeX=800,sizeY=600,matrixtype=self.matrixtype,distancemeasure=self.distancemeasure,msdata=newmsdata,poplist=self.poplist,nodeTypes=self.nodeTypes)

    def singleShuffledAlleles(self):
        newmsdata=self.msdata.copy() #make a copy of the data
        newmsdata.randomize()
        goldstein_list,unique_poplist=eden.getGoldsteinLists(self.poplist)
        newdm=newmsdata.getGroupwiseDistanceMatrix(goldstein_list,self.distancematrix.distancemeasure,groupNames=unique_poplist)
        newdm.distancemeasure=self.distancematrix.distancemeasure
	newdm.matrixtype=0
        RawDataWindow(self.parent,newdm,namestring=self.namestring+"_shuffledAlleles",sizeX=800,sizeY=600,matrixtype=self.matrixtype,distancemeasure=self.distancemeasure,msdata=newmsdata,poplist=self.poplist,nodeTypes=self.nodeTypes)

    def statsShuffled(self,shuffling):
        #first autothreshold the original network
        thNet,autoTh=autoThreshold(self.distancematrix,outputThreshold=True)

        repDialog=dialogues.AskNumberDialog(self,title="Provide a treshold leveel",bodyText="Threshold distance:",initNumber=autoTh)
        try:
            th=float(repDialog.result)
        except Exception:
            tkMessageBox.showerror(
                    "Error",
                    "Please provide a number."
                    )
            return
       
        thNet=transforms.threshold_by_value(self.distancematrix,th,accept="<=",keepIsolatedNodes=True)

        thEdges=len(thNet.edges)
        oClustering=netanalysis.globalClustering(thNet)

        #ask how many repetitions
        repDialog=dialogues.AskNumberDialog(self,title="Provide number of repetitions",bodyText="Number of repetitions:")
        try:
            reps=int(repDialog.result)
        except Exception:
            tkMessageBox.showerror(
                    "Error",
                    "Please provide a number."
                    )
            return
        
        #show progressbar
        waiter=dialogues.WaitWindow(self,title="Processing...",titlemsg="Calculating statistics...",hasProgressbar=True)

        #calculate stats
        clustering=[]
        for round in range(reps):
            newmsdata=self.msdata.copy() #make a copy of the data
            if shuffling=="nodes":
                newmsdata.shuffleNodes()
            elif shuffling=="alleles":
                newmsdata.randomize()
            else:
                raise Exception("No such shuffling method.")
            goldstein_list,unique_poplist=eden.getGoldsteinLists(self.poplist)
            newdm=newmsdata.getGroupwiseDistanceMatrix(goldstein_list,self.distancematrix.distancemeasure,groupNames=unique_poplist)
            newdm.distancemeasure=self.distancematrix.distancemeasure            
            
            newEdges=list(newdm.edges)
            newEdges.sort(key=lambda x:x[2])
            newThNet=pynet.SymmNet()
            for i in range(thEdges):
                edge=newEdges[i]
                newThNet[edge[0],edge[1]]=edge[2]

            clustering.append(netanalysis.globalClustering(newThNet))
            
            waiter.progressbar.set(float(round)/float(reps))

        #terminate the progressbar window
        waiter.ok()

        p=float(len(filter(lambda x:x>=oClustering,clustering)))/float(len(clustering))
        try: #for numpy <1.3 ? 
            temp=pylab.histogram(clustering,bins=min(reps,30),new=True)
        except TypeError:            
            temp=pylab.histogram(clustering,bins=min(reps,30))
        temp=list(temp)
        temp[1]=temp[1][:len(temp[1])-1]
        width=temp[1][1]-temp[1][0]
        t=DataWindow(self,[[temp[1],temp[0]]],self.namestring+"_shuffle_"+shuffling+"_"+str(reps),'bar','Clustering (original=%.2f,p=%g)' %(oClustering,p),'<c>','P(<c>)',addstring=",width="+str(width))


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
        if choices.result!=None:
            if linorlog=='lin':
                plotstyle='bar'
                temp=netanalysis.weight_distribution(network,'linbin',choices.result)
                width=0.75*(temp[0][-1]-temp[0][0])/float(len(temp[0]))
                addstring=',width='+str(width)+',color=\"#9e0b0f\"'

            else:

                plotstyle='loglog'
                temp=netanalysis.weight_distribution(network,'logbin',choices.result)
                addstring=''


            t=DataWindow(parent,[[temp[0],temp[1]]],namestring,plotstyle,'Distance distribution','d','P(d)',addstring=addstring)

    def percolation(self,parent,namestring='',method='fraction',reverse=False):
        edges=list(self.distancematrix.edges)
        random.shuffle(edges)

        Nedges=len(edges)
        Nnodes=len(self.distancematrix)
        edges.sort(lambda x,y:cmp(x[2],y[2]),reverse=reverse)
        ee=percolator.EvaluationList(edges)
        if method=='weight':
            ee.setStrengthEvaluations()
        else:
            ee.setLinearEvaluations(0,len(edges),100)
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
                gcctitle='Links with w<threshold added'
            else:
                gcctitle='Links with w>threshold added'
                
        else:
            data.append([threshfract,gcs])
            data.append([threshfract,susc])
            tstring='%'
            susctitle='S peaks at f=%2.2f, w=%2.2f' % (float(suscmax_fract),float(suscmax_thresh))
            if not(reverse):
                gcctitle='% weakest links added'
            else:
                gcctitle='% strongest links added'

        if len(self.coords)==0:
            mstnet=transforms.mst_kruskal(self.distancematrix,randomize=TRUE,maximum=FALSE)
            mstnet.matrixtype=0
            h=visuals.Himmeli(mstnet,treeMode=True)
            c=h.getCoordinates()
            self.coords=c
            self.mstnet=mstnet

        #Use autoThreshold results instead of max susc values
        newnet,suscmax_thresh=autoThreshold(self.distancematrix,outputThreshold=True)
        #To avoid rounding errors, not a very nice thing to do :(
        suscmax_thresh+=suscmax_thresh*0.00000001
        
        temp=dialogues.PercolationDialog(parent,title="Choose threshold distance:",titlemsg="Set threshold distance",pdata=data,suscmax_thresh=suscmax_thresh)

        if temp.result!=None:
            threshold=float(temp.result)
            newnet=transforms.threshold_by_value(self.distancematrix,threshold,accept="<=",keepIsolatedNodes=True)
            newnet.matrixtype=0
            titlestring=self.namestring+" thresholded at %2.2f" % threshold

            if len(newnet.edges)==0:
                tkMessageBox.showerror("Error while thresholding","Threshodling lead to an empty network.\nIncrease the threshold level.")
            else:
                t=EdenNetWindow(parent,newnet,titlestring,coords=self.coords,nettype=self.matrixtype,mstnet=self.mstnet,parentnet=self.distancematrix)


    def makeMstAndCoords(self):
        if len(self.coords)==0:
            mstnet=transforms.mst_kruskal(self.distancematrix,randomize=True,maximum=False)
            mstnet.matrixtype=0
            h=visuals.Himmeli(mstnet,treeMode=True)
            c=h.getCoordinates()
            self.coords=c
            self.mstnet=mstnet
        
    def generate_mst(self,network,namestring,parent,max_spanningtree):
        if max_spanningtree:
            titlestring='Max_spanning_tree_%s' % namestring
        else:
            titlestring='Min_spanning_tree_%s' % namestring

        newnet=transforms.mst_kruskal(network,randomize=True,maximum=max_spanningtree)
        newnet.matrixtype=0

        h=visuals.Himmeli(newnet,treeMode=True)
        c=h.getCoordinates()

        filename=h.netName+"_0001.eps"

        if os.path.isfile(filename):
            os.remove(filename)
        
        t=EdenNetWindow(parent,newnet,titlestring,nettype=self.matrixtype,parentnet=network,coords=c,mstnet=newnet)

    def threshold_manual(self,network,namestring,parent,mode):

         th=dialogues.AskThreshold(parent,'Choose threshold:')
         threshold=th.result

         if mode:
             titlestring='%s_edges_above_%3.2f' % (namestring,threshold)
         else:
             titlestring='%s_edges_below_%3.2f' % (namestring,threshold)

         newnet=transforms.threshold_by_value(network,threshold,mode)
	 newnet.matrixtype=0

         if len(self.coords)==0:

             # first time thresholding; generate MST visualization coords
             if self.matrixtype==1:
                 maximum=TRUE
             else:
                 maximum=FALSE
             mstnet=transforms.mst_kruskal(network,randomize=TRUE,maximum=maximum)
             h=visuals.Himmeli(mstnet)
             c=h.getCoordinates()
             self.coords=c
             self.mstnet=mstnet
           
         t=NetWindow(parent,newnet,titlestring,coords=self.coords,nettype=self.matrixtype,mstnet=self.mstnet)

    def dist_to_weight(self,network,namestring,parent):
        matrixfilename='%s_DtoW' % namestring
        m=transforms.dist_to_weights(network)
        t=MatrixWindow(parent,m,matrixfilename,matrixtype=1)
                          
    def open_network(self):
        pass

    def load_aux_data(self):
        aux_filename=tkFileDialog.askopenfilename(filetypes=[("Text files",".txt"),
                                                         ("Ascii files",".asc")],
                                              title="Choose auxiliary data file")
        if len(aux_filename)==0:
            return
        else:
            try:
                netio.loadNodeProperties(self.distancematrix,aux_filename)
            except Exception,e:
                tkMessageBox.showerror(
                    "Error importing auxiliary data:",
                    "File in wrong format:\n %s" % str(e)
                    )
        #call for method updating the property labels in the window


    def save_network(self,type="square"):
        network=self.distancematrix
        netfilename=tkFileDialog.asksaveasfilename(title="Select file for the distance matrix")
        if len(netfilename)==0:
            return
        netfile=open(netfilename,"w")
        nodes=netio.writeNet_mat(network,netfile,type=type)
        nodefilename=tkFileDialog.asksaveasfilename(title="Select file for node names")
        if len(nodefilename)>0:
            nodefile=open(nodefilename,"w")
            for node in nodes:
                nodefile.write(str(node)+"\n")
            tkMessageBox.showinfo(message='Save succesful.')
        else:
            tkMessageBox.showinfo(message='Save succesful.\nNode names not saved.')

    def exit_network(self):
        self.destroy()    

    def setAxes(self,a,c,xstring,ystring):
        setp(a,'xscale',xstring,'yscale',ystring)
        c.show()          

    def close_window(self):
        self.destroy()

class BootstrapResultsWindow(Toplevel):
    def __init__(self,parent,nodeM,nodeMRank,measureName,measureShorthand,namestring=None,sizeX=8000,sizeY=8000,bsP=None):
        Toplevel.__init__(self,parent,bg='gray90')  # creates the window itself
        self.parent=parent
        if namestring!=None:
            self.title("Bootstrapping results for "+measureName+", "+namestring+" with "+str(100*bsP)+"% of samples kept")
            self.namestring=namestring
        else:
            self.namestring=""
            
        mBar=Frame(self,relief='raised',borderwidth=2,bg='steelblue') # this will be the menu bar
        lowerhalf=Frame(self,bg='gray90')

        # --- menu bar ------------------------------
        mBar=Frame(self,relief='raised',borderwidth=2,bg='steelblue')
        FileBtn=self.makeFileMenu(mBar,self.parent)

        # --- left: bc histograms ---------------
        # generate plot frame
        plotframe=Frame(lowerhalf,relief='raised',borderwidth=2,bg='gray80',width=2000, height=2000)
        #generate the figure
        fig=self.getMFigure(nodeM,measureName,measureShorthand)
        self.figure_left=fig
        #Show the figure
        canvas=FigureCanvasTkAgg(fig,master=plotframe)
        canvas.show()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=YES, padx=10,pady=10)
        plotframe.grid(row=0,column=0,sticky=NW)

        # --- right: bc rank histograms ---------------
        # generate plot frame
        rightFrame=Frame(lowerhalf,relief='raised',borderwidth=2,bg='gray80',width=2000, height=2000)

        plotframe=Frame(rightFrame,relief='raised',borderwidth=2,bg='gray80',width=2000, height=2000)
        fig=self.getTopRankedFigure(nodeMRank,measureShorthand)
        canvas=FigureCanvasTkAgg(fig,master=plotframe)
        canvas.show()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=YES, padx=10,pady=10)
        plotframe.grid(row=0,column=0,sticky=NW)

        plotframe=Frame(rightFrame,relief='raised',borderwidth=2,bg='gray80',width=2000, height=2000)
        self.figure_upper_right=fig

        fig=self.getTopRankedFigure(nodeMRank,measureShorthand,nTop=1)
        canvas=FigureCanvasTkAgg(fig,master=plotframe)
        canvas.show()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=YES, padx=10,pady=10)
        plotframe.grid(row=1,column=0,sticky=NW)
        self.figure_lower_right=fig

        rightFrame.grid(row=0,column=1,sticky=NW)
        #rightFrame.pack(side=TOP,fill=BOTH,expand=YES)



        # -- SHOW IT ALL ---

        mBar.pack(side=TOP,fill=X,expand=YES,anchor=NW)
        lowerhalf.pack(side=TOP,fill=BOTH,expand=YES)

    def getMFigure(self,nodeBc,measureName,measureShorthand):
        #sort bc by mean
        names=nodeBc.keys()
        meanNodeBc={}
        for name in names:
            meanNodeBc[name]=pylab.mean(nodeBc[name])
        names.sort(key=lambda x:meanNodeBc[x],reverse=True)
        data=[]
        for name in names:
            data.append(nodeBc[name])

        #top 5 distributions
        nTop=5
        fig=pylab.Figure(figsize=(5,8), dpi=80)
        fig.subplots_adjust(bottom=0.08,right=0.95,top=0.95)
        for nodeIndex in range(nTop):
            axes=fig.add_subplot(nTop,1,nodeIndex+1)
            axes.hist(data[nodeIndex],100)
            axes.set_ylabel(names[nodeIndex],fontsize=8)
            for tick in axes.get_yticklabels():
                tick.set_fontsize(10)
                tick.set_fontname("Times")
            if nodeIndex==0:
                axes.set_title("Distribution of "+measureShorthand+"s for top "+str(nTop)+" locations")
        axes.set_xlabel(measureName)
            
        return fig

    def getTopRankedFigure(self,nodeBcRank,measureShorthand,nTop=5,plotAtMost=10):
        tops={}
        for node in nodeBcRank.iterkeys():
            for rank in nodeBcRank[node]:
                if rank<nTop:
                    tops[node]=tops.get(node,0)+1
        names=tops.keys()
        names.sort(key=lambda x:tops[x],reverse=True)
        data=[]
        for name in names:
            data.append(tops[name])

        fig=pylab.Figure(figsize=(5,3.8),dpi=80)
        fig.subplots_adjust(bottom=0.3)
        axes=fig.add_subplot(111)
        nums=range(min(len(data),plotAtMost))
        axes.bar(nums,data[:plotAtMost])
        axes.set_xticklabels(names,rotation=45,fontsize=8)
        axes.set_xticks(map(lambda x:x+0.5,nums))
        axes.set_title("# of times at top "+str(nTop))

        return fig

    def makeFileMenu(self,mBar,root):
        FileBtn=Menubutton(mBar,text="File",underline=0,bg='steelblue')
        FileBtn.pack(side=LEFT,padx="2m")
        FileBtn.menu=Menu(FileBtn)
        FileBtn.menu.add_command(label="Save figure on the left...",underline=0,command=self.save_left)
        FileBtn.menu.add_command(label="Save figure on the upper right...",underline=0,command=self.save_upper_right)
        FileBtn.menu.add_command(label="Save figure on the lower right...",underline=0,command=self.save_lower_right)
        FileBtn.menu.add_command(label="Close",underline=0,command=self.destroy)        
        FileBtn['menu']=FileBtn.menu
        return FileBtn

    def save_left(self):
        self.save_figure("left.png",self.figure_left)
    def save_upper_right(self):
        self.save_figure("upper_right.png",self.figure_upper_right)
    def save_lower_right(self):
        self.save_figure("lower_right.png",self.figure_lower_right)

    def save_figure(self,namestring,thefigure):
        filetypes=[("PNG file","*.png"),("PDF file","*.pdf"),('EPS file','*.eps'),("SVG file","*.svg")]
        types=["png","pdf","eps","svg"]
        newname=tkFileDialog.asksaveasfilename(initialfile=namestring,title='Save as',filetypes=filetypes)
        if len(newname)>0: #if user did not press cancel button
            ending=newname.split(".")[-1]
            if ending in types:
                thefigure.savefig(newname)
            else:
                tkMessageBox.showerror(
                    "Error saving the figure",
                    "Unknown file format. Use png, pdf, eps or svg."
                    )



def bootstrap(msdata,poplist,distancematrix,bsP,rounds,progressUpdater):
    #find the threhold for the original distance matrix
    othrNet=autoThreshold(distancematrix)
    originalThreshold=len(othrNet.edges)

    #transfer poplist to another format:
    pops={}
    for i,p in enumerate(poplist):
        if p not in pops:
            pops[p]=[]
        pops[p].append(i)

    #make node statistics containers:
    nodeBc={}
    nodeBcRank={}
    nodeDegree={}
    nodeDegreeRank={}
    for node in distancematrix:
        nodeBc[node]=[]
        nodeBcRank[node]=[]
        nodeDegree[node]=[]
        nodeDegreeRank[node]=[]

    #rounds
    for r in range(rounds):
        #first select nodes
        selected=[]
        for pop in pops:
            nPop=len(pops[pop])
            bsSize=int(bsP*nPop)
            if bsSize<1:
                bsSize=1
            selected.extend(random.sample(pops[pop],bsSize))
        selected.sort()

        #make new distance matrix
        newmsdata=msdata.getSubset(selected)
        newpoplist=[]
        for index in selected:
            newpoplist.append(poplist[index])
        newglists,newunique_poplist=eden.getGoldsteinLists(newpoplist)
        newdm=newmsdata.getGroupwiseDistanceMatrix(newglists,distancematrix.distancemeasure,groupNames=newunique_poplist)
        newdm.matrixtype=distancematrix.matrixtype
        
        #--BC statistics
        thNet=autoThreshold(newdm)
        #calculate statistics for this realization
        bc=netext.getBetweennessCentrality(thNet)
        bcRankList=list(distancematrix)
        random.shuffle(bcRankList)
        bcRankList.sort(key=lambda x:-bc[x])
        for rank,node in enumerate(bcRankList):
            nodeBcRank[node].append(rank)
        for node in distancematrix:
            nodeBc[node].append(bc[node])

        #--Degree statistics
        othNet=pynet.SymmNet()
        for node in newdm:
            othNet.addNode(node)
        edges=list(newdm.edges)
        edges.sort(lambda x,y:cmp(x[2],y[2]),reverse=False)
        for i in range(originalThreshold):
            othNet[edges[i][0],edges[i][1]]=edges[i][2]
        #calculate statistics
        degrees=netext.deg(othNet)
        for node in othNet:
            nodeDegree[node].append(degrees[node])
        degreeRankList=list(distancematrix)
        random.shuffle(degreeRankList)
        degreeRankList.sort(key=lambda x:-degrees[x])
        for rank,node in enumerate(degreeRankList):
            nodeDegreeRank[node].append(rank)
        
        progressUpdater(float(r)/float(rounds))

    return nodeBc,nodeBcRank,nodeDegree,nodeDegreeRank

def bootstrap_single(msdata,poplist,distancematrix,bsP):
    pops={}
    for i,p in enumerate(poplist):
        if p not in pops:
            pops[p]=[]
        pops[p].append(i)
    
    selected=[]
    for pop in pops:
        nPop=len(pops[pop])
        bsSize=int(bsP*nPop)
        if bsSize<1:
            bsSize=1
        selected.extend(random.sample(pops[pop],bsSize))
    selected.sort()

    #make new distance matrix
    newmsdata=msdata.getSubset(selected)
    newpoplist=[]
    for index in selected:
        newpoplist.append(poplist[index])
    newglists,newunique_poplist=eden.getGoldsteinLists(newpoplist)
    newdm=newmsdata.getGroupwiseDistanceMatrix(newglists,distancematrix.distancemeasure,groupNames=newunique_poplist)

    return newdm,newmsdata,newpoplist


def autoThreshold(net,outputThreshold=False):
    """
    Returns a thresholded copy of the net given as a parameter. 
    The threshold is determined to be the percolation threshold
    by finding a maximum value for the susceptibility.
    """
    edges=list(net.edges)
    edges.sort(lambda x,y:cmp(x[2],y[2]),reverse=False)
    ee=percolator.EvaluationList(edges)
    ee.setStrengthEvaluations()
    suscmax=0
    suscmax_thresh=0
    nNodes=len(net)
    for cs in percolator.Percolator(ee):
        s=cs.getSusceptibility(nNodes)
        if (s>=suscmax):
            suscmax=s
            suscmax_thresh=cs.threshold
            suscmax_nedges=cs.addedEdges
        
        if cs.getGiantSize()>0.95*len(net):
            break

    newNet=pynet.SymmNet()
    for node in net:
        newNet.addNode(node)

    #continue one threshold level after the threshold maximizing susceptibility
    maxReached=False
    stopAtNext=False
    for e in edges:
        if e[2]==suscmax_thresh:
            maxReached=True
        if maxReached and e[2]!=suscmax_thresh and not stopAtNext:
            stopAtNext=True
            lastW=e[2]
        if stopAtNext and lastW!=e[2]:
            break
        newNet[e[0],e[1]]=e[2]

    newNet.matrixtype=net.matrixtype	

    if outputThreshold:
        return newNet,lastW
    else:
        return newNet
