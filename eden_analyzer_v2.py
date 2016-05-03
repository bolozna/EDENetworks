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
from RawDataWindow import RawDataWindow
from EdenNetWindow import EdenNetWindow

version="2.18"

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


class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class InputError(Error):
    def __init__(self, message):
        """Exception raised for errors in the input.

        Attributes:
        message -- explanation of the error
        """
        self.message = message


#-----------------#
#                 #
#   MAIN WINDOW   #
#                 #
#-----------------#        
class MainWindow(Toplevel):
    def __init__(self,root):
        
        w=Frame(root,bg='gray80')

        #Add the labels and information about the version
        hellolabel=Label(w,text='Welcome to EDENetworks, the genetic network analyzer!',bg='gray70')
        hellolabel.pack(expand=1,fill=X)
        versionlabel=Label(w,text='version '+version,bg='gray60')
        versionlabel.pack(expand=1,fill=X)
        
        #Add the background image
        imagebody=Frame(w,bg='gray90')
        img=Image.open('infoscreen.jpg')
        self.img=img        
        s=img.size
        canvas=Canvas(imagebody,width=s[0],height=s[1])       
        tkim=ImageTk.PhotoImage(img)
        canvas.create_image(s[0]/2,s[1]/2,image=tkim,anchor=CENTER)
        canvas.tkim=tkim
        canvas.pack()

        imagebody.pack()        
        w.pack()
        self.parent=root

        #Add the menu bar
        menubar=Menu(root)
        filemenu=Menu(menubar,tearoff=0)
        filemenu.add_command(label="Begin new analysis project",command=lambda s=self,r=root: s.launchwizard(r))
        filemenu.add_separator()
        filemenu.add_command(label="Exit",command=self.quitnow)
        menubar.add_cascade(label="File",menu=filemenu)
        root.config(menu=menubar)

        self.networks=[]
        self.netobjs=[]
        self.nnets=0
        self.matrices=[]

    def launchwizard(self,root):
        # first ask if the user wants to read microsatellite repetitions, a distance matrix,
        # or a network (datatype = 'ms','dmat','net')
        dtype=dialogues.ProjectLaunchDialog(root)
        datatype=dtype.result
        if datatype==None:
            return #User pressed cancel.

        # then ask the user for filename (extension will be based on datatype)
        if datatype=='dmat':
            filet=[("Mat files",".mat"),("Text files",".txt")]
        elif datatype in ['ms_haploid',"ms_diploid",'mpop_haploid',"mpop_diploid"] :
            filet=[("Ms files",".ms"),("Text files",".txt")]
        elif datatype=="presabs" or datatype=="presabu":
            filet=[("Mat files",".mat"),("Text files",".txt")]
        else:
            filet=[("Edg files",".edg"),("GML files",".gml")]
        filename=tkFileDialog.askopenfilename(filetypes=filet,title="Step 2: Select your data file")
        if len(filename)==0:
            return #User pressed cancel.

        # to be used for all following loaders
        measuredict={}
        measuredict['nsa']='Non-shared alleles'
        measuredict['lm']='Linear Manhattan'
        measuredict['ap']='Allele parsimony'
        measuredict['hybrid']='Hybrid'
        measuredict['gs']='Goldstein'
        measuredict['jaccard_distance']='Jaccard distance'
        measuredict['czekanowski']='Czekanowski dissimilarity'
        measuredict['bc']='Bray-Curtis dissimilarity'
        measuredict['fst']='FST distance'
	measuredict['other']='Other'

        # next, the file is loaded and depending

        # --------- MICROSATELLITE DATA LOADER ------------------

        if datatype=='mpop_haploid' or datatype=="mpop_diploid":
            #First we test if the inputfile is in a format containing the locations
            #which can be used as a population labels.
            try:
                infileFormat=self.checkMSDataFormat(filename)
                inputfile=open(filename,'rU')
                if infileFormat=="loc_name_ms":
                    nodeLocations,nodeNames,inputfile=self.splitMSFile(inputfile)
                elif infileFormat=="name_ms":
                    nodeNames,inputfile=self.splitMSFile(inputfile,locations=False)
                    nodeLocations=None
                else:
                    nodeLocations,nodeNames=None,None
            except SyntaxError,e:
                tkMessageBox.showerror(
                    "Error reading microsatellites:",
                    "Microsatellite file in wrong format:\n %s" % str(e)
                    )
                return
            poplist=nodeLocations #use node locations as populations

            #Then we read the MS-data.
            try:
                if datatype=="mpop_diploid":
                    msdata=eden.MicrosatelliteData(inputfile)
                else:
                    msdata=eden.MicrosatelliteDataHaploid(inputfile)
            except Exception,e:
                tkMessageBox.showerror(
                    "Error reading microsatellites:",
                    "Microsatellite file in wrong format:\n %s" % str(e)
                    )
                return
            	    
            #If locations/populations are given we do not need to import them
            if nodeLocations==None:
                #Ask for population label file; should contain population labels/indices 
                #for each microsatellite
                filet=[("Text files",".txt")]
                auxfilename=tkFileDialog.askopenfilename(filetypes=filet,
                                                         title="Select population label file")
                if auxfilename==None:
                    return
                poplist=[]
                f=open(auxfilename,'rU')
                #TODO: add sanity checks and error windows
                for line in f:
                    poplist.append(line.strip())


            #check that node names are unique
            try:
                self.checkUniqueNames(nodeNames)
            except Exception,e:
                tkMessageBox.showerror(
                    "Error reading microsatellites:",
                    "Microsatellite file in wrong format:\n %s" % str(e)
                    )
                return



            #If node names are all integers, convert them all from str to int
            try:
                poplist=map(int,poplist)
            except ValueError:
                pass

            #Ask for the distance measure
            distanceMeasureWindow=dialogues.EDENRadioDialog(root,
                                                            title="Choose the distance measure",
                                                            titlemsg="Choose the distance measure",
                                                            options=["Goldstein","FST"])
            distanceMeasure=distanceMeasureWindow.result

            groupList,groupNames=eden.getGoldsteinLists(poplist)
            distancematrix=msdata.getGroupwiseDistanceMatrix(groupList,distanceMeasure,groupNames=groupNames)

            distancematrix.N_origsamples=msdata.getNumberOfNodes()
            distancematrix.distancemeasure=distanceMeasure

            netext.addNodeProperty(distancematrix,"population_size")
            netext.addNodeProperty(distancematrix,"genotypes")
            netext.addNodeProperty(distancematrix,"clonal_diversity")
            for i,item in enumerate(groupNames):
                samples=groupList[i]
                nSamples=len(samples)
                distancematrix.nodeProperty['population_size'][item]=nSamples
                genotypes=msdata.getSubset(samples).getUniqueSubset().getNumberOfNodes()
                distancematrix.nodeProperty['genotypes'][item]=genotypes
                if nSamples>1:
                    R=float(genotypes-1)/float(nSamples-1)
                else:
                    R=0
                distancematrix.nodeProperty['clonal_diversity'][item]=R


            # ask for metadata
            mload=dialogues.ImportMetadataYesNo(root,title="Import auxiliary node data",
                                                titlemsg="Import auxiliary data for sampling sites?",
                                                datatype='msat')
            mloadyesno=mload.result

            # ask for filename
            if mloadyesno!=None and mloadyesno:
                aux_filename=tkFileDialog.askopenfilename(filetypes=[("Text files",".txt"),
                                                                 ("Ascii files",".asc")],
                                                      title="Choose auxiliary data file")
                if len(aux_filename)==0:
                    tkMessageBox.showwarning("Auxiliary data",
                                             "Continuing without importing auxiliary data.")
                else:
                    try:
                        netio.loadNodeProperties(distancematrix,aux_filename)
                    except Exception,e:
                        tkMessageBox.showerror(
                            "Error importing auxiliary data:",
                            "File in wrong format:\n %s" % str(e)
                            )

            tempname=filename.split('/')[-1]
            windowname=tempname.split('.')[0]
            distancematrix.matrixtype=0
            distancematrix.clones='--'
            distancematrix.Nclones=0

           
            RawDataWindow(self.parent,distancematrix,windowname,msdata=msdata,poplist=poplist,nodeTypes='population')

        if datatype=='ms_haploid' or datatype=="ms_diploid":
            #infer the inputfile format and open the input file
            try:
                infileFormat=self.checkMSDataFormat(filename)
                inputfile=open(filename,'rU')
                if infileFormat=="loc_name_ms":
                    nodeLocations,nodeNames,inputfile=self.splitMSFile(inputfile)
                elif infileFormat=="name_ms":
                    nodeNames,inputfile=self.splitMSFile(inputfile,locations=False)
                else:
                    nodeLocations,nodeNames=None,None
            except Exception,e:
                tkMessageBox.showerror(
                    "Error reading microsatellites:",
                    "Microsatellite file in wrong format:\n %s" % str(e)
                    )
                return

            #if only ms-data file given, ask for node names.
            if nodeNames==None:
                nodeNameLoader=dialogues.ImportMetadataYesNo(root,title="Import node names",
                                                    titlemsg="Import node names?",
                                                    datatype='msat')
                loadNodeNames=nodeNameLoader.result

                if loadNodeNames==None: #cancel
                    loadNodeNames=False #maybe should return?
                # ask for filename
                if loadNodeNames:
                    node_filename=tkFileDialog.askopenfilename(filetypes=[("Text files",".txt")],
                                                          title="Choose node names file")
                    if len(node_filename)==0:
                        tkMessageBox.showwarning("Error importing node names",
                                                 "Continuing without importing node names.")
                    else:
                        nodeNameFile=open(node_filename,"rU")
                        nodeNames=map(lambda x:x.strip(),nodeNameFile)

            if nodeNames!=None:
                #check that node names are unique
                try:
                    self.checkUniqueNames(nodeNames)
                except Exception,e:
                    tkMessageBox.showerror(
                        "Error reading microsatellites:",
                        "Microsatellite file in wrong format:\n %s" % str(e)
                        )
                    return

                #If node names are all integers, convert them all from str to int
                try:
                    nodeNames=map(int,nodeNames)
                except ValueError:
                    pass

            # select genetic distance measure first
            mtype=dialogues.ChooseDistanceMeasure(root,title="Choose genetic distance measure",
                                                  titlemsg="Which measure should be used?")
            measuretype=mtype.result
            if measuretype==None:
                return #User pressed cancel
     
            # choose whether to remove clones
            ctype=dialogues.MsatDialog(root,title="Handling clones",
                                       titlemsg="How should clones be handled?")
            removeclones=ctype.result
            if removeclones==None:
                return #User pressed cancel
            
            # load file
            # if clones are removed, because of possible node properties to be read in soon, things have to be done differently here
            # first load the raw microsatellite data, then remove clones (if any), producing indices
            # of kept rows as well
    
            # load the data in a separate window, indicating "please wait"...
            #
            # the logic is as follows: if clones are to be removed AND attributes/metadata added,
            # first load data, generate a distance matrix, and check which rows should be removed
            # as clones. Then, load metadata, add as network attributes, and only after this remove
            # the clonal rows. Simplifies bookkeeping, yet may not be the fastest solution...
            try:                
                ms_result_temp=dialogues.MsLoadWaiter(root,inputfile,removeclones,measuretype,title="Handling microsatellite data",nodeNames=nodeNames,datatype=datatype)
            except SyntaxError,e:
                tkMessageBox.showerror(
                    "Error reading microsatellites:",
                    "Microsatellite file in wrong format:\n %s" % str(e))
                return            
            ms_results=ms_result_temp.result
            Nclones=ms_results[0]
            clones=ms_results[1] #string 'collapsed' or 'included'
            keeptheserows=ms_results[2] #list of rows
            m=ms_results[3] #distance matrix as a fullnet object
            msdata=ms_results[4] #the microsatellite data corresponding to the net
        
            tempname=filename.split('/')[-1]
            windowname=tempname.split('.')[0]
            matrixtype=0 # type 0 : distance matrix

            #if names and locations given to nodes, add locations to the network
            if infileFormat=="loc_name_ms": #this is the only format with location data
                netext.addNodeProperty(m,"location")
                for i,nodeName in enumerate(nodeNames):
                    if nodeName in m:
                        m.nodeProperty["location"][nodeName]=nodeLocations[i]

            # ask for metadata
            mload=dialogues.ImportMetadataYesNo(root,title="Import auxiliary node data",
                                                titlemsg="Import auxiliary data?",
                                                datatype='msat')
            mloadyesno=mload.result


            # ask for filename
            if mloadyesno!=None and mloadyesno:
                aux_filename=tkFileDialog.askopenfilename(filetypes=[("Text files",".txt"),
                                                                     ("Ascii files",".asc")],
                                                          title="Choose auxiliary data file")
                if len(aux_filename)==0:
                    tkMessageBox.showwarning("Auxiliary data",
                                             "Continuing without importing auxiliary data.")
                else:
                    try:
                        netio.loadNodeProperties(m,aux_filename,allowExtraData=True)
                    except Exception,e:
                        tkMessageBox.showerror(
                            "Error importing auxiliary data:",
                            "File in wrong format:\n %s" % str(e)
                            )

            # finally once metadata has (or hasn't been loaded), remove clones. This is done only now so that metadata can be inherited.
            distancematrix=m

            # add informative attributes to distancematrix
            distancematrix.clones=clones
            distancematrix.Nclones=Nclones
            distancematrix.distancemeasure=measuredict[measuretype]
            distancematrix.matrixtype=matrixtype
                                                  
            # open raw data window   
            RawDataWindow(self.parent,distancematrix,windowname,msdata=msdata,nodeTypes='individual')

        elif (datatype=='dmat'):

	    nodenames=dialogues.ChooseMatrixNodeNames(root,title="Please provide node names")	
	    nodenamelist=[]
	    nodenames_result=nodenames.result
	    if nodenames_result=='file':
		nodefilename=tkFileDialog.askopenfilename(filetypes=[("Text files",".txt"),("Ascii files",".asc")],title="Choose node label file")
		if len(nodefilename)==0:
			tkMessageBox.showwarning("Node labels","No file, nodes will be labeled 1..N")
		else:
		   	f=open(nodefilename,'rU')
			nodenamelist=[]
			for line in f:
				nodenamelist.append(line.strip())		
			f.close()
            try:
                if len(nodenamelist)>0:
			inputfile=open(filename,'rU')
                        #If node names are all integers, convert them all from str to int
                        try:
                            nodenamelist=map(int,nodenamelist)
                        except ValueError:
                            pass
                        
                        #Here we can just test reading the matrix with all possible types
                        try:
                            m=netio.loadNet_mat(inputfile,nodeNames=nodenamelist,type="square")
                        except Exception:
                            try:
                                inputfile=open(filename,'rU')
                                m=netio.loadNet_mat(inputfile,nodeNames=nodenamelist,type="upperdiag")
                            except Exception:
                                try:
                                    inputfile=open(filename,'rU')
                                    m=netio.loadNet_mat(inputfile,nodeNames=nodenamelist,type="supperdiag")
                                except Exception:
                                    try:
                                        inputfile=open(filename,'rU')
                                        m=netio.loadNet_mat(inputfile,nodeNames=nodenamelist,type="lowerdiag")
                                    except Exception:
                                        try:
                                            inputfile=open(filename,'rU')
                                            m=netio.loadNet_mat(inputfile,nodeNames=nodenamelist,type="slowerdiag")
                                        except Exception:
                                            tkMessageBox.showerror(
                                                "Error reading distance matrix:",
                                                "Distance matrix file in wrong format.")                                          
                                            return 
                                            
			inputfile.close()
		else:
                        inputfile=open(filename,'rU')
                        for i,line in enumerate(inputfile): pass
                        inputfile.close()
                        inputfile=open(filename,'rU')
                        #Here we cannot know if the triangular formats are strict or not. Square can be tested.
                        try:
                            m=netio.loadNet_mat(inputfile,nodeNames=range(1,i+2))
                            inputfile.close()
                        except Exception:
                            try:
                                inputfile=open(filename,'rU')
                                m=netio.loadNet_mat(inputfile,nodeNames=range(1,i+2),type="upperdiag")
                                if not tkMessageBox.askyesno("Input matrix", "Does the matrix contain diagonal elements?"):
                                    del m
                                    inputfile=open(filename,'rU')
                                    m=netio.loadNet_mat(inputfile,nodeNames=range(1,i+3),type="supperdiag")

                            except Exception:
                                try:
                                    inputfile=open(filename,'rU')
                                    m=netio.loadNet_mat(inputfile,nodeNames=range(1,i+2),type="lowerdiag")
                                    if not tkMessageBox.askyesno("Input matrix", "Does the matrix contain diagonal elements?"):
                                        del m
                                        inputfile=open(filename,'rU')
                                        m=netio.loadNet_mat(inputfile,nodeNames=range(1,i+3),type="slowerdiag")
                                except Exception:
                                    tkMessageBox.showerror(
                                        "Error reading distance matrix:",
                                        "Distance matrix file in wrong format.")                                          
                                    return 
                                    

            except Exception,e:
                tkMessageBox.showerror(
                    "Error reading distance matrix:",
                    "Distance matrix file in wrong format:\n %s" % str(e)
                    )
                return            
                
            # update info        
            tempname=filename.split('/')[-1]
            windowname=tempname.split('.')[0]
            matrixtype=0 # type 0 : distance matrix

            distancematrix=m

            #addinfo=dialogues.LoadMatrixDialog(root,title="Step 4: Please provide additional information")
            #addinfo_result=addinfo.result
	    #if addinfo_result==None: return
            #distancematrix.clones=addinfo_result[1]
            #distancematrix.distancemeasure=measuredict[addinfo_result[0]]

            distancematrix.matrixtype=0 # for distance matrices


            # ask for metadata
            mload=dialogues.ImportMetadataYesNo(root,titlemsg="Import auxiliary node data",datatype='msat')
            mloadyesno=mload.result
            #if mloadyesno==None: return

            if mloadyesno:
                aux_filename=tkFileDialog.askopenfilename(filetypes=[("Text files",".txt"),("Ascii files",".asc")],title="Step 5: Choose auxiliary data file")
                if len(aux_filename)==0:
                    tkMessageBox.showwarning("Auxiliary data",
                                             "Continuing without importing auxiliary data.")
                else:

		    try:
	                    netio.loadNodeProperties(distancematrix,aux_filename)
		    except Exception, e:
			    tkMessageBox.showerror(
				"Error reading node properties:\n %s" % str(e))
			    return	
                                                  
            # open raw data window        
            RawDataWindow(self.parent,distancematrix,windowname,matrixtype=matrixtype)

        elif (datatype=='net'):
            try:
                network=netio.loadNet(filename)
            except Exception,e:
                tkMessageBox.showerror(
                    "Error reading the network:",
                    "File in wrong format:\n %s" % str(e)
                    )
                return            

            mload=dialogues.ImportMetadataYesNo(root,title="Import auxiliary node data",
                                                titlemsg="Import auxiliary data for nodes?",
                                                datatype='msat')
            mloadyesno=mload.result

            # ask for filename
            if mloadyesno!=None and mloadyesno:
                aux_filename=tkFileDialog.askopenfilename(filetypes=[("Text files",".txt"),
                                                                 ("Ascii files",".asc")],
                                                      title="Choose auxiliary data file")
                if len(aux_filename)==0:
                    tkMessageBox.showwarning("Auxiliary data",
                                             "Continuing without importing auxiliary data.")
                else:
                    try:
                        netio.loadNodeProperties(network,aux_filename)
                    except Exception,e:
                        tkMessageBox.showerror(
                            "Error importing auxiliary data:",
                            "File in wrong format:\n %s" % str(e)
                            )
            namestring=os.path.split(filename)[-1].split(".")[0]

	    h=visuals.Himmeli(network)
	    c=h.getCoordinates()	

            EdenNetWindow(self.parent,network,namestring=namestring,coords=c,parentnet=network)

        elif (datatype=='presabs' or datatype=="presabu"):
            measuretype="jaccard_distance"

            #infer the inputfile format and open the input file
            try:
                infileFormat=self.checkMSDataFormat(filename)
                inputfile=open(filename,'rU')
                if infileFormat=="loc_name_ms":
                    nodeLocations,nodeNames,inputfile=self.splitMSFile(inputfile)
                elif infileFormat=="name_ms":
                    nodeNames,inputfile=self.splitMSFile(inputfile,locations=False)
                else:
                    nodeLocations,nodeNames=None,None
            except Exception,e:
                tkMessageBox.showerror(
                    "Error reading the inputfile:",
                    "File in wrong format:\n %s" % str(e)
                    )
                return

            #if only ms-data file given, ask for node names.
            if nodeNames==None:
                nodeNameLoader=dialogues.ImportMetadataYesNo(root,title="Import node names",
                                                    titlemsg="Import node names?",
                                                    datatype='msat')
                loadNodeNames=nodeNameLoader.result

                if loadNodeNames==None: #cancel
                    loadNodeNames=False #maybe should return?
                # ask for filename
                if loadNodeNames:
                    node_filename=tkFileDialog.askopenfilename(filetypes=[("Text files",".txt")],
                                                          title="Choose node names file")
                    if len(node_filename)==0:
                        tkMessageBox.showwarning("Error importing node names",
                                                 "Continuing without importing node names.")
                    else:
                        nodeNameFile=open(node_filename,"rU")
                        nodeNames=map(lambda x:x.strip(),nodeNameFile)

            if nodeNames!=None:
                #check that node names are unique
                try:
                    self.checkUniqueNames(nodeNames)
                except Exception,e:
                    tkMessageBox.showerror(
                        "Error reading microsatellites:",
                        "File in wrong format:\n %s" % str(e)
                        )
                    return

                #If node names are all integers, convert them all from str to int
                try:
                    nodeNames=map(int,nodeNames)
                except ValueError:
                    pass

            try:                
                if datatype=="presabs":
                    data=eden.BinaryData()
                    data.read_file(inputfile)
                    m=dialogues.LoadWaiter(root,data.get_distance_matrix,args=[measuretype,nodeNames]).result
                else: #its presabu
                    data=eden.MicrosatelliteDataHaploid(inputfile)
                    measuretype="czekanowski"
                    m=dialogues.LoadWaiter(root,data.getDistanceMatrix,args=[measuretype,nodeNames]).result
                #m=data.get_distance_matrix(measuretype,nodeNames)
                #ms_result_temp=dialogues.MsLoadWaiter(root,inputfile,removeclones,measuretype,title="Handling microsatellite data",nodeNames=nodeNames,datatype=datatype)
            except Exception,e:
                tkMessageBox.showerror(
                    "Error reading microsatellites:",
                    "File in wrong format:\n %s" % str(e))
                return            
            """
            ms_results=ms_result_temp.result
            Nclones=ms_results[0]
            clones=ms_results[1] #string 'collapsed' or 'included'
            keeptheserows=ms_results[2] #list of rows
            m=ms_results[3] #distance matrix as a fullnet object
            msdata=ms_results[4] #the microsatellite data corresponding to the net
            """

            tempname=filename.split('/')[-1]
            windowname=tempname.split('.')[0]
            matrixtype=0 # type 0 : distance matrix

            #if names and locations given to nodes, add locations to the network
            if infileFormat=="loc_name_ms":
                netext.addNodeProperty(m,"location")
                for i,nodeName in enumerate(nodeNames):
                    if nodeName in m:
                        m.nodeProperty["location"][nodeName]=nodeLocations[i]

            # ask for metadata
            mload=dialogues.ImportMetadataYesNo(root,title="Import auxiliary node data",
                                                titlemsg="Import auxiliary data?",
                                                datatype='msat')
            mloadyesno=mload.result


            # ask for filename
            if mloadyesno!=None and mloadyesno:
                aux_filename=tkFileDialog.askopenfilename(filetypes=[("Text files",".txt"),
                                                                     ("Ascii files",".asc")],
                                                          title="Choose auxiliary data file")
                if len(aux_filename)==0:
                    tkMessageBox.showwarning("Auxiliary data",
                                             "Continuing without importing auxiliary data.")
                else:
                    try:
                        netio.loadNodeProperties(m,aux_filename,allowExtraData=True)
                    except Exception,e:
                        tkMessageBox.showerror(
                            "Error importing auxiliary data:",
                            "File in wrong format:\n %s" % str(e)
                            )

            # finally once metadata has (or hasn't been loaded), remove clones. This is done only now so that metadata can be inherited.
            distancematrix=m

            # add informative attributes to distancematrix
            #distancematrix.clones=clones
            #distancematrix.Nclones=Nclones
            distancematrix.distancemeasure=measuredict[measuretype]
            distancematrix.matrixtype=matrixtype
                                                  
            # open raw data window   
            RawDataWindow(self.parent,distancematrix,windowname,nodeTypes='individual')



    def displayBusyCursor(self):
        self.parent.configure(cursor='watch')
        self.parent.update()
        self.parent.after_idle(self.removeBusyCursor)

    def removeBusyCursor(self):
        self.parent.configure(cursor='pencil')
        
    def readnet(self):
        netfilename=tkFileDialog.askopenfilename(filetypes=[("Edg files", ".edg"),("GML files",".gml")])
        print "Loading", netfilename
        
        self.nnets+=1

        netname_temp=netfilename.split('/')[-1]
        netname=netname_temp.split('.')[0]
        self.networks.append(netname)
        n=netio.loadNet(netfilename)

        self.netobjs.append(NetWindow(self.parent,n,self.networks[self.nnets-1]))

    def readmatrix(self,root):
        matrixfilename=tkFileDialog.askopenfilename(filetypes=[("Mat files",".mat")])

        # check if weight or distance matrix
        # types: 1 = weight matrix, 0 = distance matrix
        
        m=dialogues.MatrixDialog(root)
        matrixtype=m.result

        print "Loading", matrixfilename
       # print "Matrixtype",matrixtype

        self.nmatrices+=1
        tempname=matrixfilename.split('/')[-1]
        self.matrices.append(tempname.split('.')[0])
        m=netio.loadNet(matrixfilename)
        self.matobjs.append(MatrixWindow(self.parent,m,self.matrices[self.nmatrices-1],matrixtype=matrixtype))

    def readmsat(self,root):
        matrixfilename=tkFileDialog.askopenfilename()
        m=dialogues.MsatDialog(root)
        removeclones=m.result

        inputfile=open(matrixfilename)

        m=netio.loadNet_microsatellite(inputfile,removeClones=removeclones)

        inputfile.close()
        
        self.nmatrices+=1
        tempname=matrixfilename.split('/')[-1]
        self.matrices.append(tempname.split('.')[0])
        matrixtype=0
        
        self.matobjs.append(MatrixWindow(self.parent,m,self.matrices[self.nmatrices-1],matrixtype=matrixtype))
        
    def quitnow(self):
        self.parent.destroy()

    def dummy(self):
        pass

    def checkMSDataFormat(self,infilename):
        """
        This function checks if the given file is a flat microsatellite list or
        if it contains the sampling location in the first column, sample name in the
        second column and MS marker lengths in rest of the columns.

        Parameters
        ----------
        infilename : string

        Returns
        -------
        A string containing the result of the check:
        'ms' : Only microsatellites
        'name_ms' : names and microsatellites
        'loc_name_ms' : locations names and microsatellites


        Exceptions
        ----------
        If the file is in invalid format, a SyntaxError exception is thrown, and IOError
        if the file does not exist.
        """
        f=open(infilename,'rU')
        nColumns=None
        ftype=None
        i=None
        for i,line in enumerate(f):
            fields=line.strip().split()

            #check that the number of columns is constant
            if nColumns==None:
                nColumns=len(fields)
            if nColumns!=len(fields):
                raise SyntaxError("The number of columns in line %s differs from the previous line" % str(i+1))

            #check for type of data
            if reduce(lambda x,y:x and y,map(str.isdigit,fields)): #all fields int
                line_ftype='ms'
            elif len(fields)>1 and reduce(lambda x,y:x and y,map(str.isdigit,fields[1:])):
                line_ftype='name_ms'
            elif len(fields)>2 and reduce(lambda x,y:x and y,map(str.isdigit,fields[2:])):
                line_ftype='loc_name_ms'
            else:
                raise SyntaxError("At line "+str(i+1)+": All columns after the second column should be integers.")

            #check that each line is of same type
            if ftype==None:
                ftype=line_ftype
            if ftype!=line_ftype:
                if ftype=="ms":
                    ftype=line_ftype
                elif ftype=="name_ms" and line_ftype!="ms":
                    ftype=line_ftype

        if i==None:
            raise SyntaxError("Input file is empty.")

        return ftype

    def splitMSFile(self,infile,locations=True):
        """
        Reads the inputfile and splits it into three lists. 
        First two list contains two first columns and the third list contains the rest
        of the columns as a string. That is, first two lists contain locations and node names 
        and the third contains the MS-data.
        """
        locs=[]
        names=[]
        msdata=[]
        if locations:
            for line in infile:
                fields=line.split()
                locs.append(fields[0])
                names.append(fields[1])
                msdata.append(" ".join(fields[2:]))
            return locs,names,msdata
        else:
            for line in infile:
                fields=line.split()
                names.append(fields[0])
                msdata.append(" ".join(fields[1:]))
            return names,msdata


    def checkUniqueNames(self,namelist):
        """
        This function checks that the node names are unique. It raises an Exception
        if the node names are not unique. The exeception contains information
        about the two rows which are not unique.
        """
        nameMap={} #map from node name to its row number
        for i,name in enumerate(namelist):
            if name not in nameMap:
                nameMap[name]=i
            else:
                raise Exception("Node names not unique: '"+name+"' in lines "+str(nameMap[name]+1)+" and "+str(i+1))
