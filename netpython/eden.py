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
""" Eden module for network toolbox
This module contains Eden specific classes.

"""

import collections
import pynet,random,netext
import communities
from communities import communityTree
from math import sin,cos,asin,sqrt,pi
import math
import numpy
from itertools import *

class EDENException(Exception):
    pass

class ParsingError(EDENException):
    pass

def getGoldsteinLists(poplist):
    """Transforms the list of population indices (poplist) into
    a list of lists, where each internal list contains all indices of
    a population. Useful for getting the Goldstein population-level
    distances for microsatellite data. Outputs: the population
    member list of lists (goldstein_lists), list of unique population labels (uniquepops)"""
    uniquepops=[]
    pops=set()
    for pop in poplist:
        if pop not in pops:
            uniquepops.append(pop)
            pops.add(pop)
    #uniquepops=[ uniq for uniq in poplist if uniq not in locals()['_[1]']]

    goldstein_lists=[]
    for population in uniquepops:
        thislist=[]
        for i,item in enumerate(poplist):
            if item==population:
                thislist.append(i)
        goldstein_lists.append(thislist)
    return [goldstein_lists,uniquepops]

def loadNet_microsatellite(input,removeClones=True,distance="lm"):
    """
    Reads microsatellite data as full net from iterable object
    containing strings (opened file). Look eden.MicrosatelliteData for
    the format the data is read. This function also removes the clones by default.
    """
    msData=MicrosatelliteData(input)
    if removeClones:
        msData=msData.getUniqueSubset()
    msNet=msData.getDistanceMatrix(distance)
    return msNet

class SampleFeatureData(object):
    """ A class for representing data for a set of samples.
    The features can be microsatellites, alleles, presence/absence
    or presence/abundace data.
    """
    def __init__(self,input_iterator,missingValue="999",ploidity=2):
        self.diploidity=ploidity
        self.missingValue=missingValue
        self.parse_input(input_iterator)

    def parse_input(input_iterator):
        #read the data to memory first:
        lines=list(input_iterator)

        #parse the first line to guess the format
        try:
            fields=map(int,lines[0].split())
            self.numeric=True
            missingValue=int(missingValue)
        except ValueError:
            self.numeric=False

        #construct an empty array for the data
        self._alleles=numpy.zeros() #sample, locus, allele



        #--> old code
        self._alleles=[] #The list of locus lists. Each locus list contains alleles as tuples.

        for lineNumber,line in enumerate(input):
            fields=line.split()

            #At the first line, test if data is numerical
            if lineNumber==0:
                try:
                    fields=map(int,fields)
                    self.numeric=True
                    missingValue=int(missingValue)
                except ValueError:
                    self.numeric=False

            if len(fields)%2!=0:
                raise SyntaxError("Input should have even number of columns");
	    elif lastNumberOfFields!=None and lastNumberOfFields!=len(fields):
		raise SyntaxError("The input has inconsistent number of columns")
            else:
                lastNumberOfFields=len(fields)
                if self.numeric:
                    try:
                        fields=map(int,fields)
                    except ValueError:
                        raise SyntaxError("Input contains mixed numeric and not numeric alleles.")

                #At the first line, add lists for loci
                if len(self._alleles)==0:                    
                    for dummy in range(0,len(fields)/2):
                        self._alleles.append([])

                for i in range(0,len(fields),2):                    
                    if fields[i]!=missingValue and fields[i+1]!=missingValue:
                        if fields[i]>fields[i+1]:
                            fields[i],fields[i+1]=fields[i+1],fields[i]
                        self._alleles[i/2].append((fields[i],fields[i+1]))
                    elif fields[i]==missingValue: #None comes first
                        if fields[i+1]==missingValue:
                            self._alleles[i/2].append((None,None))
                        else:
                            self._alleles[i/2].append((None,fields[i+1]))
                    else:
                        self._alleles[i/2].append((None,fields[i]))

        if lastNumberOfFields!=None:
            self.nLoci=lastNumberOfFields/2

    def _clear_cache(self):
        pass


class MicrosatelliteData:
    """ A class for parsing and using microsatellite data
    """
    def __init__(self,input,missingValue="999"):
        """
        The microsatellite data must be given as a input where each row
        has microsatellites for one node/specimen. Alleles should be given as
        integer numbers representing the number of repetitions. The data should have
        even number of columns where each pair represents two homologous alleles.
        Input variable should contain iterable object that outputs the row as a string
        at each iteration step: for example open('msfile.txt').
        """
        self.diploid=True
	lastNumberOfFields=None

        self._alleles=[] #The list of locus lists. Each locus list contains alleles as tuples.

        for lineNumber,line in enumerate(input):
            line=line.strip()
            if len(line)>0:
                fields=line.split()

                #At the first line, test if data is numerical
                if lineNumber==0:
                    try:
                        fields=map(int,fields)
                        self.numeric=True
                        missingValue=int(missingValue)
                    except ValueError:
                        self.numeric=False

                if len(fields)%2!=0:
                    raise SyntaxError("Input should have even number of columns");
                elif lastNumberOfFields!=None and lastNumberOfFields!=len(fields):
                    raise SyntaxError("The input has inconsistent number of columns")
                else:
                    lastNumberOfFields=len(fields)
                    if self.numeric:
                        try:
                            fields=map(int,fields)
                        except ValueError:
                            raise SyntaxError("Input contains mixed numeric and not numeric alleles.")

                    #At the first line, add lists for loci
                    if len(self._alleles)==0:                    
                        for dummy in range(0,len(fields)/2):
                            self._alleles.append([])

                    for i in range(0,len(fields),2):                    
                        if fields[i]!=missingValue and fields[i+1]!=missingValue:
                            if fields[i]>fields[i+1]:
                                fields[i],fields[i+1]=fields[i+1],fields[i]
                            self._alleles[i/2].append((fields[i],fields[i+1]))
                        elif fields[i]==missingValue: #None comes first
                            if fields[i+1]==missingValue:
                                self._alleles[i/2].append((None,None))
                            else:
                                self._alleles[i/2].append((None,fields[i+1]))
                        else:
                            self._alleles[i/2].append((None,fields[i]))

        if lastNumberOfFields!=None:
            self.nLoci=lastNumberOfFields/2

    def copy(self):
        return self.getSubset(range(self.getNumberOfNodes()))
                    
    def shuffleNodes(self):
        """
        Shuffles the order of nodes
        """
        nNodes=len(self._alleles[0])
        nLoci=self.getNumberofLoci()
        for i in range(nNodes):
            r=random.randint(i,nNodes-1)
            for li in range(nLoci):
                tempA=self._alleles[li][i]
                self._alleles[li][i]=self._alleles[li][r]
                self._alleles[li][r]=tempA


    def getNode(self,index):
        """
        Returns a list of alleles where every two alleles in the same loci
        are coupled together with a tuple object.
        """
        node=[]
        for allele in self._alleles:
            node.append(allele[index])
        return tuple(node)

    def getLocusforNodeIndex(self, locus, node):
        return self._alleles[locus][node]

    def getNumberofLoci(self):
        return len(self._alleles)

    def getMSDistance_hybrid(self,x,y,lm_w=None,nsa_w=None):
        if lm_w==None:
            lm_w=self.lm_w
        if nsa_w==None:
            nsa_w=self.nsa_w
        return float(sum(self.getMSDistance_vectorHybrid(x,y,lm_w=lm_w,nsa_w=nsa_w)))

    def getMSDistance_vectorHybrid(self,x,y,lm_w=None,nsa_w=None):
        if lm_w==None:
            lm_w=self.lm_w
        if nsa_w==None:
            nsa_w=self.nsa_w
            distance=numpy.zeros(len(x))
        return nsa_w*self.getMSDistance_vectorNonsharedAlleles(x,y)+lm_w*self.getMSDistance_vectorLinearManhattan(x,y)

    def getGroupwiseDistance_DyerNason(self,x,y):
        """
        Returns the distance between two populations defined by Dyer and Nason in:
        "Population graphs: The graph theoretic shape of genetic structure, Mol Ecol
        13:1713-1727 (2004)". This implementation is coded by following: 
        "M.A. Fortuna et al.: Networks of spatial genetic variation across species, PNAS
        vol. 106 no. 45 pp. 19044-19049 (2009)".

        Parameters
        ----------
        x and y are lists of sample indices correspoding to samples of two populations.
        The distance between these populations is calculated.
        """
        #calculate the centroid multivariate coding vector:
        x_cod={} # the codification vector for population x
        for nodeIndex in x:
            for locus in range(self.getNumberofLoci()):
                alleles=self.getLocusforNodeIndex(locus,nodeIndex)
                key=(locus,alleles[0])
                x_cod[key]+=x_cod.get(key,0)+1
                key=(locus,alleles[1])
                x_cod[key]+=x_cod.get(key,0)+1
        for key in x_cod.keys(): #normalize
            x_cod[key]=x_cod[key]/float(len(x))

                
        raise NotImplementedError()

    def getGroupwiseDistance_Goldstein_D1(self,x,y):
        """
        Returns the goldstein distance between two populations
        This is the ASD (average square distance) of the Goldstein distances.

        Parameters
        ----------
        x and y are lists of sample indices correspoding to samples of two populations.
        The distance between these populations is calculated.

        Example
        -------
        >>> ms=eden.MicrosatelliteData(open("../data/microsatellites/microsatellites.txt",'r'))
        >>> ms_u = ms.getUniqueSubset()
        >>> ms_u.getGroupwiseDistance_Goldstein([1,2,3,4],[1,2,3,4]) == 47.517857142857146
        True
        """
        
        distList = []
        for locus in range(self.getNumberofLoci()):
            
            # calculates allele frequences
            xdict = {}
            for nodeIndex in x:
                xlocus = self.getLocusforNodeIndex(locus,nodeIndex)
                xdict[xlocus[0]] = xdict.get(xlocus[0],0) + 1
                xdict[xlocus[1]] = xdict.get(xlocus[1],0) + 1                
                    
            ydict = {}
            for nodeIndex in y:
                ylocus = self.getLocusforNodeIndex(locus,nodeIndex)
                ydict[ylocus[0]] = ydict.get(ylocus[0],0) + 1 
                ydict[ylocus[1]] = ydict.get(ylocus[1],0) + 1


            dist = 0
            NElementsX=float(sum(xdict.itervalues())-xdict.get(None,0))
            NElementsY=float(sum(ydict.itervalues())-ydict.get(None,0))
            if NElementsX!=0 and NElementsY!=0:
                for i in xdict:
                    if i!=None:
                        for j in ydict:
                            if j!=None:
                                dist += float(i-j)**2*xdict[i]/NElementsX*ydict[j]/NElementsY

                distList.append(dist)

        return sum(distList)/len(distList)

    def getGroupwiseDistance_Goldstein(self,x,y):
        """
        Returns the goldstein distance between two populations. 
        For each allele this is the square of the averages. This function
        returns the average of the values for each allele.

        Parameters
        ----------
        x and y are lists of sample indices correspoding to samples of two populations.
        The distance between these populations is calculated.

        """


        distList=[]
        for locus in range(self.getNumberofLoci()):
            #Calculate the averages
            mx=0.0
            nx=0.0
            for nodeIndex in x:
                xlocus = self.getLocusforNodeIndex(locus,nodeIndex)
                if xlocus[0]!=None:
                    mx+=xlocus[0]
                    nx+=1
                if xlocus[1]!=None:
                    mx+=xlocus[1]
                    nx+=1
            if nx!=0.0:
                mx=mx/float(nx)
            else:
                mx=None

            my=0.0
            ny=0.0
            for nodeIndex in y:
                ylocus = self.getLocusforNodeIndex(locus,nodeIndex)
                if ylocus[0]!=None:
                    my+=ylocus[0]
                    ny+=1
                if ylocus[1]!=None:
                    my+=ylocus[1]
                    ny+=1
            if ny!=0:
                my=my/float(ny)
            else:
                my=None

            if mx != None and my != None:
                distList.append((mx-my)**2)

        if len(distList)!=0:
            return sum(distList)/len(distList)
        else:
            return 0.0

    def getGroupwiseDistanceMatrix(self,groups,distance,groupNames=None):
        """
        Returns a distance matrix in form of a full network (pynet.SymmFullNet). The groups
        argument must be an iterable object where each element is also iterable object containing
        the indices of the nodes belonging to each group.
        """
        distance=distance.lower() #any case is ok
        grouplist=list(groups)
        ngroups=len(grouplist)
        matrix=pynet.SymmFullNet(ngroups)
        if groupNames==None:
            groupNames=range(ngroups)


        if distance in ["goldstein","goldstein_d1"]:
            #only distance measure implemented so far:
            if distance=="goldstein":
                getGroupwiseDistance=self.getGroupwiseDistance_Goldstein
            elif distance=="goldstein_d1":
                getGroupwiseDistance=self.getGroupwiseDistance_Goldstein_D1

            for i in range(0,ngroups):
                for j in range(i+1,ngroups):
                    matrix[groupNames[i],groupNames[j]]=getGroupwiseDistance(grouplist[i],grouplist[j])
            return matrix   
        elif distance in ["fst"]: #allele frequency table based distances
            afTable=AlleleFrequencyTable()
            afTable.init_msData(self,grouplist,groupNames)
            if distance=="fst":
                return afTable.getFST()
        else:
            raise NotImplementedError("Distance '"+distance+"' is not implemented.")
    #--- Distances between individuals

    def getMSDistance_linearManhattan(self,x,y):
        """
        Returns the distance between two nodes/specimen
        """
        return self.getMSDistanceByVector(self.getMSDistance_vectorLinearManhattan(x,y))

    def getMSDistance_nonsharedAlleles(self,x,y):
        return self.getMSDistanceByVector(self.getMSDistance_vectorNonsharedAlleles(x,y))

    def getMSDistance_alleleParsimony(self,x,y):
        return self.getMSDistanceByVector(self.getMSDistance_vectorAlleleParsimony(x,y))


    #-> Distance specific distance vectors
    """
    def getMSDistance_vectorLinearManhattan(self,x,y):
        distance=numpy.zeros(len(x))
        for locus in range(0,len(x)):
            if x[locus]!=None and y[locus]!=None:
                distance[locus]=abs(x[locus][0]-y[locus][0])+abs(x[locus][1]-y[locus][1])
            else:
                distance[locus]=numpy.nan
        return distance
    """
    def getMSDistance_vectorLinearManhattan(self,x,y):
        return self.getMSDistanceVectorByAlleles(x,y,self.getMSDistance_singleLocus_LinearManhattan)
    def getMSDistance_vectorNonsharedAlleles(self,x,y):
        return self.getMSDistanceVectorByAlleles(x,y,self.getMSDistance_singleLocus_NonsharedAlleles)
    def getMSDistance_vectorAlleleParsimony(self,x,y):
        return self.getMSDistanceVectorByAlleles(x,y,self.getMSDistance_singleLocus_AlleleParsimony)
    #<-

    #-> Distance specific single locus distances
    def getMSDistance_singleLocus_NonsharedAlleles(self,x,y):
        alleles=x+y
        distance=0
        for allele in alleles:
            if allele not in x or allele not in y:
                distance+=1
        return distance
            
    def getMSDistance_singleLocus_AlleleParsimony(self,x,y):
        first=len(set([x[0],y[0]]))+len(set([x[1],y[1]]))
        second=len(set([x[0],y[1]]))+len(set([x[1],y[0]]))
        return float(min([first,second]))-2.0

    def getMSDistance_singleLocus_LinearManhattan(self,x,y):
        return abs(x[0]-y[0])+abs(x[1]-y[1])
    #<-

    #->General functions for distance calculations
    def getMSDistanceVectorByAlleles(self,x,y,distance_singleLocus):
        distance=numpy.zeros(len(x))
        for locus in range(0,len(x)):
            if x[locus][0]!=None and y[locus][0]!=None:
                distance[locus]=distance_singleLocus(x[locus],y[locus])
            else:
                distance[locus]=numpy.nan
        return distance

    def getMSDistanceByVector(self,distanceVector):
        size=0
        sum=0
        for i in xrange(len(distanceVector)):
            if distanceVector[i]>=-1: #distanceVector[i]!=numpy.nan
                sum+=distanceVector[i]
                size+=1
        if size!=0:
            return sum/float(size)
        else:
            return -1 #numpy.nan

    def getMSDistanceByAlleles(self,x,y,distance_singleLocus):
        return self.getMSDistanceByVector(self.getMSDistanceVectorByAlleles(x,y,distance_singleLocus))
    #<-
    
    def getDistanceMatrix(self,distance="lm",nodeNames=None,progressUpdater=None):
        """
        Computes the distance between each node and returns the corresponding
        distance matrix.
        """
        if distance=="lm":
            getMSDistance=self.getMSDistance_linearManhattan
        elif distance=="nsa":
            getMSDistance=self.getMSDistance_nonsharedAlleles
        elif distance=="ap":
            getMSDistance=self.getMSDistance_alleleParsimony
        elif distance=="hybrid":
            getMSDistance=self.getMSDistance_hybrid
        elif distance=="czekanowski":
            getMSDistance=self.get_czekanowski_dissimilarity
        else: #default
            getMSDistance=self.getMSDistance_linearManhattan
            
        numberOfSpecimens=len(self._alleles[0])

        j=0
        minUpdateInterval=1000
        minUpdateSteps=30
        totElems=numberOfSpecimens*(numberOfSpecimens-1)/2

        updateInterval=max(min(minUpdateInterval,int(totElems/float(minUpdateSteps))),1)

        elementsAdded=0
        lastUpdate=0

        matrix=pynet.SymmFullNet(numberOfSpecimens)
        if nodeNames==None:
            nodeNames=range(0,numberOfSpecimens)
        for i,iName in enumerate(nodeNames):
            if progressUpdater!=None:
                if elementsAdded-lastUpdate>updateInterval:
                    progressUpdater(float(elementsAdded)/float(totElems))
                    lastUpdate=elementsAdded
                elementsAdded+=numberOfSpecimens-i
            for j in range(i+1,numberOfSpecimens):
                jName=nodeNames[j]
                matrix[iName,jName]=getMSDistance(self.getNode(i),self.getNode(j))
        return matrix
            
    def getSubset(self,nodes):
        """
        Returns a new MicrosatelliteData object containing only nodes given
        as a input. The input is a list of indices of the nodes.
        """
        newData=self.__class__([])
        for allele in self._alleles:
            newAllele=[]
            for node in nodes:
                newAllele.append(allele[node])
            newData._alleles.append(newAllele)
        newData.nLoci=self.nLoci
        return newData

    def randomize(self,full=False):
        """
        Shuffles the homologous alleles in the whole dataset
        """
        if full:
            for j,allele in enumerate(self._alleles):
                newAlleles=[]
                alleleList=[]
                for apair in allele:
                    alleleList.append(apair[0])
                    alleleList.append(apair[1])
                random.shuffle(alleleList)
                for i in range(0,len(alleleList),2):
                    newAlleles.append((alleleList[i],alleleList[i+1]))
                self._alleles[j]=newAlleles

        #this really should shuffle all the homologous alleles
        #and not small and large alleles separately
        else:
            for allele in self._alleles:
                random.shuffle(allele)

    def getNumberOfNodes(self):
        return len(self._alleles[0])
    
    def getUniqueSubset(self,returnOldIndices=False):
        """
        Returns a new MicrosatelliteData object with all identical nodes
        removed except the first occurances of them.
        """
        numberOfNodes=self.getNumberOfNodes()
        nodeSet=set()
        uniqueNodes=[]
        for i in range(0,numberOfNodes):
            node=self.getNode(i)
            if node not in nodeSet:
                nodeSet.add(node)
                uniqueNodes.append(i)
                #print i+1 #indices for Jenni and Jari
        if returnOldIndices:
            return (self.getSubset(uniqueNodes),uniqueNodes)
        else:
            return self.getSubset(uniqueNodes)

    def __str__(self):
        theStr=""
        for nodeIndex in range(self.getNumberOfNodes()):
            node=self.getNode(nodeIndex)
            alleleList=reduce(lambda x,y:x+y,node)
            alleleStr=reduce(lambda x,y:str(x)+" "+str(y),alleleList)
            theStr+=alleleStr+"\n"
        return theStr

class MicrosatelliteDataHaploid(MicrosatelliteData):
    def __init__(self,input,missingValue="999"):
        """
        The microsatellite data must be given as a input where each row
        has microsatellites for one node/specimen. Alleles should be given as
        integer numbers representing the number of repetitions. The data should have
        even number of columns where each pair represents two homologous alleles.
        Input variable should contain iterable object that outputs the row as a string
        at each iteration step: for example open('msfile.txt').
        """
        self.diploid=False
	lastNumberOfFields=None

        self._alleles=[] #The list of locus lists. Each locus list contains alleles as tuples.

        for lineNumber,line in enumerate(input):
            line=line.strip()
            if len(line)>0:
                fields=line.split()

                #At the first line, test if data is numerical
                if lineNumber==0:
                    try:
                        fields=map(int,fields)
                        self.numeric=True
                        missingValue=int(missingValue)
                    except ValueError:
                        self.numeric=False

                if lastNumberOfFields!=None and lastNumberOfFields!=len(fields):
                    raise SyntaxError("The input has inconsistent number of columns")
                else:
                    lastNumberOfFields=len(fields)
                    if self.numeric:
                        try:
                            fields=map(int,fields)
                        except ValueError:
                            raise SyntaxError("Input contains mixed numeric and not numeric alleles.")

                    #At the first line, add lists for loci
                    if len(self._alleles)==0:                    
                        for dummy in range(0,len(fields)):
                            self._alleles.append([])

                    for i in range(0,len(fields)):                    
                        if fields[i]!=missingValue:
                            self._alleles[i].append(fields[i])
                        else:
                            self._alleles[i].append(None)

        if lastNumberOfFields!=None:
            self.nLoci=lastNumberOfFields

    def get_czekanowski_dissimilarity(self,x,y):
        up=0.0
        down=0.0
        for locus_index in range(len(x)):
            up+=min(x[locus_index],y[locus_index])
            down+=x[locus_index]+y[locus_index]
        if down!=0:
            return 1.0-2*up/float(down)
        else:
            return 0.0


    def getMSDistance_singleLocus_NonsharedAlleles(self,x,y):
        if x==y:
            return 0.0
        else:
            return 1.0
            
    def getMSDistance_singleLocus_AlleleParsimony(self,x,y):
        if x==y:
            return 0.0
        else:
            return 1.0

    def getMSDistance_singleLocus_LinearManhattan(self,x,y):
        return abs(x-y)

    def getMSDistanceVectorByAlleles(self,x,y,distance_singleLocus):
        distance=numpy.zeros(len(x))
        for locus in range(0,len(x)):
            if x[locus]!=None and y[locus]!=None:
                distance[locus]=distance_singleLocus(x[locus],y[locus])
            else:
                distance[locus]=numpy.nan
        return distance


    def getGroupwiseDistance_Goldstein(self,x,y):
        """
        Returns the goldstein distance between two populations. 
        For each allele this is the square of the averages. This function
        returns the average of the values for each allele.

        Parameters
        ----------
        x and y are lists of sample indices correspoding to samples of two populations.
        The distance between these populations is calculated.

        """
        distList=[]
        for locus in range(self.getNumberofLoci()):
            #Calculate the averages
            mx=0.0
            nx=0.0
            for nodeIndex in x:
                xlocus = self.getLocusforNodeIndex(locus,nodeIndex)
                if xlocus!=None:
                    mx+=xlocus
                    nx+=1
            if nx!=0.0:
                mx=mx/float(nx)
            else:
                mx=None

            my=0.0
            ny=0.0
            for nodeIndex in y:
                ylocus = self.getLocusforNodeIndex(locus,nodeIndex)
                if ylocus!=None:
                    my+=ylocus
                    ny+=1
            if ny!=0.0:
                my=my/float(ny)
            else:
                my=None

            if mx!=None and my!=None:
                distList.append((mx-my)**2)

        if len(distList)!=0:
            return sum(distList)/len(distList)
        else:
            return 0.0


    def getGroupwiseDistance_Goldstein_D1(self,x,y):
        """
        Returns the goldstein distance between two populations

        Parameters
        ----------
        x and y are lists of sample indices correspoding to samples of two populations.
        The distance between these populations is calculated.
        """
        
        distList = []
        for locus in range(self.getNumberofLoci()):
            
            # calculates allele frequences
            xdict = {}
            for nodeIndex in x:
                xlocus = self.getLocusforNodeIndex(locus,nodeIndex)
                xdict[xlocus] = xdict.get(xlocus,0) + 1
                    
            ydict = {}
            for nodeIndex in y:
                ylocus = self.getLocusforNodeIndex(locus,nodeIndex)
                ydict[ylocus] = ydict.get(ylocus,0) + 1 

            # calculates goldstein distance
            dist = 0
            NElementsX=float(sum(xdict.itervalues())-xdict.get(None,0))
            NElementsY=float(sum(ydict.itervalues())-ydict.get(None,0))
            if NElementsX!=0 and NElementsY!=0:
                for i in xdict:
                    if i!=None:
                        for j in ydict:
                            if j!=None:
                                dist += float(i-j)**2*xdict[i]/NElementsX*ydict[j]/NElementsY

                distList.append(dist)

        return sum(distList)/len(distList)

    def __str__(self):
        raise NotImplemented()

class AlleleDistribution(collections.defaultdict):
    def __init__(self):
        super(AlleleDistribution, self).__init__()
        #self.freqs=collections.defaultdict()
        self.default_factory=lambda:0
        self.totalFrequency=0

    #def __get__(self,item):
    #    return self.freqs[item]
    def __set__(self,item,value):
        self.totalFrequency+=value-self[item]
        super(AlleleDistribution, self).__set__(item,value)

    def get_normaized_distribution(self):
        pass



class AlleleFrequencyTable:
    def init_freqFile(self,filename):
        f=open(filename,'r')
        data=[]
        groups=set()
        loci=set()
        for line in f:
            line=line.strip()
            if not (line.startswith("%") or len(line)==0):
                group,locus,allele=line.split()
                allele = int(allele)
                data.append((group,locus,allele))
                groups.add(group)
                loci.add(locus)
        
        self.nLoci=len(loci)
        self.nGroups=len(groups)
        self.groupNames=sorted(groups)
        locusNames=sorted(loci)
        groupNameToIndex=dict(((g,i) for i,g in enumerate(self.groupNames)))        
        locusNameToIndex=dict(((l,i) for i,l in enumerate(locusNames)))        
        self.freqsTable=[[{} for l in loci] for n in groups]
        for group,locus,allele in data:
            gi=groupNameToIndex[group]
            li=locusNameToIndex[locus]
            self.freqsTable[gi][li][allele]=self.freqsTable[gi][li].get(allele,0)+1

            
    def init_msData(self,msdata,groups,groupNames=None):
        #self.msdata=msdata
        self.nLoci=msdata.getNumberofLoci()
        self.nGroups=len(groups)
        self.freqsTable=[]
        if groupNames==None:
            self.groupNames=range(self.nGroups)
        else:
            self.groupNames=groupNames

        for groupIndex,group in enumerate(groups):
            freqList=[]
            self.freqsTable.append(freqList)
            for locus in range(msdata.getNumberofLoci()):
                #freqs = collections.defaultdict()
                #freqs.default_factory=lambda:0
                freqs={}
                for nodeIndex in group:
                    allele = msdata.getLocusforNodeIndex(locus,nodeIndex)
                    if msdata.diploid and allele!=(None,None):
                        if allele[0]!=None:
                            freqs[allele[0]] = freqs.get(allele[0],0) + 1
                        if allele[1]!=None:
                            freqs[allele[1]] = freqs.get(allele[1],0) + 1
                    elif allele!=None:
                        freqs[allele] = freqs.get(allele,0) + 1
                if len(freqs)==0:
                    raise EDENException("Group %s in input data has only missing values in locus %s."%(self.groupNames[groupIndex],locus))
                freqList.append(freqs)

    def normalizedFreqs(self,group,locus):
        nfreqs = collections.defaultdict()
        nfreqs.default_factory=lambda:0
        tf=float(self.totalFreq(group,locus))
        for key,freq in self.freqsTable[group][locus].iteritems():
            nfreqs[key]=freq/tf
        return nfreqs
    
    def totalFreq(self,group,locus):
        return sum(self.freqsTable[group][locus].itervalues())

    def getFST(self):
        """ From Reynolds, J., Weir, B.S., and Cockerham, C.C. (1983) Estimation of the 
        coancestry coefficient: basis for a short-term genetic distance. _Genetics_, 
        105:767-779, p. 769.
        """
        d=pynet.SymmFullNet(self.nGroups)

        for i in range(self.nGroups):
            for j in range(i):
                num=0.0
                den=0.0
                for locus in range(self.nLoci):
                    ni=float(self.totalFreq(i,locus))
                    nj=float(self.totalFreq(j,locus))
                    if ni>0 and nj>0:
                        summ=0.0
                        ai=1.0
                        aj=1.0
                        fi=self.normalizedFreqs(i,locus)
                        fj=self.normalizedFreqs(j,locus)
                        for l in set(chain(fi.iterkeys(),fj.iterkeys())):                            
                                summ+=(fi[l]-fj[l])**2
                                ai-=fi[l]**2
                                aj-=fj[l]**2
                                        
                        num+=summ/2.-((ni+nj)*(ni*ai+nj*aj))/(4*ni*nj*(ni+nj-1))
                        den+=summ/2.+((4*ni*nj-ni-nj)*(ni*ai+nj*aj))/(4*ni*nj*(ni+nj-1))

                if den>0:
                    d[self.groupNames[i]][self.groupNames[j]] =-math.log(1-num/den)
                else:
                    d[self.groupNames[i]][self.groupNames[j]]=0.0 #not defined

        return d

    def heterozygozity(self,freqs):
        result=1.0
        total=float(sum(freqs.itervalues()))
        for freq in freqs.itervalues():
            f=freq/total
            result-=f*f

        return result
     

    def __get__(self,item):
        return self.freqTable[item]

class BinaryData(object):
    """A class for representing presence/absense data.
    """
    def __init__(self):
        self.data=[]
    def read_file(self,inputfile):
        """Read in and parse the input file. The input is expected to be in a format where
        each row represents one organism/taxa. The inputfile can be given either as a file
        name or any iterable list of strings, e.g. open file.
        """
        def to_bool(element):
            element=element.strip()
            if element in ["0","1"]:
                return bool(int(element))
            else:
                raise ParsingError("Invalid element: "+element+", should be 0 or 1.")

        if isinstance(inputfile,str):
            ifile=open(inputfilen,'rU')
        else:
            ifile=inputfile
        nElements=None
        for i,line in enumerate(ifile):
            if len(line.strip())>0: #skip lines with only whitespaces
                try:
                    elements=map(to_bool,line.split())
                except ParsingError,e:
                    raise ParsingError("Error reading row "+str(i+1)+".\n"+str(e))
                if nElements!=None and len(elements)!=nElements:
                    raise ParsingError("Row %d has %d features while previous row(s) have %d features." % (i+1,len(elements),nElements))
                nElements=len(elements)
                self.data.append(elements)
        if len(self.data)==0 or len(self.data[0])==0:
            raise ParsingError("Error reading data: Empty file.")
    def get_union(self,x,y):
        count=0
        for i in range(len(self.data[0])):
            if self.data[x][i]==True or self.data[y][i]==True:
                count +=1
        return count
    def get_intersection(self,x,y):
        count=0
        for i in range(len(self.data[0])):
            if self.data[x][i]==True and self.data[y][i]==True:
                count +=1
        return count

    def count_true(self,x):
        count=0
        for i in range(len(self.data[x])):
            if self.data[x][i]==True:
                count +=1
        return count       

    def get_bc_dissimilarity(self,x,y):
        """ Bray-Curtis dissimilarity, defined as
        d = 1.0 - 2*intersection(x,y)/(n(x)+n(y))

        If n(x)+n(y) is zero, distance of 1.0 is returned.
        """
        s=self.count_true(x)+self.count_true(y)
        if s>0:
            return 1.0 - 2*self.get_intersection(x,y)/float(s)
        else:
            return 1.0

    def get_jaccard_distance(self,x,y):
        union=float(self.get_union(x,y))
        if union!=0.0:
            return 1-self.get_intersection(x,y)/union
        else:
            return 1.0
    def get_distance_matrix(self,distance_function,node_names,progressUpdater=None):
        size=len(self.data)

        j=0
        updateInterval=1000
        totElems=size*(size-1)/2
        elementsAdded=0
        lastUpdate=0

        distance={"jaccard_distance":self.get_jaccard_distance,
                  "bc_dissimilarity":self.get_bc_dissimilarity}
        matrix=pynet.SymmFullNet(size)
        if node_names==None:
            node_names=range(0,size)
        for i,iName in enumerate(node_names):
            if progressUpdater!=None:
                if elementsAdded-lastUpdate>updateInterval:
                    progressUpdater(float(elementsAdded)/float(totElems))
                elementsAdded+=size-i
            for j in range(i+1,size):
                jName=node_names[j]
                matrix[iName,jName]=distance[distance_function](i,j)
        return matrix


class LocationData:
    def __init__(self,dir="../data/distancematrix_and_locations/"):
        self.location=list(open(dir+"locations.txt"))
        self.name=list(open(dir+"location_names.txt"))
        self.latlong=list(open(dir+"location_coords_latlong.txt"))
        self.c=list(open(dir+"location_classes.txt"))

        #parse latlong:
        newll=[]
        for llstring in self.latlong:
            ll=map(float,llstring.strip().split())
            assert len(ll)==2
            newll.append(ll)
        self.latlong=newll

    def getLocation(self,node):
        return int(self.location[node])

    def getClass(self,node):
        return int(self.c[self.getLocation(node)-1])

    def getName(self,node):
        return self.name[self.getLocation(node)-1]

    def getLatlong(self,node):
        return self.latlong[self.getLocation(node)-1]
        
    def getClassHist(self,nodes):
        locations=map(self.getClass,nodes)
        h=[0,0,0]
        for l in locations:
            h[l-1]=h[l-1]+1;
        return h

    def getClassSets(self,nodes):
        """
        Returns a NodeFamily of the node classes 
        """
        cmap={1:[],2:[],3:[]}
        for node in nodes:
            cmap[self.getClass(node)].append(node)
        return communities.NodeFamily(cmap)

    def getLocationSets(self,nodes):
        """
        Returns a NodeFamily of the node classes 
        """

        cmap={}
        for ln in range(1,len(self.name)+1):
            cmap[ln]=[]
        for node in nodes:
            cmap[self.getLocation(node)].append(node)
        return communities.NodeFamily(cmap)
            
    def getNodesAtLocation(self,location):
        """
        Returns all node indices at the given location index.
        """
        nodes=[]
        for node,nodeLocation in enumerate(self.location):
            if location==int(nodeLocation):
                nodes.append(node)
        return nodes

    def getClassByLocation(self,location):
        """
        Returns class index for the given location.
        """
        if location<1:
            raise ValueError("Invalid location index: " +str(location))
        return int(self.c[location-1])

    def getLatlongByLocation(self,location):
        if location<1:
            raise ValueError("Invalid location index: " +str(location))        
        return self.latlong[location-1]

    def getLocations(self):
        """
        Returns the location indices.
        """
        #return map(lambda x:int(x),self.location)
        return map(lambda x:x+1,range(len(self.name)))

    def getGeoDistMatrix(self,locations=None):
        """
        Returns the matrix of distances between the locations in this object. Locations can be
        speciefied as a list of locations for between which the distance is to be calculated.
        """        
        if locations==None:
            locations=self.getLocations()
        matrix=numpy.zeros([len(locations),len(locations)])
        for i in range(len(locations)):
            for j in range(len(locations)):
                matrix[i,j]=self.getGeoDistByLocation(locations[i],locations[j],lookuptable=False)
        return matrix
    
    def getGeoDistByLocation(self,l1,l2,lookuptable=True):
        if lookuptable:
            if not hasattr(self,"geotable_location"):
                self.geotable=-1*numpy.ones([len(self.latlong),len(self.latlong)])
            d=self.geotable[l1-1,l2-1]
            if d!=-1:
                return d

        #based on code by Jari
        #formula used:   Haversine Formula (from R.W. Sinnott, "Virtues of the Haversine",
        #Sky and Telescope, vol. 68, no. 2, 1984, p. 159):
        #URL: http://www.faqs.org/faqs/geography/infosystems-faq/

        c=(pi/180)
        R=6371000
        ll1=self.getLatlongByLocation(l1)
        ll2=self.getLatlongByLocation(l2)
        lat1,lon1,lat2,lon2=c*ll1[0],c*ll1[1],c*ll2[0],c*ll2[1]

        dlon=lon2-lon1
        dlat=lat2-lat1
        a=pow(sin(dlat/2),2)+cos(lat1)*cos(lat2)*pow(sin(dlon/2),2)
        c=2*asin(min(1,sqrt(a)))
        d=R*c

        if lookuptable:
            self.geotable[l1-1,l2-1]=d
        return d

    def getGeoDist(self,node1,node2,lookuptable=True):
        """
        Calculates the geographic distance of two nodes.
        """
        return self.getGeoDistByLocation(self.getLocation(node1),self.getLocation(node2),lookuptable)
            
class ClassTree(communities.communityTree):
    def getNodeClasses(self,locationData):
        classes={}
        self.cslist.reverse()
        for level in range(0,len(self.cslist)):
            for communityIndex in range(0,len(self.cslist[level])):
                classHist=locationData.getClassHist(self.cslist[level][communityIndex])
                classString=reduce(lambda x,y:str(x)+","+str(y),classHist)
                classes[self._getNameByNode((level,communityIndex))]=classString
        self.cslist.reverse()
        return classes

    def getNodeColors(self,locationData):
        cstrings=self.getNodeClasses(locationData)
        colors={}
        for node in cstrings.keys():            
            c=map(int,cstrings[node].split(','))
            if c[0]>c[1] and c[0]>c[2]:
                color="999900" #yellow
            elif c[1]>c[0] and c[1]>c[2]:
                color="000099" #blue
            else:
                color="990000" #red
            colors[node]=color
        return colors

    def getNodeHist(self,locationData):
        cstrings=self.getNodeClasses(locationData)
        nodeHist={}
        for node in cstrings.keys():            
            c=map(int,cstrings[node].split(','))
            nodeHist[node]=c
        return nodeHist
