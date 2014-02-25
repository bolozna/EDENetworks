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

Network analysis functions and some related dialog windows for pynet

"""


import pynet,os,netio,netext
import random
import heapq
import string
import percolator
import shutil
from math import ceil
from Tkinter import *

def generateLogbins(minvalue,maxvalue,factor,uselinear=True):
    '''Generates a binning vector containing bin limits
       for log-binning. Inputs: min and max values to be binned,
       factor for increasing bin size, uselinear=[True|False] for
       making the first 10 bins linear.'''

    if uselinear:
    
        bins=[-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5]
        i=12

    else:

        bins=[]
        bins.append(minvalue*2.0/(1+factor))
        i=1

    while bins[i-1]<maxvalue:
        bins.append(bins[i-1]*factor)
        i+=1
                    
    return bins

def calculateNk(network):
    '''Calculates # of nodes of degree k (N(k)) for a
       network. Input: network, output: N(k) as list'''
    
    d=netext.deg(network)
    degrees=d.values()

    Nk=[]
    for i in range(max(degrees)+1):

        Nk.append(0)

    for i in range(len(degrees)):
        Nk[degrees[i]]+=1

    return Nk

def cumulativePk(Nk):
    '''Calculates the cumulative degree distribution.
    Input: list containing N(k), i.e. # of nodes with
    degree k. Output: P>(k) as a list, indices corresponding
    to degrees (Pk[0] = P(k>0) = 1, etc'''

    Pk=[]
    N=sum(Nk)
    for i in range(len(Nk)-1):
        Pk.append(sum(Nk[(i+1):])/float(N))
    return Pk

def degreeAverages(degrees,Nk,values):
    '''Calculates the average of some quantity, averaged
    over degree. In: list of degrees, list of N(k), values
    to be averaged. Out: [degrees, avg_values]'''

    dsum=[]
    for i in range(len(Nk)):
        dsum.append(0.0)

    for i in range(len(values)):
      #  print i,degrees[i],type(degrees[i]),values[i],type(values[i])

        if not(values[i]==None):
            dsum[degrees[i]]+=values[i]
        else:
            dsum[degrees[i]]=None

    v_of_k=[]
    k=[]

    for i in range(len(Nk)):
        if not(Nk[i]==0) and not(dsum[i]==None):
            v_of_k.append(dsum[i]/float(Nk[i]))
            k.append(i)

    return [k,v_of_k]
                          

def binAverages(bins,bincenters,degrees,sumdata):
    '''returns averages of sumdata in degree bins.
       Inputs: list of bin lower limits, list of bin
       centers, list of degrees, data to be averaged'''

    Nbin=[]
    Sbin=[]
    for i in range(len(bins)-1):
        Nbin.append(0)
        Sbin.append(0)

    Nk=[]
    Sk=[]
    for i in range(max(degrees)+1):
        Nk.append(0)
        Sk.append(0)

    for i in range(len(degrees)):
        if not(sumdata[i]==None):
            Nk[degrees[i]]+=1
            Sk[degrees[i]]+=sumdata[i]

    bin=0
    for k in range(len(Nk)):
        if k>bins[bin+1]:
            bin+=1
        Nbin[bin]+=Nk[k]
        Sbin[bin]+=Sk[k]

    kvector=[]
    sumvector=[]

    for i in range(len(Nbin)):

        if (Nbin[i]>0):
            kvector.append(bincenters[i])
            sumvector.append(Sbin[i]/float(Nbin[i]))

    return [kvector,sumvector]     

def nodelevelKnn(network,weighted=False):
    '''Calculates average nearest neigh degree for
    all nodes in a network. Returns list of knn:s.
    If weighted is set to True, the avg nn degree
    is weighted by edge weights.'''

    knni=[]
    degs=[]

    for i in network:
        ksum=0.0
        for j in network[i]:
            if weighted:
                ksum+=float(network[i][j])*network[j].deg()
            else:
                ksum+=network[j].deg()

        currdeg=float(network[i].deg())
        if currdeg>0:
            if weighted:
                knni.append(ksum/float(network[i].strength()))
            else:
                knni.append(ksum/float(network[i].deg()))
            
        else:
            knni.append(0.0)
        degs.append(network[i].deg())

    return [degs,knni]
                    

def generateLinbins(minvalue,maxvalue,no_bins):
    '''Generates no_bins of even width. First
    bin centered around minvalue, last bin
    centered around maxvalue.'''

    valuerange=maxvalue-minvalue
    binwidth=valuerange/float(no_bins-1)

    bins=[]
    bins.append(minvalue-binwidth/2.0)

    currvalue=bins[0]

    while currvalue<maxvalue:
        currvalue=currvalue+binwidth
        bins.append(currvalue)

    return bins

def binDensity(bins,Nk):
    '''Returns density in bins (Nk per bin divided by bin width)'''
    Nbin=[]
    for i in range(len(bins)-1):
        Nbin.append(0)
    bin=0

    for k in range(len(Nk)):
        if k>bins[bin+1]:
            bin+=1
        Nbin[bin]+=Nk[k]

    for i in range(len(bins)-1):
        width=float(bins[i+1]-bins[i])
        Nbin[i]=Nbin[i]/width

    return Nbin

def binCounts(bins,Nk):
    '''Returns count Nk in bins'''
    Nbin=[]
    for i in range(len(bins)-1):
        Nbin.append(0)
    bin=0

    for k in range(len(Nk)):
        if k>bins[bin+1]:
            bin+=1
        Nbin[bin]+=Nk[k]

    return Nbin


def binCenters(bins):
    '''Returns bin centers for a list of bin lower limits'''
    
    bincenters=[]
    for i in range(len(bins)-1):
        bincenters.append(0.5*(bins[i+1]+bins[i]))
    return bincenters

def logPk(network,binfactor):

    Nk=calculateNk(network)
    bins=generateLogbins(1.0,len(Nk),binfactor)

    Pbin=binCounts(bins,Nk)
    bincenters=binCenters(bins)

    return [Pbin,bincenters]

def normalize(Nk):

    normalization=float(sum(Nk))

    for i in range(len(Nk)):
        Nk[i]=Nk[i]/normalization

    return Nk

def clustering_valuelist(net):
    '''Returns a list [k_i, c_i] for each node,
       where k_i is its degree and c_i the clustering coeff'''
    
    c=[]
    degs=[]
    for i in net:
        tempc=0
        for j in net[i]:
            for k in net[j]:
                if k in net[i]:
                    tempc+=1
        k=net[i].deg()
        if k>1:
           c.append(float(tempc)/float(k*(k-1)))
           degs.append(k)
        else:
            c.append(None)
            degs.append(k)
    return [degs,c]

def weight_distribution(network,style='logbin',Nbins=25):
    '''Returns the binned weight probability distribution of a network.
       Optional inputs: style = 'logbin' or 'linbin', Nbins = # of bins
       Output: list [w,P(w)]'''

    witer=network.weights.__iter__()
    weight_vector=[]
    for eachitem in witer:
        weight_vector.append(eachitem)

    minw=min(weight_vector)
    maxw=max(weight_vector)

    if style=='logbin':

        factor=(float(maxw)/minw)**(1.0/Nbins)
        bins=generateLogbins(minw,maxw,factor,False)

    else:

        bins=generateLinbins(minw,maxw,Nbins)

    bc=binCenters(bins)
    Pbin=probabilityLogbinner(bins,weight_vector)

    return [bc,Pbin]

def strength_distribution(network,style='logbin',Nbins=25):
    '''Returns the binned strength probability distribution of a network.
    Optional inputs: style = 'logbin' or 'linbin', Nbins = # of bins
    Output: list [s,P(s)]'''

    sv=netext.strengths(network)
    strength_vector=sv.values()
    minw=min(strength_vector)
    maxw=max(strength_vector)

    if style=='logbin':

        factor=(float(maxw)/minw)**(1.0/Nbins)
        bins=generateLogbins(minw,maxw,factor,False)

    else:

        bins=generateLinbins(minw,maxw,Nbins)

    bc=binCenters(bins)
    Pbin=probabilityLogbinner(bins,strength_vector)

    return [bc,Pbin]

def knn_spectrum(network,style='logbin',weighted=False,Nbins=25):
    '''Returns (binned) degree vs average nearest neighbour degree.
    Optional inputs: style='logbin' or 'linbin', Nbins = # of bins,
    weighted = False or True (True: neigh degrees weighted by link weights)
    Output: list [k,knn(k)]'''

    knni=nodelevelKnn(network,weighted) # returns 2 lists [degrees,knn]
    d=netext.deg(network)
    degrees=d.values()
    maxk=max(degrees)

    if style=='logbin':

        factor=(maxk/10.5)**(1/float(Nbins))

        bins=generateLogbins(1.0,maxk,factor)
        bc=binCenters(bins)

        temp=binAverages(bins,bc,knni[0],knni[1])

    else:

        Nk=calculateNk(network)
        temp=degreeAverages(knni[0],Nk,knni[1])

    return temp

def degree_distribution(network,style='logbin',Nbins=25.0):
    '''Calculates degree distribution. Inputs: network,
       style=('nobin'|'logbin'|'cumulative'), Nbins=number of logbins.
       Output: [k,P(k)]'''

    if style=='logbin':

        d=netext.deg(network)
        degs=d.values()
        maxk=max(degs)
        factor=(maxk/10.5)**(1/float(Nbins))

        bins=generateLogbins(1.0,maxk,factor)
        Nk=calculateNk(network)
        Nbin=binDensity(bins,Nk)
        bc=binCenters(bins)

        NNbin=normalize(Nbin)

        temp=[bc,NNbin]

    elif style=='cumulative':

        Nk=calculateNk(network)
        Pk=cumulativePk(Nk)
        k=range(len(Pk))

        temp=[k,Pk]

    else:
        
        Nk=calculateNk(network)
        k=range(len(Nk))

        NNk=normalize(Nk)

        temp=[k,NNk]


    return temp
        

def clustering_spectrum(network,style='nobin',Nbins=25):
    '''Calculates clustering coeff as function of degree.
       Inputs: network, style=('nobin'|'logbin'),Nbins=number of logbins
       Outputs: [k,c(k)]'''

    if style=='logbin':
    
        d=netext.deg(network)
        degrees=d.values()
        maxk=max(degrees)
        cvalues=clustering_valuelist(network)
        factor=(maxk/10.5)**(1.0/Nbins)

        bins=generateLogbins(1.0,maxk,factor)
        bc=binCenters(bins)

        temp=binAverages(bins,bc,cvalues[0],cvalues[1])

    else:

        Nk=calculateNk(network)
      
        cvalues=clustering_valuelist(network) # two lists - degrees & c-coeffs
        temp=degreeAverages(cvalues[0],Nk,cvalues[1])

    return temp

def probabilityLogbinner(bins,valuevector):
    '''Generic log binner. Inputs: bin limit vector,
    vector of values. Counts # of values in each bin,
    divides by bin width, and normalizes to unit sum.
    Returns normalized probability density per bin.'''

    Nbin=[]
    for i in range(len(bins)-1):
        Nbin.append(0)

    for w in valuevector:
        mybin=findbin(bins,w)
        Nbin[mybin]+=1

    for i in range(len(bins)-1):
        width=float(bins[i+1]-bins[i])
        Nbin[i]=Nbin[i]/width

    Nbin=normalize(Nbin)

    return Nbin
        

def findbin(bins,value):
    '''Finds out the bin where input param value belongs to.
    Uses halving for fast output. Inputs: bins - list of bin limits,
    value - value to be found within bin limits'''

    lowerlimit=0
    upperlimit=len(bins)-1

    while (upperlimit-lowerlimit)>1:

        halfpoint=int(ceil(0.5*(upperlimit+lowerlimit)))

        if (value>=bins[halfpoint]):

            lowerlimit=halfpoint

        else:

            upperlimit=halfpoint

    return lowerlimit

    
def clustering(net):
    c={}
    for i in net:
        c[i]=0
        for j in net[i]:
            for k in net[j]:
                if k in net[i]:
                    c[i]=c[i]+1
        k=net[i].deg()
        if k>1:
            c[i]=float(c[i])/float(k*(k-1))
        else:
            c[i]=0
    return c

def globalClustering(net):
    c=clustering(net)
    if len(c)!=0:
        return float(sum(c.values()))/float(len(c))
    else:
        return 0
        

def weightStats(net):
    """Input: network, output: [min_weight,max_weight,avg_weight]"""
    weight_vector = [w for w in net.weights]
    maxw = max(weight_vector)
    minw = min(weight_vector)
    avgw = sum(weight_vector)/len(weight_vector)
    return minw, maxw, avgw
    
def overlap(net,node1,node2):
    """
    Returns the overlap of the edge between the two nodes 
    given as input. Overlap is defined as:
    n_ij/(k_i-1+k_j-1-n_ij)
    where n_ij is the number of common neighbors of nodes i and j
    (=number of triangles) and k_i and k_j are the degrees of nodes i
    and j.
    """
    nTriangles=0
    if net[node1].deg()>net[node2].deg():
        small,large=node2,node1
    else:
        small,large=node1,node2
    for neigh in net[small]:
        if large in net[neigh]: #assume no self-links.
            nTriangles+=1
    d=net[node1].deg()+net[node2].deg()-2-nTriangles
    if d>0:
        return nTriangles/float(d)
    else:
        return 0.0

def edgeClustering(net,node1,node2):
    """
    Returns the edge-clustering of the edge between the two nodes 
    given as input. Edge clustering is defined as:
    n_ij/(min(k_i,k_j)-1)
    where n_ij is the number of common neighbors of nodes i and j
    (=number of triangles) and k_i and k_j are the degrees of nodes i
    and j.
    In case min(k_i,k_j)=1; we define it as 0.
    """
    nTriangles=0
    if net[node1].deg()>net[node2].deg():
        small,large=node2,node1
    else:
        small,large=node1,node2
    for neigh in net[small]:
        if large in net[neigh]: #assume no self-links.
            nTriangles+=1
    d=net[small].deg()-1.0
    if d>0:
        return nTriangles/float(d)
    else:
        return 0.0
