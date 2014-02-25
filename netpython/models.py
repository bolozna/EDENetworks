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
import pynet,netext
import numpy as np

def makeER(n,p):
    """
    Make a realisation of Erdos-Renyi network
    * fast for non-sparse networks
    * the sparce version should be implemented
    """
    net=pynet.SymmNet()
    for i in range(0,n):
        for j in range(0,i):
            if p > np.random.ranf():
                net[i,j]=1
    return net

def makeSparseER(n,p):
    """
    Make a realisation of Erdos-Renyi network
    * fast for sparse networks
    * 0 < p < 1
    * Algorithm: 
    Efficient generation of large random networks
    Phys. Rev. E 71, 036113 (2005) 
    """

    net=pynet.SymmNet()
    v = 1 
    w = -1
    while (v < n):
        r = np.random.ranf()
        w = w+1+int(np.floor(np.log(1-r)/np.log(1-p)))
        while ((w >= v) and (v < n)):
            w = w-v
            v = v+1
        if (v < n):
            net[v,w]=1

    return net

def linearLattice(n,r):
    """Linear lattice with periodic boundary conditions. Two nodes are connected
    if they are at most r steps away from each other in the lattice.
    """
    net=pynet.SymmNet()
    if r>=n:
        r=n-1
    for i in range(n):
        for ri in range(r):
            net[i,(i+ri+1)%n]=1
            net[i,(i-1-ri)%n]=1
    return net

def girvanNewman(communitySize,numberOfCommunities,kIn,kOut):
    """
    A network model producing equally sized communities with equal expected 
    link density inside the communities and between the communities. The model
    was first defined in the article:
    M. Girvan and M.E.J. Newman: Community structure in social and biological networks,
    PNAS 99, 7821 (2002)

    Parameters
    ----------
    communitySize : int 
      Size of a single community in nodes.
    numberOfCommunities : int
      Number of communities
    kIn : int
      Expected value for inside community node degrees. That is, the expected
      number of links from each node going to other nodes in the same community. This
      parameter is used to calculate the probability of links inside communities. If
      kIn > communitySize-1, then kIn is set to communitySize-1.
    kOut : int
      Expected value for outside community node degrees. That is, the expected
      number of links from each node going to nodes in other communities. This
      parameter is used to calculate the probability of links between communities. If
      kOut > (numberOfCommunities-1)*communitySize, then kOut is set to (numberOfCommunities-1)*communitySize.

    Return
    ------
    net : SymmNet 
      A realisation of the model network.
      
    Complexity
    ----------
    For a network with N nodes:
    Time complexity: O(N**2)
    Memory complexity: Memory used by the returned network object. 

    Time complexity can be improved for sparse networks.

    """
    
    #Calculate pIn and pOut from kIn and kOut
    if (communitySize-1)<kIn:
        kIn=communitySize-1 
    if (numberOfCommunities-1)*communitySize<kOut:
        kOut=(numberOfCommunities-1)*communitySize
    pIn=float(kIn)/float(communitySize-1)
    if numberOfCommunities>1:
        pOut=float(kOut)/float((numberOfCommunities-1)*communitySize)        
    else:
        pOut=0.0

    net=pynet.SymmNet() #the net object to be returned
    
    #First, put the internal edges:
    for communityIndex in range(numberOfCommunities):
        for node1Index in range(communitySize):
            for node2Index in range(node1Index+1,communitySize):
                if pIn > np.random.ranf():
                    net[communityIndex*communitySize+node1Index,communityIndex*communitySize+node2Index]=1

    #Second, put the external edges:
    for community1Index in range(numberOfCommunities):
        for community2Index in range(community1Index+1,numberOfCommunities):
            for node1Index in range(communitySize):
                for node2Index in range(communitySize):
                    if pOut > np.random.ranf():
                        net[community1Index*communitySize+node1Index,community2Index*communitySize+node2Index]=1

    return net
