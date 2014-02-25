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
""" A script for testing the speed of EDENetworks backend.
"""
try:
    import numpypy as numpy#for pypy
except ImportError:
    import numpy


from timeit import Timer
import pickle
import random
from netpython import eden,percolator,netext,pynet,netanalysis,transforms
import math,time

posidonia="exampleData/posidonia.ms"

def median(l):
    sl=sorted(l)
    n=len(l)
    if len(n)==1:
        return n[0]
    if n%2==0:
        return (sl[n/2]+sl[n/2-1])/2.
    else:
        return sl[math.floor(n/2)]


############ Code copied from various places #################
# All the code here is copied from places where we shouldn't #
# have it in the first place. It's there due to history of   #
# the program. This issue should be fixed in the next major  #
# revision round.                                            #
##############################################################
def splitMSFile(infile,locations=True):
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


    if outputThreshold:
        return newNet,lastW
    else:
        return newNet

############ End of copied code  #################


def create_random_data(groups=10,nodespergroup=10,loci=10,alleletype=200):
    """Create random allele data.
    """
    data=[]
    for g in range(groups):
        for n in range(nodespergroup):
            row=[]
            row.append(str(g))
            row.append(str(g*nodespergroup+n))
            for l in range(loci):
                if isinstance(alleletype,int):
                    row.append(str(random.randint(0,alleletype)))
                else:
                    raise Exception("Invalid allele type.")
            data.append(" ".join(row))
    return data

class TestDistances(object):
    tests_slow=["get_dm_groups",
                "get_dm_ind",
                "get_globalClustering_ind",
                "auto_threshold_ind"
                ]
    tests_fast=["shuffle_alleles",
                "remove_clones",
                "get_mst_ind"
                ]
    tests_groups=["auto_threshold_groups",
                "get_globalClustering_groups",
                "get_bc_groups",
                "get_mst_groups",
                ]




    test_to_name={"auto_threshold_groups" : "Automatic thresholding",
                  "get_globalClustering_groups" : "Clustering coefficient (groups)",
                  "get_bc_groups" : "Betweenness centrality (groups)",
                  "get_mst_groups" : "Minimal spanning tree (groups)",
                  "get_dm_ind" : "Distance matrix (samples)",
                  "auto_threshold_ind" : "Automatic thresholding",
                  "get_globalClustering_ind" : "Clustering coefficient",
                  "get_mst_ind" : "Minimal spanning tree",
                  "shuffle_alleles" : "Shuffle alleles",
                  "get_dm_groups" : "Distance matrix (groups)",
                  "remove_clones" : "Remove clones"
        }

    def __init__(self,data,groups=True,ind_distance="lm",group_distance="goldstein_d1",ploidity=2):
        self.tests=self.tests_slow+self.tests_fast
        if isinstance(data,str) or isinstance(data,unicode):
            data=open(data,'r')

        self.ind_distance=ind_distance
        self.group_distance=group_distance
            
        if groups:
            locs,names,data=splitMSFile(data)
            self.data=data
            self.locs=locs
            self.names=names
            self.pops,names=eden.getGoldsteinLists(locs)
        else:
            names,data=splitMSFile(data,locations=False)
            self.data=data
            self.names=names

        if ploidity==2:
            self.ms=eden.MicrosatelliteData(self.data)
        elif ploidity==1:
            self.ms=eden.MicrosatelliteDataHaploid(self.data)
        else:
            raise Exception("invalid ploidity")

        self.ms_noclones=self.ms.getUniqueSubset()
        if groups:
            self.dm_groups=self.get_dm_groups()
            self.thnet_groups=self.auto_threshold_groups()
        self.dm_ind=self.get_dm_ind()
        self.thnet_ind=self.auto_threshold_ind()


    def shuffle_alleles(self):
        cp=self.ms.copy()
        cp.randomize()
        return cp

    def get_dm_groups(self,distance=None):
        if distance==None:
            distance=self.group_distance
        return self.ms.getGroupwiseDistanceMatrix(groups=self.pops,distance=distance)

    def get_dm_ind(self,distance=None):
        if distance==None:
            distance=self.ind_distance
        return self.ms_noclones.getDistanceMatrix(distance=distance)

    def remove_clones(self):
        return self.ms.getUniqueSubset()

    def auto_threshold_ind(self):
        return autoThreshold(self.dm_ind,outputThreshold=False)

    def auto_threshold_groups(self):
        return autoThreshold(self.dm_groups,outputThreshold=False)

    def get_globalClustering_groups(self):
        return netanalysis.globalClustering(self.thnet_groups)

    def get_globalClustering_ind(self):
        return netanalysis.globalClustering(self.thnet_ind)

    def get_bc_groups(self):
        return netext.getBetweennessCentrality(self.thnet_groups)

    def get_mst_groups(self):
        return transforms.mst(self.thnet_groups)

    def get_mst_ind(self):
        return transforms.mst(self.thnet_ind)

    def test_all(self,iterations):
        return self.test_list(iterations,self.tests)

    def test_list(self,iterations,thelist):
        results=[]
        time.sleep(1)
        for test in thelist:
            testf=eval("self."+test)
            timer=Timer(testf)
            result=min(timer.repeat(iterations,number=1))
            results.append(result)
            print result
        return results

def run_tests():
    results={}

    iterations=10
    test=TestDistances(posidonia)
    results["posidonia"]=test.test_all(iterations)


    test=TestDistances(create_random_data(groups=300,nodespergroup=10,loci=10,alleletype=200))
    test.test_list(1000,test.tests_fast)
    test.test_list(10,test.tests_slow)

    iterations=1000
    ngroups=range(100,301,50)
    ngroupres=[]
    for groups in ngroups:
        test=TestDistances(create_random_data(groups=groups,nodespergroup=10,loci=10,alleletype=200))
        ngroupres.append(test.test_list(iterations,test.tests_fast))
    results["random_fast"]=(ngroups,ngroupres)
    
    iterations=1000
    ngroups=range(1000,3001,500)
    ngroupres=[]
    for groups in ngroups:
        test=TestDistances(create_random_data(groups=groups,nodespergroup=3,loci=10,alleletype=200))
        ngroupres.append(test.test_list(iterations,test.tests_groups))
    results["random_groups"]=(ngroups,ngroupres)

    iterations=10
    ngroups=range(10,301,50)
    ngroupres=[]
    for groups in ngroups:
        test=TestDistances(create_random_data(groups=groups,nodespergroup=10,loci=10,alleletype=200))
        ngroupres.append(test.test_list(iterations,test.tests_slow))
    results["random_slow"]=(ngroups,ngroupres)



    iterations=1000
    test=TestDistances("otus.txt",groups=False,ind_distance="czekanowski",ploidity=1)
    tests=["get_dm_ind",
           "get_globalClustering_ind",
           "auto_threshold_ind",
           "shuffle_alleles",
           "remove_clones",
           "get_mst_ind"]                
    results["tests_otus"]=tests
    results["otus"]=test.test_list(iterations,tests)


    pickle.dump(results,open("results.txt",'w'))


if __name__=="__main__":
    run_tests()

