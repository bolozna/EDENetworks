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
#Do not show warnings to users
from pylab import warnings
warnings.filterwarnings('ignore')

from Tkinter import *
from netpython import netio,pynet,netext,percolator,eden,visuals,netanalysis,models,transforms,communities,dialogues
import tkFileDialog

from eden_analyzer_v2 import *


#->This is a modification to the network data structure to allow 
#edge weights of 0, which would normally indicate a missing edge.
#Basically the old weight set/get functions are replaced such
#that setting weight to w always assignes w+1 as a edge weight
#and getting w always returns w-1. This means that setting edge
#weights to 0 does not delete them as usual. Edges can still be
#deleted with the __del__ method. 
temp_get=pynet.VirtualNet.__getitem__
def new_get(self,args):
    if isinstance(args, tuple):
        v=temp_get(self,args)
        if v!=0.0:
            return v-1
        else:
            return 0.0 #netio uses += operator to add edges
            #return NaN
    else:
        return temp_get(self,args)
pynet.VirtualNet.__getitem__=new_get

temp_set=pynet.VirtualNet.__setitem__
def new_set(self,key,val):
    temp_set(self,key,val+1)
pynet.VirtualNet.__setitem__=new_set

temp_del=pynet.VirtualNet.__delitem__
def new_del(self,args):
    if isinstance(args, tuple):
        self[args[0], args[1]]=-1
    else:
        temp_del(self,args)
pynet.VirtualNet.__delitem__=new_del

def new_node_contains(self,nodeName):
    if nodeName in self.net._nodes:
        myIndex=self.net._nodes[self.name]
        otherIndex=self.net._nodes[nodeName]
        return self.net._getEdge(myIndex,otherIndex)!=0
    else:
        return False
pynet.Node.__contains__=new_node_contains
#<-


#The main that just launces an instance of the MainWindow.
if __name__=='__main__':
    root=Tk()
    root.title("EDENetworks")
    if os.name=="nt":
        root.wm_iconbitmap("logo16.ico",default="logo16.ico")
    root.grab_set()
    MainWindow(root)
    mainloop()
