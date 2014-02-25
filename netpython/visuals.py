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

# -*- coding: latin-1 -*-
import pynet
import netext,percolator,netio,transforms
import os
from pylab import *
import numpy as np
import copy
import random
import shutil
import operator
import tempfile 
import subprocess

# --------------------------------------        

class Myplot(object):
    '''empty container'''
    pass

# ---------------------------------------

def ReturnColorMapPlot(colormap,figsize=(2,0.15)):

    thisfigure=Figure(figsize=figsize,dpi=100)
    axes=thisfigure.add_subplot(1,1,1)
    #myplot.axes.append(myplot.thisFigure.add_axes([0.0,0.0,3.0,1.0]))
    axes.set_axis_off()


    a=np.outer(ones(10),arange(0,1,0.01),)

    axes.imshow(a,aspect='auto',cmap=setColorMap(colormap),origin="lower")

    return thisfigure
    

# ---------------------------------------

def ReturnPlotObject(data,plotcommand='plot',titlestring='',xstring='',
                     ystring='',figsize=(5,4),fontsize=14,fontname='Times',
                     addstr='',labelMultiplier=1.2,plotcolor='r',
                     facecolor="#cacfbe",edgecolor=None):
    """Input: data [[xseries1,yseries1],[xseries2,yseries2],...]
    (float,int, whatever)

    plotcommand = matplotlib's plot command (default 'plot', could
    also be 'loglog' etc), titlestring=plot title, xstring=x-axis
    label, ystring=y-axis label.

    Outputs a container object, corresponding to a matplotlib
    plot. This can be displayed in various ways, eg. on a TkInter
    canvas:
    myplot.canvas=FigureCanvasTkAgg(myplot.thisFigure,master=plotbody)
    myplot.canvas.show() plotbody.pack() where plotbody is a Frame
    object.


    quick hack: addstr takes in arguments for plotting command,
    e.g. ",'ro'","width=0.1", etc. To be fixed.
    """
    Nplots = len(data)

    ystrings=ystring.split(':')
    
    if len(ystrings) < Nplots:
        for i in range(0, Nplots-1):
            ystrings.append(ystrings[0])

    titlestrings = titlestring.split(':')
    if len(titlestrings) < Nplots:
        for i in range(0, Nplots-1):
            titlestrings.append(titlestrings[0])

    subplotstring = str(Nplots)+'1'

    myplot = Myplot()
    myplot.thisFigure = Figure(figsize=figsize,dpi=100,facecolor=facecolor,
                               edgecolor=edgecolor)
    myplot.axes=[]

    axisfactor=1.0/Nplots
    
    for i in range(Nplots):
        leftcorner=0.2
        width=0.6
        bottom=(0.2+float(i))*axisfactor
        top=0.6*axisfactor

        font={'fontname':'Times','fontsize':fontsize+2}
        
        myplot.axes.append(myplot.thisFigure.add_axes([leftcorner,bottom,
                                                       width,top],
                                                      title=titlestrings[i],
                                                      xlabel=xstring,
                                                      ylabel=ystrings[i]))
        s = "myplot.axes[%d].%s(data[%d][0],data[%d][1]" % (i,plotcommand,i,i)
        s += addstr + ')'
        eval(s)

        myplot.axes[i].set_title(titlestrings[i],**font)
        
        for tick in myplot.axes[i].get_xticklabels():
            tick.set_fontsize(fontsize)
            tick.set_fontname(fontname)
                                    
        for tick in myplot.axes[i].get_yticklabels():
            tick.set_fontsize(fontsize)
            tick.set_fontname(fontname)

        labels = [myplot.axes[i].get_xaxis().get_label(),
                  myplot.axes[i].get_yaxis().get_label()]
        [lb.set_size(labelMultiplier*fontsize) for lb in labels]

    myplot.thisFigure.subplots_adjust(hspace=0.5)

    return myplot

def normalizeValue(value,valueLimits):
    """Transforms a numerical value to the range (0,1).

    It is intended that the user should set valueLimits such that the
    true values fall between the limits. If this is not the case,
    values above given maxval or below given minval are truncated. The
    rest of the values are transformed linearly, such that the range
    (given minval, given maxval) becomes (0,1).
     
    normalizedValue= (true_val-given minval)/(given_maxval-given_minval)
    """ 
    if (valueLimits[0]-valueLimits[1]) == 0: 
        # If given minval and maxval are the same, all values will be
        # equal.
        normalizedValue=1
    elif value < valueLimits[0]:
        # If value is smaller than given minval
        normalizedValue=0
    elif value>valueLimits[1]:
        # If value is larger than given maxval
        normalizedValue=1
    else:
        normalizedValue=(value-valueLimits[0])/float(valueLimits[1] -
                                                     valueLimits[0])
    return normalizedValue 


def setColorMap(colorMap):
    """Set a colormap for edges.

    Two options of our own ('orange' and 'primary') are available in
    addition to the 150 pylab readymade colormaps (which can be listed
    with help matplotlib.cm ).

    Usage:
        myMap = setColorMap('bone')
    """
    if hasattr(colorMap, '_segmentdata'):
        return colorMap

    known_colormaps = ('primary', 'orange', 'bluered')
    if colorMap in known_colormaps:
        if colorMap == 'primary':
            # Jari's map: yellow->blue->red 
            segmentdata={'red': ( (0,1,1),(0.5,0,0), (1,1,1)  ),
                         'green': ( (0,1,1), (0.5,0.5,0.5), (1,0,0) ),
                         'blue': ( (0,0,0), (0.5,1,1), (1,0,0) )}
        elif colorMap=='orange':
            # Riitta's color map from white through yellow and orange to red 
            segmentdata = { 'red'  : ( (0.,.99,.99), 
                                       (0.2,.98,.98), 
                                       (0.4,.99,.99), 
                                       (0.6,.99,.99), 
                                       (0.8,.99,.99), 
                                       (1.0,.92,.92) ),
                            'green': ( (0,0.99,0.99), 
                                       (0.2,.89,.89),  
                                       (0.4,.80,.80), 
                                       (0.6,.50,.50), 
                                       (0.8,.33,.33), 
                                       (1.0,.10,.10) ),
                            'blue' : ( (0,.99,.99), 
                                       (0.2,.59,.59), 
                                       (0.4,.20,.20), 
                                       (0.6,0.0,0.0), 
                                       (0.8,0.0,0.0), 
                                       (1.0,.03,.03) )  }
        elif colorMap=='bluered':
            segmentdata={'red':  ( (0,0,0), 
                                   (0.17,0.25,0.25), 
                                   (0.33,0.7,0.7), 
                                   (0.5,.87,.87), 
                                   (0.67,.97,.97),  
                                   (0.83,.93,.93), 
                                   (1,.85,.85) ),
                         'green': ( (0,0,0), 
                                    (0.1667,0.53,0.53), 
                                    (0.3333,.8,.8), 
                                    (0.5,.9,.9), 
                                    (0.6667,.7,.7),
                                    (0.8333,.32,.32), 
                                    (1,.07,.07) ),
                         'blue': ( (0,.6,.6),  
                                   (0.1667,.8,.8),    
                                   (0.3333,1,1),    
                                   (0.5,.8,.8),    
                                   (0.6667,.33,.33),    
                                   (0.8333,.12,.12),
                                   (1,.05,.05) ) }
        myMap = matplotlib.colors.LinearSegmentedColormap(colorMap, segmentdata)
    else:
        try:
            myMap=get_cmap(colorMap)
        except AssertionError:
            comment = "Could not recognize given colorMap name '%s'" % colorMap
            raise AssertionError(comment)
    return myMap


def getConstantColorMap(rgb=(0,0,0)):
    """Return a colormap with constant color.

    Parameters
    ----------
    rgb : tuple (r, g, b)
        The color as RGB tuple. Each value must be between 0 and 1.

    Return
    ------
    cm : colorMap
        The colormap that has just one constant color.
    """
    cm={
        'red':  ( (0,rgb[0],rgb[0]), 
                  (1,rgb[0],rgb[0]) ),
        'green': ( (0,rgb[1],rgb[1]), 
                   (1,rgb[1],rgb[1]) ),
        'blue': ( (0,rgb[2],rgb[2]), 
                  (1,rgb[2],rgb[2]) ) }

    return matplotlib.colors.LinearSegmentedColormap("constant colormap",cm)  


# ---------------------------------------

def isListOfColors(theList):
    """
    Returns True if each element in the given list can be converted to a color
    by Matplotlib. Otherwise returns False.
    """
    cc=matplotlib.colors.ColorConverter()
    for element in theList:
        try:
            cc.to_rgb(element)
        except ValueError:
            return False
    return True
        

def getNodeColors(net,colorwith="strength",useColorMap="orange",parentnet=[]):
    """Returns a dictionary {node:color}. The colors are set based
    on either node strengh (colorwith="strength", default) 
    or any nodeProperty. For cases where e.g. nodes which have been thresholded
    out (k=0), the input parameter parentnet can be used - parentnet should contain the original
    network *before* thresholding, i.e. containing all original nodes and
    their attributes. IF parentnet is given, i) if strength is used, its nodes
    which are NOT in net colored gray, ii) if properties
    are used, its nodes are colored similarly to those nodes in net. Also the
    dictionary which is returned contains then all nodes in parentnet"""

    myNodeColors=setColorMap(useColorMap)

    nodeColors={}

    if colorwith=="strength":

        if hasattr(net,'matrixtype'):
            if net.matrixtype==0:        
                net=transforms.dist_to_weights(net)

        strengths = netext.strengths(net)
        max_value = max(strengths.values())
        min_value = min(strengths.values())

        if len(parentnet)>0:        # if we want the dict to include nodes not in net
            for node in parentnet:
                if node in net:     # if the node is in net, use its strength for color
                    nodeColors[node]=setColor(strengths[node],(min_value,max_value),myNodeColors)
                else:               # otherwise color it gray
                    nodeColors[node]=(0.5,0.5,0.5)
        else:
            for node in net:        # if parentnet not given, just color nodes by strength
                nodeColors[node]=setColor(strengths[node],(min_value,max_value),myNodeColors)
    else:

        numeric_props=netext.getNumericProperties(net)
        # first check if colorwith is a numeric property
        if colorwith in numeric_props:
            values=[]
            if len(parentnet)>0:    # again if we want to include nodes not in net
                for node in parentnet:  # first get min and max value of property
                    values.append(parentnet.nodeProperty[colorwith][node])

                min_value=min(values)
                max_value=max(values)
                for node in parentnet: # then set colors according to property
                    nodeColors[node]=setColor(parentnet.nodeProperty[colorwith][node],(min_value,max_value),myNodeColors)
            else:                   # otherwise do the above for nodes in net
                for node in net:
                    values.append(net.nodeProperty[colorwith][node])
                
                min_value=min(values)
                max_value=max(values)

                for node in net:
                    nodeColors[node]=setColor(net.nodeProperty[colorwith][node],(min_value,max_value),myNodeColors)
        
        else:
            # colorwith is not a numeric property, so look up unique values
            # and give them integer numbers

            values={} # empty dict for values
           
            if len(parentnet)>0:# if there are nodes not in net
                props=list(set(parentnet.nodeProperty[colorwith].values()))
            else:
                props=list(set(net.nodeProperty[colorwith].values()))

            #Check if properties can be converted to colors:
            if isListOfColors(props):
                propToColor={}
                cc=matplotlib.colors.ColorConverter()
                for p in props:
                    propToColor[p]=cc.to_rgb(p)
                if len(parentnet)>0:
                    for node in parentnet:
                        nodeColors[node]=propToColor[parentnet.nodeProperty[colorwith][node]]
                else:
                    for node in net:
                        nodeColors[node]=propToColor[parentnet.nodeProperty[colorwith][node]]                    
            else:
                for i,prop in enumerate(props):
                    values[prop]=i+1
                # now all property strings have a numerical value
                min_value=1
                max_value=max(values.values())
                if len(parentnet)>0:

                    for node in parentnet:
                        nodeColors[node]=setColor(values[parentnet.nodeProperty[colorwith][node]],(min_value,max_value),myNodeColors)
                else:
                    for node in net:
                        nodeColors[node]=setColor(values[net.nodeProperty[colorwith][node]],(min_value,max_value),myNodeColors)



    if len(nodeColors)==0:  # finally if for whatever reason no nodes were colored, just set them gray
        if len(parentnet)>0:
            for node in parentnet:
                nodeColors[node]=(0.5,0.5,0.5)
        else:
            for node in net:
                nodeColors[node]=(0.5, 0.5, 0.5)

    return nodeColors

# ------------------------------------------

def getNodeSizes(net,size_by="strength",minsize=2.0,maxsize=6.0):
    """Returns a dictionary {node:size} for visualizations. The sizes
    are set using either node strength"""

    nodeSizes={}

    if size_by=="strength":

        if hasattr(net,'matrixtype'):
            if net.matrixtype==0:        
                net=transforms.dist_to_weights(net)

        strengths = netext.strengths(net)
        maxs = max(strengths.values())
        mins = min(strengths.values())           

        if maxs==mins:
            A=0
        else:
            A=(maxsize-minsize)/(maxs-mins)
        B=maxsize-A*maxs

        for node in net:
            nodeSizes[node]=A*strengths[node]+B

    elif size_by=="fixed":
        for node in net:
            nodeSizes[node]=maxsize
    else:
        numeric_props=netext.getNumericProperties(net)
        if size_by in numeric_props:
            values=[]
            for node in net:
                values.append(net.nodeProperty[size_by][node])

            minval=min(values)
            maxval=max(values)

            if maxval==minval:
                A=0
            else:
                A=(maxsize-minsize)/(maxval-minval)

            B=maxsize-A*maxval
            for node in net:
                nodeSizes[node]=A*net.nodeProperty[size_by][node]+B

    return nodeSizes
          

def setColor(value,valueLimits,colorMap):
    """Converts a numerical value to a color.

    The value is scaled linearly to the range (0...1) using the
    function normalizeValue and the limits valueLimits. This scaled
    value is used to pick a color from the given colormap. The
    colormap should take in values in the range (0...1) and produce a
    three-tuple containing an RGB color, as in (r,g,b).
    """
    if valueLimits[0] < valueLimits[1]: 
        normalizedValue = normalizeValue(value,valueLimits)
        color = colorMap(normalizedValue)
    else:
        color=(0.5,0.5,0.5)  # gray if all values are equal
    return color


def setEdgeWidth(value,weightLimits,minwidth,maxwidth):
    """Transforms edge weights to widths in the range (minwidth,
    maxwidth). If given minwidth and maxwidth are the same, simply use
    that given width.
    """
    if not(weightLimits[0]-weightLimits[1])==0:
        # Normalizes the weight linearly to the range (0,1)
        normalizedWeight=normalizeValue(value,weightLimits)  
        # Transforms the normalized weight linearly to the range
        # (minwidth,maxwidth)
        width=minwidth+normalizedWeight*(maxwidth-minwidth)   
    else:
        # If given minwidth and maxwidth are the same, simply use that width.
        width=minwidth 
    return width


def visualizeNet(net, coords=None, axes=None, frame=False,
                 scaling=True, margin=0.025,
                 nodeShapes=None, defaultNodeShape='o',
                 nodeSizes=None, defaultNodeSize=None,
                 nodeColors=None, defaultNodeColor=None,
                 nodeEdgeColors=None, defaultNodeEdgeColor='black',
                 nodeLabels=None, nodeLabelSize=None, labelAllNodes=False,
                 labelPositions=None, defaultLabelPosition='out',
                 edgeColors=None, defaultEdgeColor=None,
                 edgeWidths=None, defaultEdgeWidth=None,
                 nodeEdgeWidths=None, defaultNodeEdgeWidth=0.2,
                 edgeLabels=None, edgeLabelSize=None, labelAllEdges=False,
                 nodePlotOrders=None, defaultNodePlotOrder=1,
                 edgePlotOrders=None, defaultEdgePlotOrder=0):
    """Visualize a network.
    
    Note that all sizes (node size, link width, etc.) are given in
    points: one point equals 1/72th of an inch.

    Basic parameters
    ----------------
    net : pynet.SymmNet
        The network to visualize
    coords : dictionary of tuples {node_ID: (x,y)}
        Coordinates of all nodes. If None, the coordinates will be
        calculated by Himmeli. The x and y coordinates are assumed to
        have the same scale, so for example an edge between nodes 0
        and 1 with `coords[0]=(0,0)` and `coords[1]=(2,2)` will be at
        an angle of 45 degrees.
    axes : pylab.axes object
        If given, the network will be drawn in this axis. Otherwise a
        new figure is created for the plot and the figure handle is
        then returned.
    frame : bool
        If False, the frame will be not be shown in the plot. You can
        still use axis labels.
    scaling : bool
        If True, the coordinate axes will be scaled for best fit. If
        false, the coordinate axes will not be altered.
    margin : float (>= 0)
        The relative size of empty margin around the network. Margin
        of 0.0 means that some nodes touch the edge of the plot,
        margin of 0.2 adds 20 % on all sides etc. This parameter has
        an effect only if scaling is True.

    Defining node and edge colors
    -----------------------------

    Colors for nodes and edges are defined similarly. The following
    explains the procedure for nodes; to control edge coloring simply
    replace the word 'node' (or 'Node') with the word 'edge' (or
    'Edge') in the parameter names.

    The color of a node is defined with the dictionary `nodeColors`:
    key is the node index and the value is any valid coloring scheme
    (see below under 'Coloring schemes'). If a node index is not in
    `nodeColors`, it is colored according to `defaultNodeColor`. This
    variable can have the same values as the values in `nodeColors`.

    Coloring schemes
    ----------------

    A constant color can be defined in any way allowed by pylab. For
    example 'k', 'black' and (0,0,0) all give black color.

    Alternatively the color can be based on the node strength, degree
    or any node property. In this case the coloring definition is a
    dictionary. The following examples illustrate the idea:
    
    color_scheme = {'by':'weight', 'scale':'log', 'cmap':'winter'}
    color_scheme = {'by':'degree', 'scale':'lin', 'min':1, 'max':10}
    color_scheme = {'by':'property:myProperty', 'scale':'log'}

    The possible keys and their default values are
    
        KEY      DEFAULT VALUE       OTHER POSSIBLE VALUES
        'by'     'strength'/'weight' 'degree', 'property:<property_name>'
        'scale'  'log'               'lin'
        'cmap'   'jet'               Any colormap
        'min'    (Min value in data) Any integer x,     1 <= x <= 'max'
        'max'    (Max value in data) Any integer x, 'min' <= x         

    Any keys that are omitted are filled in with the default
    value. Note the syntax for using node properties, where the word
    'property' is followed by a semicolon and the property name.

    Node size
    ---------

    The node size is controlled with a syntax similar to that used
    with colors. Node size is defined by dictionary `nodeSizes`, and
    if a node is not in it, the default value given by
    `defaultNodeSize` is used. 

    The value can be a single integer, which gives the node size in
    pixels. Alternatively the node size can be controlled by node
    strength, degree or any property:
    
    node_size_scheme = {'by':'strength', 'scale':'log', 'min':2, 'max':10}
    node_size_scheme = {'by':'degree', 'scale':'lin'}
    node_size_scheme = {'by':'property:myProperty', 'scale':'log'}

    The possible keys and their default values are
        KEY        DEFAULT VALUE       OTHER POSSIBLE VALUES
        'by'       'strength'          'degree', 'property:<property_name>'
        'scale'    'log'               'lin'
        'min'      (Min value in data) int;          1 <= x <= 'max'
        'max'      (Max value in data) int;      'min' <= x         
        'min_size'   1                 int;          1 <= x <= 'max_size'
        'max_size'   6                 int; 'min_size' <= x         
    Again, keys that are omitted are filled with default values.

    Edge width
    ----------

    Edge width is defined by dictionary `edgeWidths`, and if an edge
    is not in it, the default value given by `defaultEdgeWidth` is
    used.

    The value can be a single integer, which gives the edge width in
    pixels. Alternatively the edge width can be controlled by edge
    weight:
    
    edge_width_scheme = {'by':'weight', 'scale':'log', 'min':1, 'max':5}

    The possible keys and their default values are
    
        KEY         DEFAULT VALUE       OTHER POSSIBLE VALUES
        'by'        'weight'          
        'scale'     'log'               'lin'
        'min'       (Min value in data) int;   1 <= x <= 'max'
        'max'       (Max value in data) int;   'min' <= x         
        'min_size'    0.2               float; 1 <= x <= 'max_width'
        'max_size'    2.0               float; 'min_width' <= x         

    Note that the 'by'-key can always be omitted since it has only one
    possible value.

    Node labels
    -----------

    Node labels can be given in `nodeLabels` dictionary, where the key
    is node index and the value is the corresponding labels. If
    `labelAllNodes` is True, also nodes not in `nodeLabels` will
    receive a label, which is the node index.

    The values in `nodeLabels` are converted to string with str().

    There are two possibilities for the positioning of node labels:
    'in' and 'out' (default). 'in' means that the label is printed at
    the exact position of the node; if the node is hollow, this
    effectively prints the node labes inside the nodes (make sure the
    node size and font sizes are compatible). 'out' prints the label
    next to the node.

    Edge labels
    -----------

    Also edges can have labels, given in `edgeLabels` dictionary,
    where the key is a tuple (i,j) of end node indices. The edge
    labels are always printed on the right side of each edge
    (with direction defined from i to j).

    If `labelAllEdges` is True, also the edges not listed in
    `edgeLabels` will be given a label. In this case the label is
    '(i,j)', where i and j are the indices of the end nodes.

    Return
    ------
    fig : pylab.Figure (None if `axes` is given.)
        Figure object with one axes containing the plotted network
        figure.

    Examples
    --------
    >>> # Construct an example network.
    >>> from netpython import pynet, visuals
    >>> net = pynet.SymmNet()
    >>> net[0][1] = 1.0
    >>> net[1][2] = 3.5
    >>> net[0][2] = 5.0

    >>> # Simplest case: get coordinates from Himmeli, plot into a
    >>> # new figure and save it to disk
    >>> fig = visuals.visualizeNet(net)
    >>> fig.savefig('myNet.eps')

    >>> # Draw the figure in the upper left subfigure, with predefined
    >>> # coordinates. Note that drawNet does not return anything.
    >>> import pylab
    >>> coords = {0:(0,0), 1:(4,0), 2:(2,3)}
    >>> fig = pylab.figure()
    >>> ax = fig.add_subplot(2,2,1)
    >>> visuals.visualizeNet(net, coords=coords, axes=ax)
    """

    #
    # DEFAULT VALUES. These will be used whenever the user has not
    # defined a given value for defaultNodeColor etc.
    # 
    internal_defaultNodeColor = {'by':'strength', 'scale':'log', 'cmap':'jet'}
    internal_defaultEdgeColor = {'by':'weight', 'scale':'log', 'cmap':'jet'}

    internal_defaultNodeSize = {'by':'strength', 'scale':'log',
                                'min_size':2, 'max_size':6}
    internal_defaultEdgeWidth = {'by':'weight', 'scale':'log',
                                 'min_size':0.2, 'max_size':2.0}

    node_label_font_color = 'k'
    node_label_font_size = (8 if nodeLabelSize==None else nodeLabelSize)
    edge_label_font_color = 'k'
    edge_label_font_size = (5 if edgeLabelSize==None else edgeLabelSize)

    #
    # PROCESS INPUT PARAMETERS
    #
    
    if coords is None:
        coords = Himmeli(net).getCoordinates()

    fig = None
    if axes is None:
        fig = figure()
        axes = fig.add_subplot(111)

    nodeShapes = (nodeShapes or {})

    nodeColors = (nodeColors or {})
    defaultNodeColor = (defaultNodeColor or {})
    if isinstance(defaultNodeColor, dict):
        for k,v in internal_defaultNodeColor.iteritems():
            if k not in defaultNodeColor:
                defaultNodeColor[k] = v

    nodeEdgeColors = (nodeEdgeColors or {})
                
    edgeColors = (edgeColors or {})
    defaultEdgeColor = (defaultEdgeColor or {})
    if isinstance(defaultEdgeColor, dict):
        for k,v in internal_defaultEdgeColor.iteritems():
            if k not in defaultEdgeColor:
                defaultEdgeColor[k] = v
    
    nodeSizes = (nodeSizes or {})
    if defaultNodeSize is None:
        defaultNodeSize = {}
    if isinstance(defaultNodeSize, dict):
        for k,v in internal_defaultNodeSize.iteritems():
            if k not in defaultNodeSize:
                defaultNodeSize[k] = v

    edgeWidths = (edgeWidths or {})
    if defaultEdgeWidth is None:
        defaultEdgeWidth = {}
    if isinstance(defaultEdgeWidth, dict):
        for k,v in internal_defaultEdgeWidth.iteritems():
            if k not in defaultEdgeWidth:
                defaultEdgeWidth[k] = v

    nodeEdgeWidths = (nodeEdgeWidths or {})

    nodeLabels = (nodeLabels or {})
    labelPositions = (labelPositions or {})

    edgeLabels = (edgeLabels or {})

    nodePlotOrders = (nodePlotOrders or {})
    edgePlotOrders = (edgePlotOrders or {})

    if margin < 0: margin = 0.0

    #
    # AUXILIARY FUNCTIONS
    #

    def points_to_data(ax, x_points):
        """Converts point length `x_points` (1/72th of inch) to length
        in data coordinates. Note that this only works (in general) if
        the axis are equal; otherwise converting x and y-coordinates
        gives different results."""
        ax_pos = ax.get_position()
        V = ax.axis()
        width_inches = ax_pos.width*ax.get_figure().get_figwidth()
        return (V[1]-V[0])*x_points/(72*width_inches)

    def scaled(scaling_type, value, value_limits, final_limits):

        def lin_scaling(value, value_limits, final_limits):
            value_span = value_limits[1] - value_limits[0]
            final_span = final_limits[1] - final_limits[0]
            if final_span == 0:
                return final_limits[0]
            if value_span == 0:
                p = 0.5
            else:
                p = float(value - value_limits[0])/value_span
            return final_limits[0]+p*final_span

        if value <= value_limits[0]:
            return final_limits[0]
        if value >= value_limits[1]:
            return final_limits[1]

        if scaling_type == 'log' or scaling_type == 'logarithmic':
            return lin_scaling(np.log(value),
                               np.log(value_limits),
                               final_limits)
        else:
            return lin_scaling(value, value_limits, final_limits)

    def determine_size(scheme, i, net, values, limits, defaults):
        if not isinstance(scheme, dict):
            return scheme
        else:
            # Determine what defines the size. Calculate the limits
            # for this property if not yet done.
            size_by = scheme.get('by', defaults['by'])
            if size_by not in limits:
                property_name = "".join(size_by.split(':')[1:])
                np_ = sorted(net.nodeProperty[property_name].values())
                limits[size_by] = (np_[0], np_[-1])
            if size_by not in values:
                property_name = "".join(size_by.split(':')[1:])
                values[size_by] = net.nodeProperty[property_name][i]
                
            scale = scheme.get('scale', defaults['scale'])
            val_min = scheme.get('min', limits[size_by][0])
            val_max = scheme.get('max', limits[size_by][1])
            size_min = scheme.get('min_size', defaults['min_size'])
            size_max = scheme.get('max_size', defaults['max_size'])

            #print size_by, scale, val_min, val_max, size_min, size_max
            return scaled(scale, values[size_by], [val_min, val_max],
                          [size_min, size_max])
                    
    def determine_color(scheme, i, net, values, limits, defaults):
        if not isinstance(scheme, dict):
            return scheme
        else:
            color_by = scheme.get('by', defaults['by'])
            if color_by not in limits:
                property_name = "".join(color_by.split(':')[1:])
                np_ = sorted(net.nodeProperty[property_name].values())
                limits[color_by] = (np_[0], np_[-1])
            if color_by not in values:
                property_name = "".join(color_by.split(':')[1:])
                values[color_by] = net.nodeProperty[property_name][i]
                
            scale = scheme.get('scale', defaults['scale'])
            cmap = scheme.get('cmap', defaults['cmap'])
            val_min = scheme.get('min', limits[color_by][0])
            val_max = scheme.get('max', limits[color_by][1])

            cval = scaled(scale, values[color_by],
                          [val_min, val_max], [0.0,1.0])
            cm = setColorMap(cmap)
            return cm(float(cval))

    def get_edge_angle(xcoords, ycoords):
        """Return the edge angle in [-pi/2, 3*pi/2]."""
        dx, dy = xcoords[1]-xcoords[0], ycoords[1]-ycoords[0]
        if dx == 0:
            if dy > 0: theta = np.pi/2
            else: theta = -np.pi/2
        elif dx > 0:
            theta = np.arctan(dy/dx)
        else:
            theta = np.arctan(dy/dx) + np.pi
        return theta

    def edge_label_pos(axes, xcoords, ycoords, edge_width, label_size, offset=1.5):
        """Return the baseline position and label rotation (in angles)
        for an edge label. The label will be on the right side of the
        edge, in proper orientation for reading, and located `offset`
        points from the edge.
        """
        theta = get_edge_angle(xcoords, ycoords)
        if theta > -np.pi/2 and theta < np.pi/2:
            # Edge goes from left to right.
            label_rotation = theta
        else:
            # Edge goes from right to left.
            label_rotation = theta - np.pi

        offset_points = offset+0.5*edge_width+0.5*label_size
        label_position = (0.5*sum(xcoords), 0.5*sum(ycoords))
        offset_dir = theta - np.pi/2
        label_offset = (offset_points*np.cos(offset_dir),
                        offset_points*np.sin(offset_dir))
        return label_position, label_offset, 180*label_rotation/np.pi

    def draw_edge(axes, xcoords, ycoords, width, color, symmetric, zorder,
                  nodesize):
        if width == 0: return
        if symmetric:
            axes.plot(xcoords, ycoords, '-', lw=width,
                      color=color, zorder=zorder)
        else:
            dx, dy = xcoords[1]-xcoords[0], ycoords[1]-ycoords[0]
            theta = get_edge_angle(xcoords, ycoords)
                
            x_diff = points_to_data(axes, 0.5*nodesize*np.cos(theta))
            y_diff = points_to_data(axes, 0.5*nodesize*np.sin(theta))
            x, y = xcoords[0]+x_diff, ycoords[0]+y_diff
            dx, dy = dx-2*x_diff, dy-2*y_diff

            #print ("Start (%.4f, %.4f), End (%.4f, %.4f)" 
            #       % (x, y, xcoords[0]+dx, ycoords[0]+dy))
            
            arrow_width = points_to_data(axes, width)

            # See pylab.matplotlib.patches.FancyArrow for
            # documentation of the arrow command.
            axes.arrow(x, y, dx, dy,
                       color=color, width=arrow_width,
                       head_width=12*arrow_width,
                       head_length=18*arrow_width,
                       shape='left', overhang=0.2,
                       length_includes_head=True, 
                       head_starts_at_zero=False,
                       zorder=zorder)

    def draw_node(axes, x, y, shape, color, size, edgecolor, edgewidth, zorder):
        if size == 0: return
        axes.plot([x], [y], shape,
                  markerfacecolor=color,
                  markeredgecolor=edgecolor,
                  markeredgewidth=edgewidth,
                  markersize=size,
                  zorder=zorder)

    #
    # INITIALIZE SOME DATA STRUCTURES
    #
        
    # Find out the minimum and maximum value for strength and degree.
    strengths = netext.strengths(net)
    smin, smax = min(strengths.values()), max(strengths.values())
    degrees = netext.deg(net)
    dmin, dmax = min(degrees.values()), max(degrees.values())

    limits = {"strength":(smin, smax), "degree":(dmin,dmax)}

    #
    # INITIALIZE AXES SIZE
    #

    node_diameters = {};
    for nodeIndex in net:
        values = {"strength": strengths[nodeIndex],
                  "degree": degrees[nodeIndex]}
        node_diameters[nodeIndex] = (determine_size(nodeSizes.get(nodeIndex, defaultNodeSize),
                                                    nodeIndex, net, values, limits,
                                                    defaultNodeSize) +
                                     determine_size(nodeEdgeWidths.get(nodeIndex, 
                                                                       defaultNodeEdgeWidth),
                                                    nodeIndex, net, values, limits,
                                                    defaultNodeEdgeWidth))

    if scaling:
        # Make axis equal making sure the nodes on the edges are not
        # clipped. We cannot use `axis('equal')` because nothing has been
        # drawn yet, but we still need to do this so the arrows will be
        # draw properly in the plotting phase. This is a bit tricky
        # because the axes might not be square.

        max_node_diameter = max(node_diameters.values())
        y_coords = sorted(map(operator.itemgetter(1), coords.values()))
        x_coords = sorted(map(operator.itemgetter(0), coords.values()))
        ax_pos = axes.get_position()
        x_span = x_coords[-1] - x_coords[0]
        y_span = y_coords[-1] - y_coords[0]
        ax_width_inches = ax_pos.width*axes.get_figure().get_figwidth()
        ax_height_inches = ax_pos.height*axes.get_figure().get_figheight()

        if (x_span*ax_height_inches > y_span*ax_width_inches):
            # The x-span dictates the coordinates. Calculate the margin
            # necessary to fit in the nodes on the edges.
            rad_frac = 0.5*max_node_diameter/(72*ax_width_inches)
            x_margin = x_span*rad_frac/(1-2*rad_frac)
            x_min, x_max = x_coords[0]-x_margin, x_coords[-1]+x_margin
            y_mid = 0.5*(y_coords[-1] + y_coords[0])
            y_axis_span = (x_max-x_min)*ax_height_inches/ax_width_inches
            y_min, y_max = y_mid - 0.5*y_axis_span, y_mid + 0.5*y_axis_span,
        else:
            # The y-span dictates the coordinates. Calculate the margin
            # necessary to fit in the nodes on the edges.
            rad_frac = 0.5*max_node_diameter/(72*ax_height_inches)
            y_margin = y_span*rad_frac/(1-2*rad_frac)
            y_min, y_max = y_coords[0]-y_margin, y_coords[-1]+y_margin
            x_mid = 0.5*(x_coords[-1] + x_coords[0])
            x_axis_span = (y_max-y_min)*ax_width_inches/ax_height_inches
            x_min, x_max = x_mid - 0.5*x_axis_span, x_mid + 0.5*x_axis_span,

        y_span, x_span = y_max - y_min, x_max - x_min
        axes.set_ylim(ymin=y_min-margin*y_span, ymax=y_max+margin*y_span)
        axes.set_xlim(xmin=x_min-margin*x_span, xmax=x_max+margin*x_span)
    prev_autoscale = axes.get_autoscale_on()
    axes.set_autoscale_on(False)

    #
    # DRAW EDGES
    #

    edges = list(net.edges)
    if edges:
        # Sort by edge weight.
        edges.sort(key=operator.itemgetter(2))
        limits['weight'] = (edges[0][2], edges[-1][2])

        # DEBUGGING: Print statistics of edge weights.
        #import data_utils
        #weights_ = map(operator.itemgetter(2), edges) 
        #print "Edges:", limits['weight'], " 10/25/50/75/90:",
        #mp_ = len(edges)/2
        #print "%d, %d, %d, %d, %d" % (int(data_utils.percentile(weights_, 0.1)),
        #                              int(data_utils.percentile(weights_, 0.25)),
        #                              int(data_utils.percentile(weights_, 0.5)),
        #                              int(data_utils.percentile(weights_, 0.75)),
        #                              int(data_utils.percentile(weights_, 0.9)))

        for i,j,w in edges:
            values = {"weight": w,
                      "strength": strengths[j],
                      "degree": degrees[j]}
            
            # Determine edge width.
            if (i,j) in edgeWidths:
                width = determine_size(edgeWidths[(i,j)], (i,j), net,
                                       values, limits, defaultEdgeWidth)
            elif (j,i) in edgeWidths:
                width = determine_size(edgeWidths[(j,i)], (j,i), net,
                                       values, limits, defaultEdgeWidth)
            else:
                width = determine_size(defaultEdgeWidth, (i,j), net,
                                       values, limits, defaultEdgeWidth)

            # Determine edge color.
            if (i,j) in edgeColors:
                color = determine_color(edgeColors[(i,j)], (i,j), net,
                                       values, limits, defaultEdgeColor)
            elif (j,i) in edgeColors:
                color = determine_color(edgeColors[(j,i)], (j,i), net,
                                       values, limits, defaultEdgeColor)
            else:
                color = determine_color(defaultEdgeColor, (j,i), net,
                                       values, limits, defaultEdgeColor)

            if (i,j) in edgePlotOrders:
                zorder = edgePlotOrders[(i,j)]
            elif (j,i) in edgePlotOrders:
                zorder = edgePlotOrders[(j,i)]
            else:
                zorder = defaultEdgePlotOrder

            # FOR DEBUGGING:
            #print "Edge (%d,%d) : %.1f %s %f" % (i,j,width,str(color),zorder)
            xcoords, ycoords = (coords[i][0], coords[j][0]), (coords[i][1], coords[j][1])
            draw_edge(axes, xcoords, ycoords, width, color,
                      net.isSymmetric(), zorder, node_diameters[j])

            # Add edge label.
            if (labelAllEdges or (i,j) in edgeLabels or 
                (net.isSymmetric() and (j,i) in edgeLabels)):
                if (i,j) in edgeLabels:
                    label = str(edgeLabels[(i,j)])
                elif (net.isSymmetric() and (j,i) in edgeLabels):
                    label = str(edgeLabels[(j,i)])
                else:
                    label = "(%d,%d)" % (i,j)

                lpos, loffset, lrot = edge_label_pos(axes, xcoords, ycoords,
                                                     width, edge_label_font_size)
                axes.annotate(label, lpos, xytext=loffset,
                              textcoords='offset points',
                              color=edge_label_font_color,
                              size=edge_label_font_size,
                              horizontalalignment='center',
                              verticalalignment='center',
                              rotation=lrot,
                              zorder=zorder+0.5)

    #
    # DRAW NODES
    #

    max_node_diameter = 0;

    for nodeIndex in net:
        values = {"strength": strengths[nodeIndex],
                  "degree": degrees[nodeIndex]}

        # Determine node shape.
        shape = nodeShapes.get(nodeIndex, defaultNodeShape)

        # Determine node size (and update max size).
        size = determine_size(nodeSizes.get(nodeIndex, defaultNodeSize),
                              nodeIndex, net, values, limits,
                              defaultNodeSize)

        # Determine node edge width.
        edgewidth = determine_size(nodeEdgeWidths.get(nodeIndex, defaultNodeEdgeWidth),
                                   nodeIndex, net, values, limits,
                                   defaultNodeEdgeWidth)
        max_node_diameter = max(max_node_diameter, size+edgewidth)

        # Determine node color
        color = determine_color(nodeColors.get(nodeIndex, defaultNodeColor),
                                nodeIndex, net, values, limits,
                                defaultNodeColor)

        # Determine node edge color
        edgecolor = determine_color(nodeEdgeColors.get(nodeIndex, defaultNodeEdgeColor),
                                    nodeIndex, net, values, limits,
                                    defaultNodeEdgeColor)

        # Determine z-order.
        zorder = nodePlotOrders.get(nodeIndex, defaultNodePlotOrder)
    
        # FOR DEBUGGING:
        #print "Node %d : %f %s %f %s" % (nodeIndex, size, str(color), edgewidth, str(edgecolor))

        draw_node(axes, coords[nodeIndex][0], coords[nodeIndex][1],
                  shape, color, size, edgecolor, edgewidth, zorder)

        # Add node label.
        if nodeIndex in nodeLabels or labelAllNodes:
            if nodeIndex in nodeLabels:
                label = str(nodeLabels[nodeIndex])
            else:
                label = str(nodeIndex)

            label_pos = labelPositions.get(nodeIndex, defaultLabelPosition)
            if label_pos == 'out':
                nodeLabel_offset = int(np.ceil(float(size)/2))+1
                axes.annotate(label,
                              (coords[nodeIndex][0],coords[nodeIndex][1]),
                              textcoords='offset points',
                              xytext=(nodeLabel_offset, nodeLabel_offset),
                              color=node_label_font_color,
                              size=node_label_font_size,
                              zorder=zorder+0.5)
            elif label_pos == 'in':
                axes.annotate(label,
                              (coords[nodeIndex][0],coords[nodeIndex][1]),
                              color=node_label_font_color,
                              size=max(5, min(size-1, 0.7*size)),
                              horizontalalignment='center',
                              verticalalignment='center',
                              zorder=zorder+0.5)
 

    # Remove frame.
    if not frame:
        #Using 'axes.set_axis_off()' would also turn of the axis
        #labels, which is too much. The following lines are required
        #to turn of the frame and tick labels while keeping axis
        #labels.
        axes.set_frame_on(False)
        axes.set_xticklabels([])
        axes.xaxis.set_ticks_position('none')
        axes.set_yticklabels([])
        axes.yaxis.set_ticks_position('none')

    # Return autoscaling to the original value.
    axes.set_autoscale_on(prev_autoscale)

    # Return figure. Note that if `axes` was given as a input
    # argument, the returned value is None.
    return fig
    

def VisualizeNet(net, xy, figsize=(6,6), coloredNodes=True, equalsize=False,
                 equalshape=True, labels=None, fontsize=7, showAllNodes=True,
                 nodeColor=None, nodeShape='o', nodeEdgeColor='k', nodeSize=1.0,
                 minnode=2.0, maxnode=6.0, nodeColors=None, nodeSizes=None,
                 nodeShapes=None, bgcolor='white', maxwidth=2.0, minwidth=0.2,
                 uselabels='none', edgeColorMap='winter', weightLimits=None,
                 setNodeColorsByProperty=None, nodeColorMap='winter',
                 nodePropertyLimits=None, nodeLabel_xOffset=None, coloredvertices=None,
                 vcolor=None, vsize=None, frame=False, showTicks=False, 
                 axisLimits=None, baseFig=None,interactive=False,keepAspectRatio=False): 
    """Visualizes a network.

    The coloring of the nodes is decided as follows:
      a) If dictionary `nodeColors` is given, it is used.  If it does
         not contain a color for every node, the rest are colored 
           1) according to property `setNodeColorsByProperty`, if it is
              given, or else
           2) by `nodeColor` if it is given, or
           3) white if neither of the above is given.
      b) If dictionary `nodeColors` is not given, but `nodeColor`
         is given, all nodes are colored with `nodeColor`.
      c) If none of `setNodeColorsByProperty`,`nodeColors` and
         `nodeColor` is given, nodes are colored by strength using the
         colormap `nodeColorMap` (by default 'winter').

    Parameters
    ----------
    net : pynet.SymmNet
        The network to visualize
    xy : dictionary of tuples {node_ID: (x,y,z)}
        Coordinates of all nodes. These usually originate from
        visuals.Himmeli, e.g. 
          h = visuals.Himmeli(net, ...)
          xy = h.getCoordinates()
    figsize : (x,y)
        Size of the output figure in inches. dpi is set to 100.
    coloredNodes : bool
        If True, nodes are colored. Otherwise all nodes are white.
    nodeColors : dict
        Dictionary of node colors by node index.
    nodeSizes : dict
        Dictionary of node sizes by node index.
    nodeShapes : dict
        Dictonary of node shapes by node index.
    minnode : minimum node size, if autoscaling used
    maxnode : maximum node size, if autoscaling used
    nodeColor : RGB color tuple
        Default color of a node. Three values between 0 and 1, for
        example (1.0, 0, 0) is red and (0.5, 0.5, 0.5) is middle gray.
    nodeEdgeColor : any valid matplotlib color (default 'k')
        The color of the edges of nodes. Default is black.
    setNodeColorByProperty : sequence of node indices
        If `setNodeColorsByProperty` is specified, any node not
        appearing in the dictionary `nodeColors` will be colored
        according to the given property (using `nodeColorMap` and
        `nodePropertyLimits`). Option `nodeColors` overrides the
        'setNodeColorsByProperty' option.
    nodeColorMap : str
        A colormap used to color the nodes listed in
        `setNodeColorsByProperty`.
    nodePropertyLimits : [min_val max_val]
        If nodes are coloured according to a nodeProperty, these
        are the min and max values of said property.
    equalsize : bool
        If True, all nodes are of size `nodeSize`. If False, node size
        is based on node strength.
    equalshape : bool
        If True, all nodes are of shape 'nodeShape'. If False, node shape is
        based on 'nodeShapes'
    showAllNodes : bool
        If True, displays disconnected components and nodes which have
        no edges left after e.g. thresholding. (quick hack?)
    bgcolor : sequence (r, g, b)
        Background color as RGB tuple. Default is black.
    minwidth : float
        Minimum width of plotted edges.
    maxwidth : float
        Maximum width of plotted edges.
    labels : dict {nodename:labelstring}
        Dictionary of node labels.
    uselabels : either 'some' (default), 'none' or 'all'
        If 'some', the label is shown for the nodes whose label is
        given in `labels`. If 'all', the node index is shown also for
        those nodes that are not listed in `labels`. If 'none', no
        node labels are printed.
    fontsize : int
        Sets font size for labels.
    edgeColorMap : str or cmap
        Allows the user to set color scheme for edges. Edges are
        always colored according to edge weights, which are first
        normalized to the range (0,1) and then transformed to colors
        using edgeColorMap. There are 150 colormaps available in
        pylab; for a full listing, please see help(pylab.cm) (and look
        for DATA). Or try, for example, edgeColorMap='orange' or
        edgeColorMap='primary', two colormaps of our own that are not
        available in pylab. To make all edges have the same color,
        create a constant color map with getConstantColorMap(rgb).
    weightLimits : tuple (minWeight, maxWeight)
        Provides the minimum and maximum value for weights. If not are
        given, (nearly) the true min and max weights in the network
        will be used. The weightLimits are used for setting edge
        colors and width. They enable the user to plot several
        networks (which may have different min and max weights) so
        that a certain color and width always correspond to a certain
        edge weight. Thus, the color and width in the visualization
        can be used to infer edge weight. If the network turns out to
        contain weights above the given maxWeight (below minWeight)
        these will be rounded downwards (upwards) to the given
        limit.
    nodeLabel_xOffset : float (default nodeSize/40)
        Amount for moving node labels the right so that the text does
        not fall on the nodes.
    frame : bool
        If True, draws a box around the figure.
    showTicks : bool
        If True, adds ticks in the frame. Setting `showTicks` to True
        will always set `frame` to True also.
    axisLimits : tuple ((minX, maxX),(minY, maxY))
        Sets tick limits if `showTicks` is True.
    baseFig : FigureCanvasBase
        If None, the network is drawn on an empty figure, otherwise
        baseFig is used as a starting point.
    interactive : bool
        If True, the nodes can be moved around in the figure. Interaction needs to be
        started by calling fig.startInteraction(). It can be stopped by calling
        fig.stopInteraction(). Note that the coordinates in the xy dict are modified
        when user move the nodes around.
    keepAspectRatio : bool
        If axisLimits are not given, they are inferred from the xy-data. If keepAspectRatio 
        is true, the aspect ratio in the aspect ratios is kept same in axis limits as in the
        given figsize argument. Otherwise the axis are streched to the whole figure.
         

    Return
    ------
    fig : pylab.Figure
        The plotted network figure.


    Examples
    --------
    >>> from netpython import pynet, visuals
    >>> net = pynet.SymmNet()
    >>> net[0][1] = 1.0
    >>> net[1][2] = 3.5
    >>> net[0][2] = 5.0

    >>> # Here are the coordinates, a dictionary that contains 2-tuples 
    >>> xy = {0:(0,0), 1:(4,0), 2:(2,3)}
    >>> # With larger network use Himmeli to calculate coordinates.
    >>> h = visuals.Himmeli(net)
    >>> xy = h.getCoordinates()

    >>> f = FigureCanvasBase(visuals.VisualizeNet(net,xy))
    >>> f.print_eps("myPlot_1.eps", dpi=80.0)

    >>> f2 = FigureCanvasBase(visuals.VisualizeNet(other_net,xy,baseFig=f))

    >>> f=FigureCanvasBase(visuals.VisualizeNet(net,xy,edgeColorMap='orange'))
    >>> f.print_eps("myPlot_2.eps", dpi=80.0)

    >>> f=FigureCanvasBase(visuals.VisualizeNet(net,xy,edgeColorMap='orange',
                                                equalsize=True, nodeSize=16))
    >>> f.print_eps("myPlot_3.eps", dpi=80.0)

    (General questions: Is there a neater way to output the figures
    than using FigureCanvasBase? How can I have a look at the figures
    from within python, without saving them to .eps files?)
    """

   
    def plot_edge(plotobject, xcoords, ycoords, width=1.0, colour='k',
                  symmetric=True):
        if symmetric:
            return plotobject.plot(xcoords,ycoords,'-',lw=width,color=colour)[0]
        else:
            arr = Arrow(xcoords[0], ycoords[0], xcoords[1]-xcoords[0], 
                        ycoords[1]-ycoords[0], edgecolor='none',
                        facecolor=colour,linewidth=width)
            return plotobject.add_patch(arr)[0]


    def plot_node(plotobject,x,y,shape='o',color='w',size=8.0,edgecolor='w'):
        return plotobject.plot([x], [y], 'yo', marker=shape,markerfacecolor=color,
                        markeredgecolor=edgecolor,markersize=size)[0]


    if interactive:
        edgeObjectIndices=pynet.SymmNet()
        edgeObjects=[]
        nodeObjects={}
        nodeLabelObjects={}

    # Warn about obsolete input arguments
    if coloredvertices!=None or vcolor!=None or vsize!=None:
        warnings.warn("\n\n The options \n"
                      "\t coloredvertices, vcolor, and vsize \n"
                      "are now obsolete. Please use instead \n"
                      "\t coloredNodes, nodeColor, and nodeSize.\n")
    if coloredvertices != None:
        coloredNodes = coloredvertices
    if vcolor != None and nodeColor == None:
        nodeColor = vcolor 
    if vsize != None:
        nodeSize = vsize

    if nodeColors is None:
        nodeColors = {}
    if nodeSizes is None:
        nodeSizes = {}
    if nodeShapes is None:
        nodeShapes = {}
    if labels is None:
        labels = {}

    # The following is for the EDEN software, where "nets" or nets
    # derived from matrices can have edge distances instead of weights.
    if hasattr(net,'matrixtype'):
        if net.matrixtype == 0:        
            net=transforms.dist_to_weights(net)

    if baseFig==None:
        thisfigure = Figure(figsize=figsize,dpi=100,facecolor=bgcolor)
        axes = thisfigure.add_subplot(111)
    else:
        thisfigure=baseFig.figure
        thisfigure.set_facecolor(bgcolor)
        axes=thisfigure.gca()

    axes.set_axis_bgcolor(bgcolor)
    
    if frame == False and showTicks == False:
        axes.set_axis_off()

    # Set the color for node labels
    fontcolor='w'
    if bgcolor=='white':
        fontcolor='k'
        
    # First draw all edges, if there are any
    edges=list(net.edges)
    if len(edges)>0:
        wlist=[]
        for edge in edges:
            wlist.append(edge[2])

        wmin=min(wlist)
        wmax=max(wlist)

        # If weightLimits were not given, use (almost) the true min
        # and max weights in the network. Note: using a value slightly
        # below wmin, because otherwise when normalizing the weights,
        # the minimum weights would be transformed to zero and the
        # edges not visible at all.  - Riitta
        if weightLimits==None:
            if wmin==0:
                weightLimits=(wmin,wmax)
            else:
                weightLimits=(wmin-0.00001,wmax) 

        myEdgeColorMap=setColorMap(edgeColorMap)

        # Plot edges according to weight, beginning with small weight
        sortedEdges=list(net.edges)
        sortedEdges.sort(key=lambda x: x[2])
        for edge in sortedEdges:

            width=setEdgeWidth(edge[2],weightLimits,minwidth,maxwidth)
            colour=setColor(edge[2],weightLimits,myEdgeColorMap)
            xcoords=[xy[edge[0]][0],xy[edge[1]][0]]
            ycoords=[xy[edge[0]][1],xy[edge[1]][1]]
            edgeObject=plot_edge(axes, xcoords, ycoords, width=width, colour=colour,
                      symmetric=net.isSymmetric())
            if interactive:
                #0.1 is because 0 is not added as an edge. The offset of 0.1 is dropped in int()
                edgeObjectIndices[edge[0],edge[1]]=len(edgeObjects)+0.1 
                edgeObjects.append(edgeObject)


    # Then draw nodes, depending on given options showAllNodes
    # displays also nodes who do not have any edges left after
    # e.g. thresholding
    nodelist=[]
    if showAllNodes:
        nodelist = [node for node in xy.keys()]
    else:
        nodelist = [node for node in net]



    strengths = netext.strengths(net)
    maxs = max(strengths.values())
    mins = min(strengths.values())           

    if not(equalsize):
        A = (0 if maxs == mins else (maxnode-minnode)/(maxs-mins))
        B = maxnode-A*maxs    


    myNodeColorMap=setColorMap(nodeColorMap)

    # If nodes will be colored by setNodeColorsByProperty but
    # nodePropertyLimits were not given, use the true min and max
    # property values in the network.
    if setNodeColorsByProperty != None:
        if nodePropertyLimits == None:
            np_=[net.nodeProperty[setNodeColorsByProperty][node] for node in net]
            nodePropertyLimits=(min(np_),max(np_))

    for node in nodelist:

        # Define the shape of the node
        if equalshape:                    # If all the shapes are equal
            nodeshape = nodeShape         # Use a default shape
        elif len(nodeShapes) > 0:         # If nodeShapes are given
            if node in nodeShapes.keys(): # If shape of the node is provided
                nodeshape = nodeShapes[node]
            else :
                nodeshape = nodeShape
        else :
            nodeshape = nodeShape
        
        # First define size
        if equalsize:
            nodesize=nodeSize
            if (nodesize<1.0):          # hack: Himmeli wants size <1.0
                nodesize=nodesize*maxnode  # if Himmeli-type size used, scale up

        elif len(nodeSizes)>0:          # if nodeSizes are given, use it
            
            if node in nodeSizes.keys(): # if this node is in nodeSizes, use the value
                nodesize=nodeSizes[node]
            else:
                nodesize=minnode        # otherwise use min value (e.g. if nodeSizes has only k>0 nodes,
                                        # and k=0 nodes from thresholding are shown too.
        else:
            if node in net:
                nodesize=A*strengths[node]+B
            else:
                nodesize=minnode

        if node in net:
            nodestrength=strengths[node]
        else:
            # This is for nodes which appear in MST coords (i.e. have
            # zero links) and are thus not included in net, but should
            # yet be displayed when visualizing a thresholded network
            nodestrength=mins     

        # Then determine color
        if coloredNodes:
            if setNodeColorsByProperty != None:
                # If setNodeColorsByProperty is given, use it initially
                value = net.nodeProperty[setNodeColorsByProperty][node]
                color = setColor(value,nodePropertyLimits,myNodeColorMap)

            if len(nodeColors)>0:
                # If dict nodeColors is given, it overrides
                # setNodeColorsByProperty
                if not nodeColors.get(node): 
                    # If node is not contained in dict nodeColors 
                    if setNodeColorsByProperty == None:
                        # Use setNodeColorsByProperty if it was given,
                        # otherwise use nodeColor and if it is not
                        # given use white.
                        color = (nodeColor or (1,1,1))

                else:
                    # If node IS contained in dict nodeColors, use
                    # nodeColors[node]
                    ctemp = nodeColors[node]
                    if len(ctemp)==6: 
                        # Recognize as Himmeli-type string ('999999')
                        rc=float(ctemp[0:2])/99.0
                        gc=float(ctemp[2:4])/99.0
                        bc=float(ctemp[4:6])/99.0

                        # this is a stupid hack; sometimes rounding
                        # errors result in rc=1.0 + epsilon and
                        # matplotlib complains...
                        rc = max(0.0, min(rc, 1.0))
                        bc = max(0.0, min(bc, 1.0))
                        gc = max(0.0, min(gc, 1.0))

                        color = (rc,gc,bc)
                    else:
                        # Otherwise assume it is an RGB tuple
                        color = nodeColors[node] 

            elif setNodeColorsByProperty is None and nodeColor is not None:
                # If neither setNodeColorsByProperty or dict
                # nodeColors is given, but nodeColor is, use nodeColor.
                if len(nodeColor)==6:

                    rc=float(nodeColor[0:2])/99.0
                    gc=float(nodeColor[2:4])/99.0
                    bc=float(nodeColor[4:6])/99.0

                    rc = max(0.0, min(rc, 1.0))
                    bc = max(0.0, min(bc, 1.0))
                    gc = max(0.0, min(gc, 1.0))

                    color=(rc,gc,bc)
                else:
                    color=nodeColor         

            elif setNodeColorsByProperty == None:
                # Set color by node strength
                color = setColor(nodestrength,(mins,maxs),myNodeColorMap) 
        else:
            # If coloredNodes is False, use white.
            color=(1.0,1.0,1.0)

        # Move node labels slightly to the right so that they
        # don't coincide on the nodes
        nodeLabel_xOffset = (nodeLabel_xOffset or float(nodesize)/40)

        nodeObject=plot_node(axes, x=xy[node][0], y=xy[node][1], shape=nodeshape,
                  color=color, size=nodesize,edgecolor=nodeEdgeColor)

        if uselabels!='none' and (node in labels or uselabels == 'all'):
            if node in labels:
                if isinstance(labels[node],float):
                    showthislabel="%2.2f" % labels[node]
                else:
                    showthislabel=labels[node]
            else:
                showthislabel = str(node)

            #axes.annotate(showthislabel,(xy[node][0]+nodeLabel_xOffset,xy[node][1]),
            #              color=fontcolor,size=fontsize)

            nodeLabel_offset = int(np.ceil(float(nodesize)/2))+1
            nodeLabelObject=axes.annotate(showthislabel,(xy[node][0],xy[node][1]),
                          textcoords='offset points',
                          xytext=(nodeLabel_offset, nodeLabel_offset),
                          color=fontcolor,size=fontsize)
            if interactive:
                nodeLabelObjects[node]=nodeLabelObject
        elif interactive:
            nodeLabelObjects[node]=None

        if interactive:
            nodeObjects[node]=nodeObject


    xylist = xy.values()
    xlist=[]
    ylist=[]
    for elem in xylist:
        xlist.append(elem[0])
        ylist.append(elem[1])

    minx=min(xlist)
    maxx=max(xlist)
    miny=min(ylist)
    maxy=max(ylist)

    if keepAspectRatio:
        ywidth=maxy-miny
        xwidth=maxx-minx
        if xwidth!=0 and figsize[1]/float(figsize[0])>ywidth/float(xwidth): #if x is limiting
            maxy=miny+ywidth/2.+figsize[1]/float(figsize[0])*xwidth/2.
            miny=miny+ywidth/2.-figsize[1]/float(figsize[0])*xwidth/2.
            ywidth=maxy-miny
        else:
            maxx=minx+xwidth/2.+figsize[0]/float(figsize[1])*ywidth/2.
            minx=minx+xwidth/2.-figsize[0]/float(figsize[1])*ywidth/2.
            xwidth=maxx-minx

    xdelta=0.05*(maxx-minx)
    ydelta=0.05*(maxy-miny)

    if frame==True and showTicks==False:
        setp(axes,'xticks',[],'xticklabels',[],'yticks',[],'yticklabels',[])
    if not axisLimits==None:
        # If limits are given, use them, whatever the values of
        # showTicks or frame ...
        setp(axes,
             'xlim', (axisLimits[0][0],axisLimits[0][1]),
             'ylim', (axisLimits[1][0],axisLimits[1][1]))
    else:
        setp(axes,
             'xlim', (minx-xdelta,maxx+xdelta),
             'ylim', (miny-ydelta,maxy+ydelta))

    if interactive:
        # Save everything that is needed later
        thisfigure.nodeObjects=nodeObjects
        thisfigure.edgeObjectIndices=edgeObjectIndices
        thisfigure.edgeObjects=edgeObjects
        thisfigure.nodeLabelObjects=nodeLabelObjects
        thisfigure.coords=xy
        thisfigure.selectedNode=None

        def on_press(fig, event):
            """Find the right node and animate it.
            """
            #if event.inaxes != fig.axes: return
            if fig.selectedNode!=None: return
            if event.xdata==None or event.ydata==None: return 

            candidateNodes=[]
            xlim=fig.axes[0].get_xlim()
            ylim=fig.axes[0].get_ylim()
            xScalingFactor=1./72./fig.get_figwidth()*(xlim[1]-xlim[0])
            yScalingFactor=1./72./fig.get_figheight()*(ylim[1]-ylim[0])
            for node in fig.nodeObjects:
                xy=np.array([[nodeObjects[node].get_xdata()[0],nodeObjects[node].get_ydata()[0]]])
                d=(xy-np.array([event.xdata,event.ydata]))[0]
                d[0]=d[0]/xScalingFactor
                d[1]=d[1]/yScalingFactor
                d=np.sqrt(sum(d*d))
                if d<max(nodeObjects[node].get_markersize(),5.0):
                    candidateNodes.append((node,d))
            if len(candidateNodes)==0:
                return

            sd=None
            closestNode=None
            for node,d in candidateNodes:
                if sd==None or d<sd:
                    sd=d
                    closestNode=node

            fig.selectedNode=closestNode
            nodeObject=fig.nodeObjects[closestNode]
            nodeLabelObject=fig.nodeLabelObjects[closestNode]
            
            x0,y0=nodeObject.get_xdata()[0],nodeObject.get_ydata()[0]
            if nodeLabelObject!=None:
                xl0,yl0=nodeLabelObject.xy
            else:
                xl0,yl0=None,None
            
            fig.press=x0,y0,event.xdata,event.ydata,xl0,yl0

            # draw everything but the selected objects and store the pixel buffer
            canvas = fig.canvas
            axes = fig.axes[0]
            bbox=fig.bbox

            #set the node and all its edges animated
            nodeObject.set_animated(True)
            if nodeLabelObject!=None:
                nodeLabelObject.set_animated(True)
            for edgeIndex in fig.edgeObjectIndices[closestNode].weights:
                edgeObject=fig.edgeObjects[int(edgeIndex)]
                edgeObject.set_animated(True)

            canvas.draw()
            fig.background=canvas.copy_from_bbox(bbox)

            # now redraw the animated objects
            if nodeLabelObject!=None:
                axes.draw_artist(nodeLabelObject)
            for edgeIndex in fig.edgeObjectIndices[closestNode].weights:
                edgeObject=fig.edgeObjects[int(edgeIndex)]
                axes.draw_artist(edgeObject)
            axes.draw_artist(nodeObject)

            # and blit the redrawn area
            canvas.blit(bbox)

        def on_motion(fig, event):
            'on motion we will move the rect if the mouse is over us'
            #if event.inaxes != fig.axes: return
            if fig.selectedNode==None: return
            if event.xdata==None: return
            
            x0, y0, xpress, ypress, xl0,yl0 = fig.press
            dx = event.xdata - xpress
            dy = event.ydata - ypress
            
            nodeObject=fig.nodeObjects[fig.selectedNode]
            nodeLabelObject=fig.nodeLabelObjects[fig.selectedNode]
            for edgeIndex in fig.edgeObjectIndices[fig.selectedNode].weights:
                edgeObject=fig.edgeObjects[int(edgeIndex)]
                #select the right end of the edge
                if nodeObject.get_xdata()[0]==edgeObject.get_xdata()[0] and nodeObject.get_ydata()[0]==edgeObject.get_ydata()[0]:
                    thisEnd=0
                else:
                    thisEnd=1
                x,y=list(edgeObject.get_xdata()),list(edgeObject.get_ydata())
                x[thisEnd]=x0+dx
                y[thisEnd]=y0+dy
                edgeObject.set_xdata(x)
                edgeObject.set_ydata(y)
            nodeObject.set_xdata([x0+dx])
            nodeObject.set_ydata([y0+dy])
            if nodeLabelObject!=None:
                nodeLabelObject.xy=(xl0+dx,yl0+dy)

            canvas = fig.canvas
            axes = fig.axes[0]
            # restore the background region
            canvas.restore_region(fig.background)

            # redraw just the current rectangle
            for edgeIndex in fig.edgeObjectIndices[fig.selectedNode].weights:
                edgeObject=fig.edgeObjects[int(edgeIndex)]
                axes.draw_artist(edgeObject)
            axes.draw_artist(nodeObject)
            if nodeLabelObject!=None:
                axes.draw_artist(nodeLabelObject)

            # blit just the redrawn area
            canvas.blit(axes.bbox)

        def on_release(fig, event):
            'on release we reset the press data'
            if fig.selectedNode==None: return

            # turn off the rect animation property and reset the background
            nodeObject=fig.nodeObjects[fig.selectedNode]
            nodeLabelObject=fig.nodeLabelObjects[fig.selectedNode]
            nodeObject.set_animated(False)
            if nodeLabelObject!=None:
                nodeLabelObject.set_animated(False)
            for edgeIndex in fig.edgeObjectIndices[fig.selectedNode].weights:
                edgeObject=fig.edgeObjects[int(edgeIndex)]
                edgeObject.set_animated(False)

            fig.coords[fig.selectedNode]=tuple([nodeObject.get_xdata()[0],nodeObject.get_ydata()[0]])

            fig.background = None

            # redraw the full figure
            fig.canvas.draw()

            fig.press = None
            fig.selectedNode = None



        def connect(fig):
            'connect to all the events we need'
            fig.cidpress = fig.canvas.mpl_connect(
                'button_press_event', lambda event:on_press(fig,event))
            fig.cidrelease = fig.canvas.mpl_connect(
                'button_release_event', lambda event:on_release(fig,event))
            fig.cidmotion = fig.canvas.mpl_connect(
                'motion_notify_event', lambda event:on_motion(fig,event))

        def disconnect(fig):
            'disconnect all the stored connection ids'
            fig.canvas.mpl_disconnect(fig.cidpress)
            fig.canvas.mpl_disconnect(fig.cidrelease)
            fig.canvas.mpl_disconnect(fig.cidmotion)



        thisfigure.startInteraction=lambda:connect(thisfigure)
        thisfigure.stopInteraction=lambda:disconnect(thisfigure)

    return thisfigure
            
# ---------------------------------------
              
class Himmeli:
    """ This class uses the executable Himmeli, which produces an .eps
    file AND outputs x-y-coordinates of nodes for visualization.
    """

    #---------------------------------------
    
    # First we have to find this executable
    if sys.platform=='win32':
        # For Win use direct path (probably works better with the
        # executable network toolbox...)
        netext_path = os.path.dirname(netext.__file__)
        himmeliExecutableAlternatives=[os.path.join(netext_path,"himmeli_3.0.1","win32","himmeli.exe"),
                                       os.path.join(netext_path,"..","himmeli_3.0.1","win32","himmeli.exe"),
                                       os.path.abspath(os.path.join("himmeli_3.0.1","win32","himmeli.exe")),
                                       os.path.abspath(os.path.join("himmeli_3.0.1","himmeli.exe"))
                                       ]
        himmeliExecutable=None
        for alt in reversed(himmeliExecutableAlternatives):
            if os.path.isfile(alt):
                himmeliExecutable=alt
        if himmeliExecutable==None:
            raise Exception("Cannot find Himmeli executable!")
    else:
        # Trick: find out where netext.py is (must be in the netpython
        # directory), then add the rest of the path
        #himmeliExecutable = (os.path.dirname(netext.__file__)
        #                     +"/Himmeli/himmeli.exe")
        netext_path = os.path.dirname(netext.__file__) 
        himmeliExecutable = ("%s%s../himmeli_3.0.1/himmeli.exe" % 
                             (netext_path, ("/" if netext_path else "")))

        if not(os.path.isfile(himmeliExecutable)):
            # Just in case Himmeli was compiled without the .exe:
            #himmeliExecutable = (os.path.dirname(netext.__file__)
            #                     + "/Himmeli/himmeli")
            himmeliExecutable = ("%s%s../himmeli_3.0.1/himmeli" % 
                                 (netext_path, ("/" if netext_path else "")))

    # Directly complain if Himmeli not found.
    if not(os.path.isfile(himmeliExecutable)):
        complaint = ("Cannot find Himmeli! This is where it should be: "
                     + himmeliExecutable)
        raise Exception(complaint)

    # ----------------------------------------

    epsilon=0.0001 #hack, find the real thing

    def __init__(self, inputnet, time=20, configFile=None, threshold=None,
                 useMST=False, wmin=None, wmax=None, coloredNodes=True,
                 equalsize=True, nodeColor="999999", nodeSize="1.0",
                 coordinates=None, labels={}, distanceUnit=1, showAllNodes=True,
                 edgeLabels=False, nodeColors={}, treeMode=False,tempdir=None):
        """
        inputs:
        time - time limit (secs) for the Himmeli optimization of layout; for large nets, use higher values
        configFile - Himmeli .cfg file, if you want to use a pre-existing one
        threshold - for weighted nets; use only edges above this
        useMST (true/false) : true - uses pre-calculated MST coords for the (weighted) net
        tempdir - Temporary directory for Himmeli input and output files. If None, temporary
                  directory is set as tempfile.gettempdir().
        """
        if tempdir==None:
            tempdir=tempfile.gettempdir()

        # Checking that the given net is valid and not empty
        #if net.__class__!=pynet.Net and net.__class__!=pynet.SymmNet:
        if not isinstance(inputnet,pynet.VirtualNet):
            raise AttributeError("Unknown net type "+str(inputnet.__class__))
        if len(inputnet._nodes) == 0:
            raise AttributeError("The net cannot be empty.")

        # Another EDEN-specific piece of code: if the network is a
        # distance matrix/network, first transform distances to
        # weights.
        if hasattr(inputnet,'matrixtype'):
            if inputnet.matrixtype==0:
                inputnet = transforms.dist_to_weights(inputnet)

        if wmin is None:
            # Finds out smallest and largest weight
            witer = inputnet.weights.__iter__()
            wmin = wmax = witer.next()

            for wt in witer:
                if wt < wmin:
                    wmin = wt
                if wt > wmax:
                    wmax = wt

        if showAllNodes:
            # This is especially designed for thresholded nets, where
            # there may be nodes with zero links in addition to
            # disconnected components. IF coordinates have been
            # calculated for the network WITHOUT THRESHOLDING,
            # augmentNet adds every node mentioned in coordinates to
            # the network (with epsilon-weight "ghost"
            # links). Furthermore, augmentNet joins all disconnected
            # components to the largest component with ghost links.
            net = self.augmentNet(inputnet,useMST,coordinates,wmin)

        else:
            net=inputnet

        if (showAllNodes or useMST):
            threshold2 = ['abs', wmin, wmax]

        #First we need to generate names for this net and its files
        rgen = random.Random()
        netName = str(rgen.randint(1,10000))
        edgFileName = os.path.join(tempdir,"himmeli_tmp"+netName+".edg")
        confFileName = os.path.join(tempdir,"himmeli_tmp"+netName+".cfg")
        vtxFileName = os.path.join(tempdir,"himmeli_tmp"+netName+".vtx")
        coordFileName = os.path.join(tempdir,netName+".vertices.txt")
        legendFileName = os.path.join(tempdir,netName+".legend.eps")
        output_confFileName = os.path.join(tempdir,netName+'.config.txt')
        output_edgFileName = os.path.join(tempdir,netName+'.edges.txt')
        output_psFileName = os.path.join(tempdir,netName+'.ps')

        #Then we make config for Himmeli or read it from a file
        if configFile==None:
            config = ("EdgeHeadVariable\tHEAD\n"
                      "EdgeTailVariable\tTAIL\n"
                      "EdgeWeightVariable\tWEIGHT\n"
                      "FigureLimit\t1\n"
                      "DecorationMode\ton\n"
                      "IncrementMode\ton\n")
            if treeMode:
                config += "TreeMode\ton\n"
            else:
                config += "TreeMode\toff\n"
            if coordinates==None:
                config += "TimeLimit\t"+str(time)+"\n"
            else:
                config += "TimeLimit\t0\n"
            if net.isSymmetric():
                config += "ArrowMode\toff\n"
            else:
                config += "ArrowMode\ton\n"
            #config += "PageSize\tauto\tauto\n"
            config += "PageSize\tauto\ta4\n"
            config += "DistanceUnit\t"+str(distanceUnit)+"\n"
            if len(labels)>0 and edgeLabels==False:
                config += "LabelMode\tvertex\n"
            elif len(labels)>0 and edgeLabels==True:
                config += "LabelMode\ton\n"
            else:
                config += "LabelMode\toff\n"
            #config += "PageOrientation\tportrait\n"
            if threshold != None:
                # Not tested properly (and everything else is?)
                filterType=threshold[0]
                minedge=threshold[1]
                maxedge=threshold[2]
                config += "EdgeWeightFilter\t"+filterType+"\t"
                config += str(0.99*minedge)+"\t"+str(1.01*maxedge)+"\n"
            if useMST or showAllNodes:
                filterType=threshold2[0]
                minedge=threshold2[1]
                maxedge=threshold2[2]
                config += "EdgeWeightMask\t"+filterType+"\t"
                config += str(0.99*minedge)+"\t"+str(1.01*maxedge)+"\n"

        else:
            configTemplate = open(configFile)
            config = configTemplate.read()
            configTemplate.close()

        #---These are specific for this Himmeli run
        config += "GraphName\t"+netName+"\n"
        config += "EdgeFile\t"+edgFileName+"\n"

        config += "VertexFile\t"+vtxFileName+"\n"
        config += "VertexNameVariable\tVNAME"+"\n"
        config += "VertexLabelVariable\tVLABEL"+"\n"

        if coloredNodes:
            config += "VertexColorVariable\tVCOLOR"+"\n"
        if equalsize:
            config += "VertexSizeVariable\tVSIZE\n"
        if coordinates!=None:
            config += "VertexXVariable\tVX\n"
            config += "VertexYVariable\tVY\n"

            #if showAllNodes:
            #    config += ("EdgeWeightMask\tabs\t"+str(2*self.epsilon)+"\t"
            #               +str(1/self.epsilon)+"\n")
            #    config += ("EdgeWeightFilter\tabs\t"+str(2*self.epsilon)+"\t"
            #               +str(1/self.epsilon)+"\n")

        #---

        # Now we write all nessesary files for Himmeli to disk edge
        # file:
        netio.writeNet(net,edgFileName,headers=True)

        # Vertex file:
        self.writeVertexFile(net, vtxFileName, coloredNodes, equalsize,
                             nodeColor, nodeSize, coordinates=coordinates,
                             labels=labels, nodeColors=nodeColors)
        # Config file:
        confFile=open(confFileName,'w')
        confFile.write(config)
        confFile.close()

        # All is set for running Himmeli
        himmeli = subprocess.Popen([self.himmeliExecutable,confFileName],
                                   cwd=tempdir,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   creationflags=subprocess.SW_HIDE,
                                   shell=True)
        output,errors=himmeli.communicate()

        #print himmeli #for short debug

        # Save coordinates produced by Himmeli
        self.coordinates = self._parseCoordinates(coordFileName)

        # Remove files that are not needed anymore
        os.remove(edgFileName)
        os.remove(confFileName)
        os.remove(coordFileName)
        if os.path.isfile(legendFileName):
            os.remove(legendFileName)
        os.remove(output_confFileName)
        os.remove(output_edgFileName)
        os.remove(output_psFileName)
        os.remove(vtxFileName)

        # Finally we save all the information needed later
        self.netName=netName

    def getCoordinates(self):
        return self.coordinates

    def saveEps(self,filename):
        shutil.copyfile(self.netName+"_0001.eps",filename)

    def _parseCoordinates(self, coordFileName):
        coordFile = open(coordFileName, 'r')
        coordFile.readline()
        coordinates={}
        for line in coordFile:
            columns=line.split("\t")
            if len(columns) == 13:
                # Previous version:
                #[comp,name,x,y,z,degree_in,dg_out,strengt_in,strength_out,
                # color,shape,size,label]=line.split("\t")
                name, x, y, z = line.split("\t")[1:5]
                if len(name) > 0:
                    try:
                        name = int(name)
                    except ValueError:
                        pass
                    coordinates[name] = tuple(map(float, (x,y,z)))
        return coordinates
        
    def _writeEpsilonEdges(self,net,fileName):
        """Obsolete! Replaced by augmentNet, which generates a new
        network with the epsilon edges directly added. No need for
        rewriting to any file; EdgeWeightFilter is always used.
        """
      #  with open(fileName,'a') as f:
      #      last=None
      #      for node in net:
      #          if last!=None:
      #             if net[last,node]==0:
      #                 file.write(str(node)+"\t"+str(last)+"\t"
      #                            +str(self.epsilon)+"\n")
      #          last=node
        pass

    def __del__(self):
        if os.path.isfile(self.netName+"_0001.eps"):
            os.remove(self.netName+"_0001.eps")
        
    def augmentNet(self, net, useMST=False, coords=None, wmin=None):
        """Make a network singly connected.

        If the network is singly connected, returns the network
        untouched. If not, connects one node of each component to the
        giant component with a very small weight and produces a list
        of these ghost edges. Alternatively, if `useMST` is True, all
        nodes which are keys in the coordinate list are added to the
        network, joined with epsilon edges.  This is useful for
        thresholding networks with various thresholds, so that even
        all k=0 nodes are always visible. 

        Returns augmented net with epsilon edges added.
        """
        # If min weight is not given, find it out
        if wmin is None:
            witer = net.weights.__iter__()
            wmin = witer.next()
            for wt in witer:
                if wt < wmin:
                    wmin = wt

        epsilon_factor = 25.0

        init_comp = percolator.getComponents(net)
        sizelist = [len(c) for c in init_comp]

        mc = max(sizelist)
        maxcomponent_index = 0
        for i, v in enumerate(sizelist):
            if v == mc:
                maxcomponent_index = i

        giantmembers = init_comp[maxcomponent_index]

        Ncomponents = len(init_comp)
        maxcomponent = len(giantmembers)

        # MAKE A COPY OF THE ORIGINAL NETWORK;
        #
        # deepcopy generated segfaults probably due to C++ interfacing
        # so this is simple and stupid (and slow). 

        N = len(net._nodes)

        if (isinstance(net,pynet.SymmFullNet)):
            newnet=pynet.SymmFullNet(N)
        else:
            newnet=pynet.SymmNet()
    
        for n_i, n_j, w_ij in net.edges:
            newnet[n_i][n_j] = w_ij

        # If MST not used, just connect all disjoint components.

        if useMST==False:
            if Ncomponents != 1:
                # Find out which component is the giant one
                # list its members for ghost edge targets.
                gianttemp = copy.deepcopy(giantmembers)

                giantmembers = []
                for i in range(0,len(gianttemp)):
                    giantmembers.append(gianttemp.pop())

                for index,c in enumerate(init_comp):
                    if index != maxcomponent_index:
                        ghostsource = c.pop()
                        tindex = int(math.ceil(random.random()*
                                               float(len(giantmembers)-1)))
                        ghosttarget = giantmembers[tindex]
                        newnet[ghostsource][ghosttarget] = wmin/epsilon_factor

        # Next, if useMST=True, use keys of coordinates to insert all
        # original nodes to the net, again with epsilon weights

        else:
            prevNode=None
            for newNode in coords.keys():
                if prevNode is not None:
                   if newnet[prevNode][newNode] == 0:
                       # JS 290509 changes [newNode,prevNode] to [newNode][prevNode]
                       newnet[newNode][prevNode] = wmin/epsilon_factor   
                prevNode = newNode

        return newnet


    def writeVertexFile(self, net, fileName, coloredNodes=True, equalsize=True,
                        singlecolor="999999", vcolors=0, singlesize="0.3",
                        coordinates=None, labels={}, nodeColors={}):
        with open(fileName, 'w') as f:

            if len(nodeColors) > 0:
                coloredNodes = True

            headers = ["VNAME", "VLABEL"]
            if coloredNodes:
                headers.append("VCOLOR")
            if equalsize:
                headers.append("VSIZE")
            if coordinates is not None:
                headers.append("VX")
                headers.append("VY")
            f.write("\t".join(headers) + "\n")

            for i in net:
                f.write(str(i))
                try:
                    f.write("\t"+str(labels[i]))
                except KeyError:
                    f.write("\t"+str(i))
                if coloredNodes:
                    if len(nodeColors)>0:
                        f.write("\t"+str(nodeColors[i]))
                    else:
                        f.write("\t"+singlecolor)
                if equalsize:
                    f.write("\t"+singlesize)
                if coordinates!=None:
                    f.write("\t"+str(coordinates[i][0]))
                    f.write("\t"+str(coordinates[i][1]))
                f.write("\n")

            
# ---------------------------------------
    
def drawNet(net,labels={},coordinates=None,showAllNodes=False):
    """Display a picture of the network using Himmeli
    """
    from Tkinter import Tk,TOP,BOTH,YES
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

    if coordinates==None:
        h=Himmeli(net, labels=labels, coordinates=coordinates,
                  showAllNodes=showAllNodes)
        coordinates=h.getCoordinates()


    f=VisualizeNet(net,coordinates,interactive=True,uselabels='all',labels=labels)
    w=Tk()
    canvas=FigureCanvasTkAgg(f,master=w)
    canvas.show()
    canvas._tkcanvas.pack(side=TOP,fill=BOTH,expand=YES)
    f.startInteraction()


# ---------------------------------------        

# def shiftCoordinates(xy,nodelist,shift):
def shiftCoordinates(xy,nodelist, xshift=0, yshift=0, zshift=0):
    """Translate coordinates of given nodes. 
    
    Parameters
    ----------
    xy : dict
        A dictionary in which item 'node' is a tuple containing the
        coordinates of 'node'. They must contain either two or three
        elements, as in (x, y) or (x, y, z).
    nodelist : sequence
        Listing the subset of keys in xy that need to be translated
    xshift, yshift, zshift : float
        These values indicate how much to shift the coordinates. Each
        defaults to zero. If the coordinate list contains tuples of
        length two, zshift will be ignored.

    Return
    ------
    xy : dict
        The shifted coordinates.
    """
    
    for node in nodelist:
        coords = xy[node]
        if len(coords) == 2:
            xy[node] = (coords[0]+xshift, coords[1]+yshift)
        elif len(coords) == 3:
            xy[node]=(coords[0]+xshift, coords[1]+yshift, coords[2]+zshift)
        else:
            raise ValueError("The coordinate tuples must contain two or "
                             "three elements.") 
    return xy


def getWheelCoords(net, node, N_trys=1):
    """Return coordinates for a friend wheel."""
    
    def calculate_cost(loc, net):
        cost = 0
        for ni,nj,wij in net.edges:
            dist = np.abs(loc[ni] - loc[nj])
            cost += (min(dist, len(loc)-dist)-1)*wij
        return cost
            
    # Get the subnetwork spanned by the neighbours of `node`.
    neighbours = list(net[node])
    N = len(neighbours) # Includes also non-connected nodes.
    neighbour_net = transforms.getSubnet(net, neighbours)

    # There is nothing to optimize if N <= 3. Just return the obvious
    # answer.
    if N_all == 0:
        return {}

    if N >= 3:
        curr_res = (-1, {})
        for try_count in range(N_trys):
            # Go through all neighbouring nodes in a random order, switching
            # each to the best location, until no more change occurs.
            locs = dict(zip(neighbours, range(N)))
            rand_order = neighbours[:]
            changed = True
            while changed:
                changed = False
                np.random.shuffle(rand_order)
                current_cost = calculate_cost(locs, neighbour_net)
                for i, n_rand in enumerate(rand_order):
                    # Find the best location for n_rand.
                    current_loc = locs[n_rand]
                    best = (current_cost, n_rand)
                    for other in rand_order[:i]+rand_order[i+1:]:
                        new_locs = locs.copy()
                        new_locs[n_rand], new_locs[other] = new_locs[other], new_locs[n_rand]
                        new_cost = calculate_cost(new_locs, neighbour_net)
                        if new_cost < best[0]:
                            best = (new_cost, other)
                    if best[0] < current_cost:
                        locs[n_rand], locs[best[1]] = locs[best[1]], locs[n_rand]
                        current_cost = best[0]
                        changed = True

            if current_cost < curr_res[0] or curr_res[0] == -1:
                curr_res = (current_cost, locs)

        locs = curr_res[1] # Use the best result.

    # Optimal configuration found, calculate the coordinates on a
    # circle.
    coords = {}
    for n, loc in locs.iteritems():
        coords[n] = (np.cos(2*np.pi*loc/N_all), np.sin(2*np.pi*loc/N_all), 0.0)

    return coords

if __name__ == '__main__':
    """Run unit tests if called."""
    from tests.test_visuals import *
    unittest.main()
