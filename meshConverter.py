#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 14:06:41 2020

@author: carlos
"""


import gmsh
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

"""
Line elements:
    type, phisical marker, vertice1, vertice2
    
Tri elements:
    type, phisical marker, vertice1, vertice2, vertice3

"""


def calculateNodes():
    
    nodeTags, nodeCoords, _ = gmsh.model.mesh.getNodes()        
    
    nodesAux = zip(nodeCoords[0::3], nodeCoords[1::3], nodeCoords[2::3])

    N = len(nodeTags)
    nodes = np.zeros((N, 3))
    for ii in range(0, N):
        nodes[ii, 0] = nodesAux[ii][0]
        nodes[ii, 1] = nodesAux[ii][1]
        nodes[ii, 2] = nodesAux[ii][2]

    return nodes

class block():
    
    def __init__(self, eList, n):
        
        self.eList = eList
        self.n = n
        self.N = 0
        self.markerDict = dict()
        self._calculate()
        
        
    def _calculate(self):
        
        elemTags = []
        elemNodeTags = []
        elemP = []
        self.elemMarkers = []
        
        for e in self.eList:
            elemTypes, elemTags0, elemNodeTags0 = gmsh.model.mesh.getElements(e[0], e[1])
            elemTags += elemTags0[0]
            elemNodeTags += elemNodeTags0[0]
            marker = gmsh.model.getPhysicalGroupsForEntity(e[0], e[1])[0]
            markerName = gmsh.model.getPhysicalName(e[0], marker)
            
            if (marker, markerName) not in self.elemMarkers:
                self.elemMarkers.append((marker, str(markerName)))
                    
            self.markerDict[str(markerName)] = marker
            
            for ii in range(0, len(elemTags0[0])):
                elemP.append(marker)

        self.N = len(elemTags)

        self.elem = np.zeros((self.N, self.n), dtype=long)
        
        for ii in range(0, self.N):
            self.elem[ii, 0] = self.n
            self.elem[ii, 1] = elemP[ii]
            self.elem[ii, 2] = elemNodeTags[0 + (self.n-2)*ii]
            self.elem[ii, 3] = elemNodeTags[1 + (self.n-2)*ii]
            if self.n > 4:
                self.elem[ii, 4] = elemNodeTags[2 + (self.n-2)*ii]
            if self.n > 5:
                self.elem[ii, 5] = elemNodeTags[3 + (self.n-2)*ii]
                
        return None

    def printBlock(self, f):
        
        for ii in range(0, self.N):
            f.write((' %li,'*self.n + '\n') % tuple(self.elem[ii, :]))
            
        return None


if __name__=='__main__':

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    fileName = os.path.abspath(sys.argv[1])
    gmsh.open(fileName)
    
    print("Model name: " + gmsh.model.getCurrent())
    
    # get all elementary entities in the model
    entities = gmsh.model.getEntities()
    
    e1List = []
    e2List = []
    
    if not entities:
        print("Error: No entities.")
        raise
    
    for e in entities:
        aux = gmsh.model.getType(e[0], e[1])
        print("Entity " + str(e) + " of type " + aux)
        if e[0] == 1:
            e1List.append(e)
        elif e[0] == 2:
            e2List.append(e)

    blk1 = block(e1List, 4)
    blk2 = block(e2List, 5)
    
        
    elemMarkers = blk1.elemMarkers + blk2.elemMarkers
    print(elemMarkers)
    
    nodes = calculateNodes()

    gmsh.finalize()
    
    # Open output file
    out = open(fileName+"sim", 'w')
    
    # Print Numbers
    N = blk1.elem.shape[0]+blk2.elem.shape[0]
    out.write('%i, %i, %i,\n' % (len(elemMarkers), nodes.shape[0], N))
    
    # Print markers and markers names
    for ii in range(0, len(elemMarkers)):
        out.write('%i, %s,\n' % elemMarkers[ii])

    # Print nodes
    for ii in range(0, nodes.shape[0]):
        out.write('%16.10lf, %16.10lf,\n' % (nodes[ii][0], nodes[ii][1]))

    # Print line elements
    blk1.printBlock(out)
    # Print tri elements
    blk2.printBlock(out)

    # Close output file
    out.close()
