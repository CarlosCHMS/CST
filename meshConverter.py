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


class triangleBlock():

    def __init__(self, eSurf):
        
        self.eSurf = eSurf
        
    def calculateInterior(self):
    
        self.elemPT = []
        
        eSurf = self.eSurf
        
        elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(eSurf[0], eSurf[1])    
        
        p = gmsh.model.getPhysicalGroupsForEntity(eSurf[0], eSurf[1])[0]
        self.elemPT.append(p)
        
        pn = gmsh.model.getPhysicalName(eSurf[0], p)
        print((str(pn), (eSurf[0], eSurf[1]), p))
         
        
        #boundary = gmsh.model.getBoundary(elemTags[0])
        #n1 = dict(zip(elemTags, elemNodeTags[0][0::3]))
        elem = np.zeros((len(elemTags[0]), 8), dtype=long)
        
        for ii in range(0, len(elemTags[0])):
            elem[ii, 0] = elemTags[0][ii]
            elem[ii, 1] = elemNodeTags[0][0 + 3*ii]
            elem[ii, 2] = elemNodeTags[0][1 + 3*ii]
            elem[ii, 3] = elemNodeTags[0][2 + 3*ii]
        
            elem[ii, 7] = gmsh.model.getPhysicalGroupsForEntity(eSurf[0], eSurf[1])[0]
        
    #    aux = []
    #    for ii in range(0, elem.shape[0]):
    #        aux.append([])
        
        for ii in range(0, elem.shape[0]):
            for jj in range(0, elem.shape[0]):
                if ii != jj:
                    isBou1 = False
                    isBou2 = False
                    isBou3 = False
                    if elem[ii, 1] == elem[jj, 1]:
                        isBou1 = True
                    elif elem[ii, 1] == elem[jj, 2]:
                        isBou1 = True
                    elif elem[ii, 1] == elem[jj, 3]:
                        isBou1 = True
                    
                    if elem[ii, 2] == elem[jj, 1]:
                        isBou2 = True
                    elif elem[ii, 2] == elem[jj, 2]:
                        isBou2 = True
                    elif elem[ii, 2] == elem[jj, 3]:            
                        isBou2 = True
                        
                    if isBou1 and isBou2:
                        elem[ii, 4] = elemTags[0][jj]
                        
                    elif isBou1 or isBou2:
                        
                        if elem[ii, 3] == elem[jj, 1]:
                            isBou3 = True
                        elif elem[ii, 3] == elem[jj, 2]:
                            isBou3 = True
                        elif elem[ii, 3] == elem[jj, 3]:
                            isBou3 = True
                            
                        if isBou2 and isBou3:
                            elem[ii, 5] = elemTags[0][jj]
                        if isBou3 and isBou1:
                            elem[ii, 6] = elemTags[0][jj]
    
        self.elem = elem
        return None

    def cloupleLine(self, line):
        
        for ii in range(0, self.elem.shape[0]):
            for jj in range(0, line.shape[0]):
                    isBou1 = False
                    isBou2 = False
                    isBou3 = False
                    if self.elem[ii, 1] == line[jj, 1]:
                        isBou1 = True
                    elif self.elem[ii, 1] == line[jj, 2]:
                        isBou1 = True
                    
                    if self.elem[ii, 2] == line[jj, 1]:
                        isBou2 = True
                    elif self.elem[ii, 2] == line[jj, 2]:
                        isBou2 = True
                        
                    if isBou1 and isBou2:
                        self.elem[ii, 4] = line[jj, 0]
                        line[jj, 3] = self.elem[ii, 0]
                        
                    elif isBou1 or isBou2:
                        
                        if self.elem[ii, 3] == line[jj, 1]:
                            isBou3 = True
                        elif self.elem[ii, 3] == line[jj, 2]:
                            isBou3 = True
                            
                        if isBou2 and isBou3:
                            self.elem[ii, 5] = line[jj, 0]
                            line[jj, 3] = self.elem[ii, 0]
                            
                        if isBou3 and isBou1:
                            self.elem[ii, 6] = line[jj, 0]
                            line[jj, 3] = self.elem[ii, 0]

        return line

    def calcEdge(self):
        
        self.edgeList = []
        for ii in range(0, self.elem.shape[0]):
            a = self.elem[ii, 1]
            b = self.elem[ii, 2]
            c = self.elem[ii, 3]
            
            t1 = self.orderTuple(a, b)
            t2 = self.orderTuple(a, c)
            t3 = self.orderTuple(b, c)
            
            if t1 not in self.edgeList:
                self.edgeList.append(t1)
            if t2 not in self.edgeList:
                self.edgeList.append(t2)
            if t3 not in self.edgeList:
                self.edgeList.append(t3)

        return None

    def orderTuple(self, a, b):

        if a > b:
            ans = (b, a)
        else:
            ans = (a, b)
            
        return ans

    def plotMesh(self, nodes):
        
        for e in self.edgeList:
            x1 = nodes[e[0]-1, 0]
            x2 = nodes[e[1]-1, 0]
            y1 = nodes[e[0]-1, 1]
            y2 = nodes[e[1]-1, 1]
            
            plt.plot([x1, x2], [y1, y2])
            
    def plotNodes(self, nodes):
        
        for ii in range(0, len(nodes)):
            x1 = nodes[ii, 0]
            y1 = nodes[ii, 1]
            
            
            plt.scatter(x1, y1)
            plt.text(x1, y1, str(ii+1))
        plt.axis('equal')
        

class lineBlock():
    
    def __init__(self, eList):
        
        self.eList = eList
        
    def calculateInterior(self):
        
        elemTags = []
        elemNodeTags = []
        elemP = []
        self.elemPT = []
        elemPN = []
        
        #print(dir(gmsh.model))
        #print(gmsh.model.getPhysicalName(1, 1))
        
        for e in self.eList:
            elemTypes, elemTags0, elemNodeTags0 = gmsh.model.mesh.getElements(e[0], e[1])
            elemTags += elemTags0[0]
            elemNodeTags += elemNodeTags0[0]
            p = gmsh.model.getPhysicalGroupsForEntity(e[0], e[1])[0]
            if p not in self.elemPT:
                self.elemPT.append(p)
            pn = gmsh.model.getPhysicalName(e[0], p)
            elemPN.append((str(pn), (e[0], e[1]), p))
            for ii in range(0, len(elemTags0[0])):
                elemP.append(p)
        print(elemPN)
        #raise

        self.elem = np.zeros((len(elemTags), 5), dtype=long)
        
        for ii in range(0, len(elemTags)):
            self.elem[ii, 0] = elemTags[ii]
            self.elem[ii, 1] = elemNodeTags[0 + 2*ii]
            self.elem[ii, 2] = elemNodeTags[1 + 2*ii]
            self.elem[ii, 4] = elemP[ii]
           
        return None


if __name__=='__main__':

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    #gmsh.open(sys.argv[1])
    fileName = os.path.abspath(sys.argv[1])
    gmsh.open(fileName)
    
    print("Model name: " + gmsh.model.getCurrent())
    
    # get all elementary entities in the model
    entities = gmsh.model.getEntities()
    
    e1List = []
    elemPT = []
    
    if not entities:
        print("Error: No entities.")
        raise
    
    for e in entities:
        aux = gmsh.model.getType(e[0], e[1])
        print("Entity " + str(e) + " of type " + aux)
        if e[0] == 1:
            e1List.append(e)
        elif e[0] == 2:
            eSurf = e

    blk1 = lineBlock(e1List)
    blk1.calculateInterior()
    

    blk2 = triangleBlock(eSurf)
    blk2.calculateInterior()
    
    
    elemPT = blk1.elemPT + blk2.elemPT
    print(elemPT)
    
    blk1.elem = blk2.cloupleLine(blk1.elem)
            
    nodes = calculateNodes()

    blk2.calcEdge()
    
    #plt.figure()
    #blk2.plotNodes(nodes)
    #blk2.plotMesh(nodes)
    #plt.show()
    
    gmsh.finalize()
    
    out = open(fileName+"sim", 'w')
    
    N = blk1.elem.shape[0]+blk2.elem.shape[0]
    out.write('%i, %i, %i,\n' % (len(elemPT), nodes.shape[0], N))
    for ii in range(0, len(elemPT)):
        out.write('%i,\n' % (elemPT[ii]))

    for ii in range(0, nodes.shape[0]):
        out.write('%16.10lf, %16.10lf,\n' % (nodes[ii][0], nodes[ii][1]))

    for ii in range(0, blk1.elem.shape[0]):
        out.write('4, %li, %li, %li,\n' % (blk1.elem[ii, 4], blk1.elem[ii, 1], blk1.elem[ii, 2]))
    for ii in range(0, blk2.elem.shape[0]):
        out.write('5, %li, %li, %li, %li,\n' % (blk2.elem[ii, 7], blk2.elem[ii, 1], blk2.elem[ii, 2], blk2.elem[ii, 3]))

    out.close()
