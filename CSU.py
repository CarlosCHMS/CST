# Programed in Python3 but run in python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 15:39:49 2020

@author: carlos
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import configparser
import sys
import subprocess as sp
import os


class configuration():

    def __init__(self, inputName):
        
        self.inputName = os.path.abspath(inputName)
        self.config = configparser.ConfigParser()
        self.config.optionxform = str
        self.config.read(self.inputName)
        
        self._read()
        

    def _read(self):
        
        self.var = dict()
                
        head, tail = os.path.split(self.inputName)
        self.var['head'] = head
        self.var['tail'] = tail

        # Read input values
        
        self.section = 'input'
        
        floatList = ['rho', 'Cp', 'k', 't', 'dt', 'Tref', 'Lref']
        
        for n in floatList:
            self.var[n] = self.config.getfloat(self.section, n)
            
        intList = ['Nsave']

        for n in intList:
            self.var[n] = self.config.getint(self.section, n)

        strList = ['meshName', 'outName']

        for n in strList:
            self.var[n] = self.config.get(self.section, n)

        # Read boundary values

        self.var['markers'] = []
        self.var['fInput'] = []
        self.var['types'] = []

        aux = ''
        for section in self.config.sections():
            if "boundary" in section:
                aux = section.split('=')               
                self.var['markers'].append(int(aux[1]))
                typ = self.config.get(section, 'type')
                typ = typ.strip()
                print(typ)
                if typ == 'constant temperature':
                    self.var['types'].append('T')
                    inp = self.config.getfloat(section, 'temperature')
                    self.var['fInput'].append(inp)
                elif typ == 'constant heat flux':
                    self.var['types'].append('Q')
                    inp = self.config.getfloat(section, 'heat flux')
                    self.var['fInput'].append(inp)                    
                elif typ == 'domain':
                    self.var['types'].append('D')
                    inp = self.config.getfloat(section, 'initial temperature')
                    self.var['fInput'].append(inp)
                    
        self.var['Nmarkers'] = len(self.var['markers'])
        
        self.var['meshName'] = self.var['head']+'/'+self.var['meshName']        
        self.var['outName'] = self.var['head']+'/'+self.var['outName']
        

        return None

    def printVar(self):
        
        for n in self.var.keys():
            print(n, ' = ', self.var[n])
            
        return None
    
    def dimensionless(self):
        
        self.var['c'] = self.var['k']*self.var['dt']/(self.var['rho']*self.var['Cp']*(self.var['Lref']**2))
        self.var['N'] = int(self.var['t']/self.var['dt'])
        self.var['saveStep'] = int(self.var['N']/self.var['Nsave'])
        
        qref = (self.var['rho']*self.var['Cp']*self.var['Tref']*self.var['Lref']/self.var['dt'])
        
        for ii in range(0, self.var['Nmarkers']):
            if self.var['types'][ii] == 'Q':
                self.var['fInput'][ii] /= qref
            else:
                self.var['fInput'][ii] /= self.var['Tref']
                              
        return None

    def printVarAuxFile(self):
        
        print('c, %.10f,' % (self.var['c']))
        print('N, %i,' % (self.var['N']))
        print('Nmarkers, %i,' % (self.var['Nmarkers']))
        print(('markers,'+self.var['Nmarkers']*' %i,'+'') % 
                tuple(self.var['markers']))
        print(('fInput,'+self.var['Nmarkers']*' %.10f,'+'') % 
                tuple(self.var['fInput']))
        print(('types,'+self.var['Nmarkers']*' %s,'+'') % 
                tuple(self.var['types']))
        print('meshName, %s,' % (self.var['meshName']))
        print('saveStep, %i,' % (self.var['saveStep']))
        print('outName, %s,' % (self.var['outName']))
        print('Lref, %f,' % (self.var['Lref']))       
            
        return None
    
    def writeAuxFile(self):
        
        f = open('aux.csv', 'w')
        
        f.write('c, %.10f,\n' % (self.var['c']))
        f.write('N, %i,\n' % (self.var['N']))
        f.write('Nmarkers, %i,\n' % (self.var['Nmarkers']))
        f.write(('markers,'+self.var['Nmarkers']*' %i,'+'\n') % 
                tuple(self.var['markers']))
        f.write(('fInput,'+self.var['Nmarkers']*' %.10f,'+'\n') % 
                tuple(self.var['fInput']))
        f.write(('types,'+self.var['Nmarkers']*' %s,'+'\n') % 
                tuple(self.var['types']))
        f.write('meshName, %s,\n' % (self.var['meshName']))
        f.write('saveStep, %i,\n' % (self.var['saveStep']))
        f.write('outName, %s,\n' % (self.var['outName']))
        f.write('Lref, %f,\n' % (self.var['Lref']))
        
        f.close()        
        
        return None

class outputClass():
    
    def __init__(self, gridFileName, outFileName, var):
        
        
        self.var = var
        self.Nm = 0
        self.Np = 0
        self.Ne = 0
        self.triangles = []
        self.type2 = 5
        aloc = True
        
        # Read grid data
        f1 = open(gridFileName, 'r')
        ii = 0
        for row in f1:
            if ii == 0:
                aux = row.split(',')
                self.Nm = int(aux[0])
                self.Np = int(aux[1])
                self.Ne = int(aux[2])
            elif ii < (self.Nm + 1):
                if aloc:
                    self.aloc()
                    aloc = False

                iii = ii - 1
                aux = row.split(',')
                self.marker[iii] = float(aux[0])
                
            elif ii < (self.Nm + self.Np + 1):
                    
                iii = ii - (self.Nm + 1)
                aux = row.split(',')
                self.x[iii] = float(aux[0])
                self.y[iii] = float(aux[1])
                
            elif ii < (self.Nm + self.Np + self.Ne + 1):
                
                iii = ii - (self.Nm + self.Np + 1)
                aux = row.split(',')
                for jj in range(0, len(aux)-1):
                    self.elem[iii][jj] = int(aux[jj])
                
            
            ii += 1
        
        f1.close()
           
        # Read output data
        f2 = open(outFileName, 'r')
        ii = 0
        self.Ns = []
        self.solution = []
        auxSol = np.zeros(self.Np, dtype="float")
        for row in f2:
            if ii == 0:
                aux = row.split(',')
                Np2 = int(aux[0])
                if Np2 != self.Np:
                    print("Grid and output incompatibles.")
                self.Ns.append(int(aux[1]))
                
            elif ii < (self.Np + 1):
                iii = ii - 1
                aux = row.split(',')
                auxSol[iii] = float(aux[0])*self.var['Tref']
                
                if ii == self.Np:
                    self.solution.append(auxSol.copy())
                    ii = -1

            ii += 1

        if len(self.Ns) > 10:
            print("Too much graphics")
            raise        
            
        self._calcTriangles()
        
        self.triangulation = mtri.Triangulation(self.x, self.y, self.triangles)
        
        f2.close()
                
    def aloc(self):
        
        self.marker = np.zeros(self.Nm, dtype="int")
        self.x = np.zeros(self.Np, dtype="float")
        self.y = np.zeros(self.Np, dtype="float")
        self.elem = np.zeros((self.Ne, 5), dtype="int")
        
        return None
            
    def printNumbers(self):
        
        print("Numbers")
        print("%i, %i, %i," % (self.Nm, self.Np, self.Ne))
            
        return None

    def printMarkers(self):
        
        print("Markers:")
        for ii in range(0, self.Nm):
            print("%i," % (self.marker[ii]))
            
        return None
    
    def printCoords(self):
        
        print("Coordinates:")
        for ii in range(0, self.Np):
            print("%10.6f, %10.6f," % (self.x[ii], self.y[ii]))
            
        return None

    def printElements(self):
        
        print("Elements:")
        for ii in range(0, self.elem.shape[0]):
            aux = ''
            for jj in range(0, self.elem[ii, 0]):
                aux += " %5i," % (self.elem[ii, jj])
            print(aux)
            
        return None
    
    def plotPoints(self):
        
        for ii in range(0, self.Np):
        
            plt.text(self.x[ii], self.y[ii], "%i" % (ii))
        
        return None
       
    def _calcTriangles(self):
    
        for ii in range(0, self.Ne):
            
            if self.elem[ii][0] == self.type2:
                
                self.triangles.append([self.elem[ii][2]-1, self.elem[ii][3]-1, 
                                       self.elem[ii][4]-1])

        return None                
        
    def plot(self, grid=False):
        
        TmaxList = []
        TminList = []
        for ii in range(0, len(self.Ns)):
            TmaxList.append(max(self.solution[ii]))
            TminList.append(min(self.solution[ii]))
            
        Tmin = min(TminList)
        Tmax = max(TmaxList)
        for ii in range(0, len(self.Ns)):    
            self.plotii(ii, Tmin, Tmax, grid)
        
        return None
    
    def plotii(self, ii, Tmin, Tmax, grid):
        
        plt.figure()
        if grid:
            plt.triplot(self.triangulation, 'k-')
            
        plt.title("time: %f [s]" % (self.Ns[ii]*self.var['dt']))
        plt.tricontourf(self.triangulation, self.solution[ii], 
                        vmin=Tmin, vmax=Tmax)
    
        cbar = plt.colorbar()
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel('Temperature [K]', rotation=270)
    
        plt.axis("equal")
        plt.show()
        
        return None
    

if __name__=="__main__":

    External = True
    
    if External:
    
        if len(sys.argv) < 1:
            print("Insert the input file")
            exit(0)
        
        path = sys.argv[1]
        
    else:
        
        path = "./cases/case7/test7.cs"
        #path = "./test7.cs"
        
    conf1 = configuration(path)
    
    print("\nCSU - Code of Simulation for Unstructured meshs:")    
    conf1.printVar()
    conf1.dimensionless()
    print("\nAux File:")
    conf1.printVarAuxFile()
    conf1.writeAuxFile()
    print()
       
    sp.call("./bin/Debug/CSU", shell=True)
    
    out1 = outputClass(conf1.var['meshName'], conf1.var['outName'], conf1.var)
    out1.plot(grid=False)


