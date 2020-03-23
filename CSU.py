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

class material():
    
    def __init__(self):
        
        self.rho = 0
        self.Cp = 0
        self.k = 0
        

class configuration():

    def __init__(self, inputName):
        
        self.inputName = os.path.abspath(inputName)
        self.config = configparser.ConfigParser()
        self.config.optionxform = str
        
        # Stefan-Boltzmann constant
        self.sig =  5.670374419e-8

        if not os.path.isfile(inputName):
            
            raise Exception('This .csi file does not exist')
                
        self.config.read(self.inputName)        
        self._read()
        

    def _read(self):
        
        self.var = dict()
                
        head, tail = os.path.split(self.inputName)
        self.var['head'] = head
        self.var['tail'] = tail
        flagMaterial = False

        # Read input values
        
        self.section = 'input'
        
        floatList = ['t', 'dt']
        
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
        self.var['fInput1'] = []
        self.var['fInput4'] = []
        self.var['types'] = []
        self.var['material'] = dict()
        self.var['markerMat'] = []

        aux = ''
        for section in self.config.sections():
            if "material" in section:
                flagMaterial = True
                aux = section.split('=')
                name = str(aux[1]).strip()
                
                self.var['material'][name] = dict()
                                
                self.var['material'][name]['rho'] = self.config.getfloat(section, 'rho')
                self.var['material'][name]['Cp'] = self.config.getfloat(section, 'Cp')
                self.var['material'][name]['k'] = self.config.getfloat(section, 'k')
          
        print(self.var['material'])
        for section in self.config.sections():        
        
            if "boundary" in section:
                aux = section.split('=')
                self.var['markers'].append(int(aux[1]))
                typ = self.config.get(section, 'type')
                typ = typ.strip()
                #print(typ)
                
                if typ == 'constant temperature':
                    self.var['types'].append('T')
                    inp = self.config.getfloat(section, 'temperature')
                    self.var['fInput'].append(inp)
                    self.var['fInput1'].append(0.0)
                    self.var['fInput4'].append(0.0)
                    self.var['markerMat'].append(None)
                    
                elif typ == 'heat flux':
                    
                    heatDict = {'heat flux' : 0.0, 
                                'convection coefficient' : 0.0,
                                'recuperation temperature' : 0.0, 
                                'emissivity' : 0.0, 
                                'environment temperature' : 0.0}
                    
                    self.var['types'].append('Q')

                    inputDict = dict(self.config.items(section))
                    #print(inputDict)
                    #print(heatDict)
                    for inp in inputDict.keys():
                        if inp != 'type':
                            heatDict[inp] = float(inputDict[inp])

                    q = heatDict['heat flux']
                    h = heatDict['convection coefficient']
                    Tr = heatDict['recuperation temperature']
                    e = heatDict['emissivity']
                    Te = heatDict['environment temperature']

                    self.var['fInput'].append(q + h*Tr + e*self.sig*(Te**4))
                    self.var['fInput1'].append(h)
                    self.var['fInput4'].append(e*self.sig)
                    self.var['markerMat'].append(None)

                elif typ == 'domain':
                    self.var['types'].append('D')
                    inp = self.config.getfloat(section, 'initial temperature')
                    self.var['fInput'].append(inp)
                    self.var['fInput1'].append(0.0)
                    self.var['fInput4'].append(0.0)
                    
                    if flagMaterial:
                        inp = self.config.get(section, 'material')
                        self.var['markerMat'].append(self.var['material'][inp])
                    else:
                        self.var['markerMat'].append(None)

        #print(self.var['markerMat'])
                    
        self.var['Nmarkers'] = len(self.var['markers'])
        
        self.var['C'] = []
        self.var['K'] = []
        if flagMaterial:
            for ii in range(0, self.var['Nmarkers']):
                if self.var['types'][ii] == 'D':
                    rho = self.var['markerMat'][ii]['rho']
                    Cp = self.var['markerMat'][ii]['Cp']
                    self.var['C'].append(rho*Cp/self.var['dt'])
                    self.var['K'].append(self.var['markerMat'][ii]['k'])
                else:
                    self.var['C'].append(0.0)
                    self.var['K'].append(0.0)     
                
        self.var['meshName'] = self.var['head']+'/'+self.var['meshName']        
        self.var['outName'] = self.var['head']+'/'+self.var['outName']
        
        #self.var['c'] = self.var['rho']*self.var['Cp']/self.var['dt']
        self.var['N'] = int(self.var['t']/self.var['dt'])
        self.var['Lref'] = 1 # Deprecated. Must be removed.

        return None

    def printVar(self):
        
        for n in self.var.keys():
            print(n, ' = ', self.var[n])
            
        return None

    def printVarAuxFile(self):
        
        #print('c, %.10f,' % (self.var['c']))
        print('vague,')
        print('N, %i,' % (self.var['N']))
        print('Nmarkers, %i,' % (self.var['Nmarkers']))
        print(('markers,'+self.var['Nmarkers']*' %i,'+'') % 
                tuple(self.var['markers']))
        print(('fInput,'+self.var['Nmarkers']*' %.10f,'+'') % 
                tuple(self.var['fInput']))
        print(('fInput1,'+self.var['Nmarkers']*' %.10f,'+'') % 
                tuple(self.var['fInput1']))
        print(('types,'+self.var['Nmarkers']*' %s,'+'') % 
                tuple(self.var['types']))
        print('meshName, %s,' % (self.var['meshName']))
        print('Nsave, %i,' % (self.var['Nsave']))
        print('outName, %s,' % (self.var['outName']))
        print('Lref, %f,' % (self.var['Lref']))
        #print('k, %.10f,' % (self.var['k']))
        print('vague,')
        print(('fInput4,'+self.var['Nmarkers']*' %.10f,'+'') % 
                tuple(self.var['fInput4']))
        
        print(('C,'+len(self.var['C'])*' %.10f,'+'') % 
                tuple(self.var['C']))
        print(('K,'+len(self.var['K'])*' %.10f,'+'') % 
                tuple(self.var['K']))


            
        return None
    
    def writeAuxFile(self):
        
        f = open('aux.csv', 'w')
        
        #f.write('c, %.10f,\n' % (self.var['c']))
        f.write('vague,\n')
        f.write('N, %i,\n' % (self.var['N']))
        f.write('Nmarkers, %i,\n' % (self.var['Nmarkers']))
        f.write(('markers,'+self.var['Nmarkers']*' %i,'+'\n') % 
                tuple(self.var['markers']))
        f.write(('fInput,'+self.var['Nmarkers']*' %.10f,'+'\n') % 
                tuple(self.var['fInput']))
        f.write(('types,'+self.var['Nmarkers']*' %s,'+'\n') % 
                tuple(self.var['types']))
        f.write('meshName, %s,\n' % (self.var['meshName']))
        f.write('Nsave, %i,\n' % (self.var['Nsave']))
        f.write('outName, %s,\n' % (self.var['outName']))
        #f.write('Lref, %f,\n' % (self.var['Lref']))
        f.write('vague,\n')
        f.write(('fInput1,'+self.var['Nmarkers']*' %.10f,'+'\n') % 
                tuple(self.var['fInput1']))
        #f.write('k, %.10f,\n' % (self.var['k']))
        f.write('vague,\n')
        f.write(('fInput4,'+self.var['Nmarkers']*' %.10f,'+'\n') % 
                tuple(self.var['fInput4']))        
        
        f.write(('C,'+len(self.var['C'])*' %.10f,'+'\n') % 
                tuple(self.var['C']))
        f.write(('K,'+len(self.var['K'])*' %.10f,'+'\n') % 
                tuple(self.var['K']))

        
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
                auxSol[iii] = float(aux[0])
                
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
        
        path = "./testCases/case1/case1.csi"
        
    conf1 = configuration(path)
    
    print("\nCSU - Code of Simulation for Unstructured meshs:")    
    conf1.printVar()
    print("\nAux File:")
    conf1.printVarAuxFile()
    conf1.writeAuxFile()
    print()

    sp.call("./bin/Debug/CSU", shell=True)
    
    out1 = outputClass(conf1.var['meshName'], conf1.var['outName'], conf1.var)
    out1.plot(grid=False)


