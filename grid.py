# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 16:55:17 2018

@author: HOS
"""
from numpy import *
import matplotlib.pyplot as plt

class Grid:
    #constructor
    def __init__(self,p_length,g_size, prot1):
        self.N = g_size
        self.grid = zeros((self.N,self.N),dtype = int)
        self.prot = prot1
    #methods
    def setProt(self,prot1):
        self.prot = prot1
        #if size too large, give error
    def update(self):
        self.grid[self.grid > 0] = 0
        pos = self.prot.posArray
        for i in range(len(pos[0])):
            self.grid[pos[0,i],pos[1,i]] = i+1
    
    def getGrid(self):
        return self.grid
    
    def showOff(self):
        plt.plot()
        X, Y = grid.nonzero()  # Gir to matriser som sammen gir indeksene til non-zero elements
        # X=[4 4 4 5 5 6 6 6] Y=[4 5 6 4 6 4 5 6] (origo oppe til venstre)

        # Any or all of X, Y, s, and c may be masked arrays, in which case all masks will be combined and only unmasked points will be plotted.

        plt.scatter(Y, X, c=grid[X, Y], s=100, marker="o",
                    cmap=plt.cm.viridis)  # Marker size is scaled by s and marker color is mapped to c.
        plt.axis([0, arr.shape[1], 0, grid.shape[0]])
        plt.gca().invert_yaxis()
        # plt.rcParams['axes.axisbelow'] = True
        plt.show()
        
    def easyPlot(self):
        plt.imshow(self.getGrid());
        #plt.colorbar()
        plt.show()
        

    
    """def getProt(self):
        return 
    def setGrid(self):
        return"""
