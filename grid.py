# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 16:55:17 2018

@author: HOS
"""
from numpy import *

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

    
    """def getProt(self):
        return 
    def setGrid(self):
        return"""