# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 08:39:17 2018

@author: HOS
"""
from numpy import *
#from random import randint
import matplotlib.pyplot as plt

class Protein:
    #methods
    def __init__(self,p_length = 0, g_size = 0):
        self.d= int(p_length)
        self.N= int(g_size)
        self.node = self.d//2
        self.posArray = array([[self.N/2,i] for i in range((self.N-self.d)//2, (self.N+self.d)//2)],dtype = int) 
    
    def pRotate(self,pivot, direction):
        if direction == 'counterwise':
            A = array([[0, -1], [1, 0]])
        elif direction == 'clockwise':
            A = array([[0,1],[-1,0]])
        #deler opp posArray i tre deler, prepivot[0], pivot[1] og postpivot[2]
        temp = vsplit(self.posArray,array([pivot-1,pivot]))
        if pivot > self.node: #pivotmonomer større enn noden
            temp[2] -= temp[1]
            print(temp[2])
            #hvis pivot høyere enn noden, roterer ---*--- <-'høyre-hale'
            temp[2] = (A@temp[2].transpose()).transpose()
            temp[2] += temp[1]
        elif pivot < self.node: #pivotmonomer mindre enn noden
            temp[0] -= temp[1]
            temp[0] = (A@temp[0].transpose()).transpose()
            temp[0] += temp[1]
        self.posArray = vstack([temp[0], temp[1], temp[2]])
        
    def getPos(self):
        return self.posArray
    

class Grid:
    #constructor
    def __init__(self,p_length,g_size):
        self.N = g_size
        self.grid = zeros((N,N),dtype = int)
    #methods
    def setProt(self,prot1):
        self.prot = prot1
        #if size too large, give error
    def updateGrid(self):
        self.grid[self.grid > 0] = 0
        pos = self.prot.posArray
        for i in range(len(pos)):
            self.grid[pos[i,0],pos[i,1]] = i+1
    """def getProt(self):
        return 
    def setGrid(self):
        return"""
    def getGrid(self):
        return self.grid

#utility
def d(a):
    if a == 1:
        return 'counterwise'
    return 'clockwise'

N = 12
d = 10
Prt = Protein(d,N)
#print(Prt.posArray)
Prt.pRotate(6,'counterwise')
Prt.pRotate(8,'counterwise')
Prt.pRotate(2,'clockwise')
Prt.pRotate(5,'counterwise')
Prt.pRotate(4,'clockwise')

Grd = Grid(d,N)
Grd.setProt(Prt)
Grd.updateGrid()

print(Grd.getGrid())

print(Prt.posArray)
print("test")