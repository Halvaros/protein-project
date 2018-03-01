# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 16:57:16 2018

@author: HOS
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 08:39:17 2018

@author: HOS
"""
from numpy import *
import matplotlib.pyplot as plt
from grid import Grid

T = 0
T_end = 5000
N_T = 100


class Protein:
    kb = 1.38064852e-23
    #methods
    def __init__(self,p_length = 0, g_size = 0):
        self.d= int(p_length)                           #antall monomerer
        self.N= int(g_size)                             #størrelse grid
        self.node = self.d//2                                #noden
        self.posArray = array([[self.N/2,i] for i in range((self.N-self.d)//2, (self.N+self.d)//2)],dtype = int).T
        self.U = 0
        self.L = int(p_length)                      #lengden til et bestemt protein
        self.T = T
        self.U_ij = random.random_sample((self.d,self.d))*(-3.47e-21 + 10.4e-21) - 10.4e-21
        #trenger ikke være symm i vårt tilfelle
        
    def pRotate(self,pivot, direction):
        if direction == 'counterwise':
            A = array([[0, -1], [1, 0]])
        elif direction == 'clockwise':
            A = array([[0,1],[-1,0]])
        #deler opp posArray i tre deler, prepivot[0], pivot[1] og postpivot[2]
        temp = hsplit(self.posArray.copy(),array([pivot-1,pivot]))
        if pivot > self.node: #pivotmonomer større enn noden
            temp[2] -= temp[1]
            temp[2] = A@temp[2]
            temp[2] += temp[1]
        elif pivot < self.node: #pivotmonomer mindre enn noden
            temp[0] -= temp[1]
            temp[0] = A@temp[0]
            temp[0] += temp[1]
        return hstack([temp[0], temp[1], temp[2]])
        
    def legalTwist(self,newArray):
        return unique(newArray, axis = 1).shape == newArray.shape
    
    def getUL(self, newArray):
        U = 0
        L = 0
        for i in range(newArray.T.shape[0]-2):          #utelater de siste to monomerene
            for j in range(i+2,newArray.T.shape[0]):    #utelater de to første; seg selv og nærmeste nabo
                relDist = linalg.norm(newArray.T[j] - newArray.T[i],2)
                if relDist == 1:
                    U += self.U_ij[i,j]                 #legger U(r_i - r_j), der j>i (upper triangular)
                if relDist > L:                           
                    L = relDist                         #finner største relative avstand mellom to monomerer
        return U, L
    """
    def getU(self,newArray):
        U = 0
        for i in range(newArray.T.shape[0]-2): #utelater de siste to monomerene
            for j in range(i+2,newArray.T.shape[0]): #utelater de to første; seg selv og nærmeste nabo
                if linalg.norm(newArray.T[j] - newArray.T[i],1) == 1:
                    U += self.U_ij[i,j] #legger U(r_i - r_j), der j>i (upper triangular)
        return U
    """
    
    def checkU(self,newEnergy):
        if self.U > newEnergy:
            return True
        return random.uniform() < exp(-(newEnergy - self.U)/(self.kb*self.T)) #antar kb*T > 0
    
    def tryRotate(self):
        piv = random.randint(2,self.d) #utelukker endemonomerer
        rot= dire(random.randint(0,2))
        newArray = self.pRotate(piv,rot)
        if self.legalTwist(newArray):
            newEnergy, newLength = self.getUL(newArray)
            if self.checkU(newEnergy):
                self.posArray = newArray
                self.U = newEnergy
                self.L = newLength
                 
    def getPos(self):
        return self.posArray.copy()
    
#utility
def dire(a):
    if a == 1:
        return 'counterwise'
    return 'clockwise'
