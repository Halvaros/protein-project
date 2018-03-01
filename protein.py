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
# from random import randint
import matplotlib.pyplot as plt
from grid import Grid
N_T = 100

class Protein:
    kb = 1.38064852e-23

    # methods
    def __init__(self, p_length=0, g_size=0,Temp=1E-2):
        self.d = int(p_length)
        self.N = int(g_size)
        self.node = self.d // 2
        self.posArray = array([[self.N / 2, i] for i in range((self.N - self.d) // 2, (self.N + self.d) // 2)],
                              dtype=int).T
        self.U=float(0)
        self.T = Temp
        self.U_ij = random.random_sample((self.d, self.d)) * (-3.47e-21 + 10.4e-21) - 10.4e-21
        self.twists = 0
        # trenger ikke være symm i vårt tilfelle

    def pRotate(self, pivot, direction):
        if direction == 'counterwise':
            A = array([[0, -1], [1, 0]])
        elif direction == 'clockwise':
            A = array([[0, 1], [-1, 0]])
        # deler opp posArray i tre deler, prepivot[0], pivot[1] og postpivot[2]
        temp = hsplit(self.posArray.copy(), array([pivot - 1, pivot]))
        if pivot > self.node:  # pivotmonomer større enn noden
            temp[2] -= temp[1]
            temp[2] = A @ temp[2]
            temp[2] += temp[1]
        elif pivot < self.node:  # pivotmonomer mindre enn noden
            temp[0] -= temp[1]
            temp[0] = A @ temp[0]
            temp[0] += temp[1]
        return hstack([temp[0], temp[1], temp[2]])

    def legalTwist(self, newArray):
        temp = set(tuple(row) for row in newArray.T)  # kolonner i newArray
        return len(temp) == newArray.shape[1]

    def getU(self, newArray):
        U = 0
        for i in range(newArray.T.shape[0] - 2):  # utelater de siste to monomerene
            for j in range(i + 2, newArray.T.shape[0]):  # utelater de to første; seg selv og nærmeste nabo
                if linalg.norm(newArray.T[j] - newArray.T[i], 1) == 1:
                    U += self.U_ij[i, j]  # legger U(r_i - r_j), der j>i (upper triangular)
        return U

    def checkU(self, newEnergy):
        if self.U > newEnergy:
            return True
        else:
            return random.uniform() < exp(-(newEnergy - self.U) / (self.kb * self.T))

    def tryRotate(self):
        piv = random.randint(1, self.d + 1)
        rot = dire(random.randint(0, 2))
        tryArray = self.pRotate(piv, rot)
        if self.legalTwist(tryArray):
            #print("kjører")
            tryEnergy = self.getU(tryArray)
            if self.checkU(tryEnergy):
                self.posArray = tryArray
                self.U+=tryEnergy
                self.twists += 1
                return True

        return False

    def getPos(self):
        return self.posArray.copy()


# utility
def dire(a):
    if a == 1:
        return 'counterwise'
    return 'clockwise'
