# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 17:00:23 2018

@author: HOS
"""
from numpy import *
from protein import Protein
from grid import Grid
import matplotlib.pyplot as plt

N = 18
d = 15
Prt = Protein(d,N)
Grid = Grid(d,N,Prt)

#folding 
for i in range(30):
    Prt.tryRotate()
    Grid.update()
    #Grid.easyPlot()
    Grid.ShowOff()

