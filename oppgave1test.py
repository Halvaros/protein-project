# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 17:00:23 2018

@author: HOS
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 16:57:44 2018

@author: HOS
"""
from numpy import *
from protein import Protein
from grid import Grid
import matplotlib.pyplot as plt

N = 16
d = 10

def showIt(Prt,Grid):
    Grid = Grid(d,N,Prt)
    Grid.update()
    Grid()

Prt = Protein(d,N)
showIt(Prt,Grid)
while Prt.twists<20:
    if Prt.tryRotate():
        showIt(Prt,Grid)


Prt.pRotate(6,'counterwise')
Prt.pRotate(8,'counterwise')
Prt.pRotate(2,'clockwise')
Prt.pRotate(5,'counterwise')
#2 bonds
Prt.pRotate(9,'counterwise')
#doom
Prt.pRotate(4,'clockwise')

#current workbench

"""
#folding 
for i in range(100):
    Prt.tryRotate()
    Grid.update()
    plt.imshow(Grid.getGrid());
    #plt.colorbar()
    plt.show()

#easy plot
plt.imshow(Grid.getGrid());
plt.colorbar()
plt.show()
"""




