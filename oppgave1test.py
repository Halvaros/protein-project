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

N = 18
d = 15
Prt = Protein(d,N)
Grid = Grid(d,N,Prt)
#print(Prt.posArray)

"""
Prt.pRotate(6,'counterwise')
Prt.pRotate(8,'counterwise')
Prt.pRotate(2,'clockwise')
Prt.pRotate(5,'counterwise')
#2 bonds
Prt.pRotate(9,'counterwise')
#doom
Prt.pRotate(4,'clockwise')
"""
#current workbench


for i in range(100):
    Prt.tryRotate()
    Grid.update()
    plt.imshow(Grid.getGrid());
    #plt.colorbar()
    plt.show()
    
    

print(Prt.U)

print(Prt.posArray)
"""


print(Prt.getU(Prt.posArray))
print(Prt.getU(Prt.posArray))
print(Prt.getU(Prt.posArray))
"""
#print(Prt.legalTwist(Prt.posArray))
