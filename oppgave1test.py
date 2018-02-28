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

"""
#folding 
for i in range(100):
    Prt.tryRotate()
    Grid.update()
    plt.imshow(Grid.getGrid());
    #plt.colorbar()
    plt.show()
"""

#easy plot
plt.imshow(Grid.getGrid());
plt.colorbar()
plt.show()

#plott
arr = Grid.getGrid()
plt.plot()
X, Y = arr.nonzero()  #Gir to matriser som sammen gir indeksene til non-zero elements
                      #X=[4 4 4 5 5 6 6 6] Y=[4 5 6 4 6 4 5 6] (origo oppe til venstre)

#Any or all of X, Y, s, and c may be masked arrays, in which case all masks will be combined and only unmasked points will be plotted.

plt.scatter(Y, X, c=arr[X, Y], s=100, marker="o", cmap=plt.cm.viridis) #Marker size is scaled by s and marker color is mapped to c.
plt.axis([0, arr.shape[1], 0, arr.shape[0]])
plt.gca().invert_yaxis()
print(arr[X,Y])
print(X,Y)
#plt.grid()
#plt.rcParams['axes.axisbelow'] = True
plt.show()

