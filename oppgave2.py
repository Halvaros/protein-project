from numpy import *
from protein import Protein
from grid import Grid
import matplotlib.pyplot as plt


def showIt(Prt, Grid):
    Grid = Grid(d, N, Prt)
    Grid.update()
    Grid()


def dFunk(s, T):
    return dMax * exp(-s * T)

N = 20
L= 15
dMax = 1e4
Nt = 100
convergenceLimit=200

sArr = linspace(1E-4,1E-2, 3)
TArr = linspace(1E-2, 1500, Nt)
sDic = {}

for s in sArr:
    epsilon=zeros(Nt)
    for i in range(Nt):
        totalE=0
        d = int(round(dFunk(s, TArr[i])))
        #print(d)
        Prt = Protein(L, N, TArr[i])
        for j in range(d):
            twist = Prt.tryRotate()
            if twist:
                totalE+=Prt.U
        if d==0:
            break
        else:
            epsilon[i]=totalE/Prt.twists

    sDic[s] = epsilon
    plt.plot(TArr, epsilon)
    plt.gca().invert_yaxis()
    plt.show()

"""    
showIt(Prt, Grid)
while Prt.twists < 20:
    if Prt.tryRotate():
        showIt(Prt, Grid)
"""