from numpy import *
from protein import Protein
from grid import Grid
import matplotlib.pyplot as plt


def showIt(Grid):

    Grid.update()
    Grid()


def dFunk(s, T):
    return dMax * exp(-s * T)


gridN = 20
N = 15
dMax = 10000
Nt = 10

sArr = linspace(1E-4, 1E-2, 3)
TArr = linspace(1, 1500, Nt)
sDic = {}

for s in sArr:
    epsilon = zeros((Nt, dMax))
    L = zeros((Nt, dMax))
    dArr = zeros(Nt)
    L[1, 1] = N
    for i in range(Nt):
        d = int(round(dFunk(s, TArr[i])))
        dArr[i] = d
        print(d)
        Prt = Protein(N, gridN, TArr[i])
        Grid = Grid(d, gridN, Prt)
        for j in range(d):
            twist = Prt.tryRotate()
            showIt(Grid)
            while not twist:
                twist = Prt.tryRotate()
            epsilon[i, j] = Prt.U
            #print(Prt.U)

    epsilon = sum(epsilon, 1)
    epsilon = divide(epsilon, dArr, out=zeros_like(epsilon), where=dArr != 0)
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

