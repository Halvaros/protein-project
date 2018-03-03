from numpy import *
import pickle
from protein import Protein
from grid import Grid
import matplotlib.pyplot as plt


def showIt(Grid):
    Grid.update()
    Grid()


def dFunk(s, T):
    return dMax * exp(-s * T)


def legalTwist(protein):  # legal i forstanden gyldig struktur,
    twist = protein.tryRotate()  # ikke nødvendigvis lavere energi
    while not twist:
        twist = protein.tryRotate()


def plotIt(X, Y, ax=None, T=None):
    if T == None:
        plt.figure()
        plt.plot(X, Y)
        plt.grid()
    else:
        ax.set_title(str(T) + ' K')
        ax.set_ylim(auto=True)  # ([-3.47e-21,-10.4e-21])
        ax.plot(X, Y)


def markerPlot(xMatrix, yMatrix):
    xMarkers = xMatrix[:, -1]
    yMarkers = yMatrix[:, -1]
    otherXMarkers=xMatrix[:,0]
    otherYmarkers=yMatrix[:,0]
    xMatrix = xMatrix.flatten()
    yMatrix = yMatrix.flatten()
    plt.figure()
    plt.plot(xMatrix, yMatrix,linewidth=0.2)
    #plt.plot(otherXMarkers,otherYmarkers,marker='o',color='y')
    #plt.plot(xMarkers, yMarkers,marker='o', color='g')
    plt.plot((xMarkers+otherXMarkers)/2,(yMarkers+otherYmarkers)/2,marker='o',color="r")
    plt.grid()


def averageIt(arr, dArr):
    arr = sum(arr, 1)
    arr = divide(arr, dArr, out=zeros_like(arr), where=dArr != 0)
    return arr


#################################################
gridN = 20
N = 15
dMax = 10000

# s=1E-6 #for Nt=15
s = 1E-6
# sArr = array([1E-6])  # linspace(1E-6, 1E-4, 3)

#################################################

# Gradvis avkjøling

d = 600
TArr = arange(1500, -30, -30)
Nt = len(TArr)
epsilon = zeros((Nt, d))
dArr = arange(1, d * Nt + 1).reshape((Nt, d))
L = zeros((Nt, d))
L[1, 1] = N
Prt = Protein(N, gridN, TArr[0])

for i in range(Nt):
    for j in range(d):
        legalTwist(Prt)
        # showIt(Grd)
        epsilon[i, j] = Prt.U
        L[i, j] = Prt.L
        # print(Prt.U)
    Prt.T = TArr[i]

# print(dArr)

meanEpsilon = averageIt(epsilon.copy(), d)
meanL = averageIt(L, d)
markerPlot(dArr, epsilon)
plotIt(TArr, meanEpsilon)
plt.gca().invert_xaxis()
plotIt(TArr, meanL)
plt.gca().invert_xaxis()

#Grd = Grid(d, gridN, Prt)
#while True:
#    showIt(Grd)
plt.show()