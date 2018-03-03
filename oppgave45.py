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


def plotIt(X, Y, ax=None, N=None):
    if ax == None:
        plt.figure()
        plt.plot(X, Y)
        plt.grid()
    else:
        ax.set_title(str(N))
        ax.grid()
        ax.invert_xaxis()
        ax.plot(X, Y)


def markerPlot(xMatrix, yMatrix, ax):
    xMarkers = xMatrix[:, -1]
    yMarkers = yMatrix[:, -1]
    otherXMarkers = xMatrix[:, 0]
    otherYmarkers = yMatrix[:, 0]
    xMatrix = xMatrix.flatten()
    yMatrix = yMatrix.flatten()
    ax.plot(xMatrix, yMatrix, linewidth=0.2)
    # ax.plot(otherXMarkers,otherYmarkers,marker='o',color='y')
    # ax.plot(xMarkers, yMarkers,marker='o', color='g')
    ax.plot((xMarkers + otherXMarkers) / 2, (yMarkers + otherYmarkers) / 2, marker='o', color="r")


def averageIt(arr, dArr):
    arr = sum(arr, 1)
    arr = divide(arr, dArr, out=zeros_like(arr), where=dArr != 0)
    return arr


def chillax(N, gridN):
    epsilon = zeros((Nt, d))
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

    meanEpsilon = averageIt(epsilon.copy(), d)
    meanL = averageIt(L, d)
    return epsilon, meanEpsilon, meanL


NArr = array([15, 30])
gridNArr = array([20, 35])
# s=1E-6 #for Nt=15
s = 1E-6
d = 600
TArr = arange(1500, -30, -30)
Nt = len(TArr)
dArr = arange(1, d * Nt + 1).reshape((Nt, d))

fig1, axArr1 = plt.subplots(1, 2, sharex=True)
fig2, axArr2 = plt.subplots(1, 2, sharex=True)
fig3, axArr3 = plt.subplots(1, 2, sharex=True)
for i in range(2):
    epsilon, meanEpsilon, meanL = chillax(NArr[i], gridNArr[i])
    markerPlot(dArr, epsilon, axArr1[i])
    plotIt(TArr, meanEpsilon, axArr2[i], NArr[i])
    plotIt(TArr, meanL, axArr3[i], NArr[i])

    # Grd = Grid(d, gridN, Prt)
    # while True:
    #    showIt(Grd)

plt.show()
