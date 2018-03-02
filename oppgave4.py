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

        # plt.show()


def averageIt(arr, dArr):
    arr = sum(arr, 1)
    arr = divide(arr, dArr, out=zeros_like(arr), where=dArr != 0)
    return arr


#################################################
gridN = 20
N = 15
dMax = 10000
Nt = 20
# s=1E-6 #for Nt=15
s = 1E-6
# sArr = array([1E-6])  # linspace(1E-6, 1E-4, 3)
TArr = linspace(1E-2, 1500, Nt)
#################################################

# Midlere E  og  L som funksjon av T
# Gradient E som funksjon av T

epsilon = zeros((Nt, dMax))
L = zeros((Nt, dMax))
dArr = zeros(Nt)
L[1, 1] = N
for i in range(Nt):
    d = int(round(dFunk(s, TArr[i])))
    dArr[i] = d
    print(d)
    Prt = Protein(N, gridN, TArr[i])
    # Grd = Grid(d, gridN, Prt)
    for j in range(d):
        legalTwist(Prt)
        # showIt(Grd)
        epsilon[i, j] = Prt.U
        L[i, j] = Prt.L
        # print(Prt.U)

#print(dArr)
dEpsilon = sum(epsilon, 1).copy()
dEpsilon = diff(dEpsilon) / diff(TArr)
epsilon = averageIt(epsilon, dArr)
L = averageIt(L, dArr)

# epsilon = divide(epsilon, dArr, out=zeros_like(epsilon), where=dArr != 0)

plotIt(TArr, epsilon)
plotIt(TArr, L)
TArr = linspace(1E-2, 1500, Nt - 1)
plotIt(TArr, dEpsilon)
plt.show()

# Bindingsenergi for T=0K og T=500K
'''

d = 5000
dArr = arange(1, d + 1)
TArr = array([1E-2, 500]).astype('f8')

fig, axArr = plt.subplots(1, 2, sharex=True)
for i in range(len(TArr)):
    epsilon = zeros_like(dArr).astype('f8')
    Prt = Protein(N, gridN, TArr[i])
    for j in range(d):
        legalTwist(Prt)
        epsilon[j]=Prt.U
        #showIt(Grd)
        #print(epsilon[j])

    plotIt(dArr, epsilon, axArr[i], TArr[i])
plt.show()

for i in range(len(TArr)):
    epsilon = zeros_like(dArr).astype('f8')
    Prt = Protein(N, gridN, TArr[i])
    for j in range(d):
        legalTwist(Prt)
        epsilon[j]=Prt.U
   
    plt.figure(i)
    plt.plot(dArr,epsilon)
plt.show()
'''

# Midlere E ved T=0
'''
reps=20
epsilon = zeros((reps, dMax))
T=1E-2
d = int(round(dFunk(s, T)))
for i in range(reps):
    Prt = Protein(N, gridN, T)
    #Grd = Grid(d, gridN, Prt)
    for j in range(d):
        legalTwist(Prt)
        # showIt(Grd)
        epsilon[i, j] = Prt.U
        # print(Prt.U)

epsilon = sum(epsilon, 1)
epsilon =epsilon/d
plotIt(arange(1,reps+1), epsilon)
plt.show()
'''
