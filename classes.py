# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 08:39:17 2018

@author: HOS
"""
from numpy import *
#from random import randint
import matplotlib.pyplot as plt

T = 1
T_end = 5000
N_T = 100

class Protein:
    kb = 1.38064852e-23

    # methods
    def __init__(self, p_length=0, g_size=0):
        self.d = int(p_length) #15
        self.N = int(g_size)   #18
        self.node = self.d // 2 #7
        self.posArray = array([[self.N / 2, i] for i in range((self.N - self.d) // 2, (self.N + self.d) // 2)],
                              dtype=int).T

        # array([9,i] for i in range (1,15)])
        #[[ 9  9  9  9  9  9  9  9  9  9  9  9  9  9  9]
        # [ 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15]]

        self.U = 0 #Energi
        self.T = T #Temp
        self.U_ij = random.random_sample((self.d, self.d)) * (-3.47e-21 + 10.4e-21) - 10.4e-21
        #plt.imshow(self.U_ij)
        #plt.colorbar()
        #plt.show()
        # trenger ikke være symm i vårt tilfelle

    def pRotate(self, pivot, direction):  #brekker ved pivot (8), node (7) er fast
        posCopy = self.posArray.copy()
        if direction == 'counterwise':
            A = array([[0, -1], [1, 0]])
        elif direction == 'clockwise':
            A = array([[0, 1], [-1, 0]])
        # deler opp posArray i tre deler, prepivot[0], pivot[1] og postpivot[2]
        temp = hsplit(posCopy, array([pivot - 1, pivot]))

        #[[9 9 9 9 9]  [[9]  [[ 9  9  9  9  9  9  9  9  9]
        #[1 2 3 4 5]],  [6]], [ 7  8  9 10 11 12 13 14 15]]

        if pivot > self.node:  # pivotmonomer større enn noden
            temp[2] -= temp[1]

            #[[9  9  9  9  9  9  9  9  9]  -  [[9   =   [0  ... 0]    9 relativvektorer [[y pekende til høyre
             #[7  8  9 10 11 12 13 14 15]]     [8]]     [1,2,..9]                         x]]

            temp[2] = A @ temp[2]   #2x2 x 2x9 =2x9
            #[[-1 -2 -3 -4 -5 -6 -7]
            #[ 0  0  0  0  0  0  0]]  A gir 9 relativvektorer pekende oppover (y går nedover)

            temp[2] += temp[1]
            print(temp[2])
        #[[-1 - 2 - 3 - 4 - 5 - 6 - 7] + [[9   =[8,7,6..2
        # [0  0  0  0  0  0  0]]           6]]  [8 8 ..8]

        elif pivot < self.node:  # pivotmonomer mindre enn noden
            temp[0] -= temp[1]
            temp[0] = A @ temp[0]
            temp[0] += temp[1]
        return hstack([temp[0], temp[1], temp[2]])  #Stack arrays in sequence horizontally (column wise)

    def legalTwist(self, newArray):
        temp=set(tuple(row) for row in newArray.T)  #kolonner i newArray
        return len(temp)==newArray.shape[1]

        #gir de unike kolonnevektorene (som set), ved overlapp vil newArray ha to eller flere like kolonnevektorer slik at == gir false

    def getU(self, newArray):
        U = 0
        #For hver i- monomer, sjekk de resterende j-monomerene
        for i in range(newArray.T.shape[0] - 2):  # utelater de siste to monomerene som er sjekket ved tidligere iterasjoner
            for j in range(i + 2, newArray.T.shape[0]):  # utelater de to første; seg selv og nærmeste nabo (sterk binding)
                if linalg.norm(newArray.T[j] - newArray.T[i], 1) == 1:  #avstand mellom i-monomer og j-monomer
                    U += self.U_ij[i, j]  # legger U(r_i - r_j), der j>i (upper triangular)
        return U

    def checkU(self, newEnergy):
        if self.U > newEnergy:
            return True
        if self.kb * self.T == 0:
            return False
        else:
            return random.uniform() < exp(-(newEnergy - self.U) / (self.kb * self.T))

    def tryRotate(self):
        piv = random.randint(1, self.d + 1)
        rot = dire(random.randint(0, 2))
        tryArray = self.pRotate(piv, rot)
        if self.legalTwist(tryArray):
            tryEnergy = self.getU(tryArray)
            if self.checkU(tryEnergy):
                self.posArray = tryArray
                self.U = tryEnergy

    def getPos(self):
        return self.posArray


# utility
def dire(a):
    if a == 1:
        return 'counterwise'
    return 'clockwise'


class Grid:
    # constructor
    def __init__(self, p_length, g_size, prot1):
        self.N = g_size
        self.grid = zeros((self.N, self.N), dtype=int)
        self.prot = prot1

    # methods
    def setProt(self, prot1):
        self.prot = prot1
        # if size too large, give error

    def update(self):   #Legg til initial protein
        self.grid[self.grid > 0] = 0
        pos = self.prot.posArray                        #[[ 9  9  9  9  9  9  9  9  9  9  9  9  9  9  9]  rad     (y)
        print(pos)                                      # [ 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15]] kolonne (x)
        for i in range(len(pos[0])):
            self.grid[pos[0, i], pos[1, i]] = i + 1     #i=0: grid[pos[0,0],pos[1,0]]=grid[9,1]=0 blir til 0+1=1
                                                        #i=1: grid[pos[0,1],pos[1,1]]=grid[9,2]=0 blir til 1+1=2
    def getGrid(self):                                  #setter hvert punkt i grid-proteinet til indeksverdi+1
        return self.grid



# utility
def d(a):
    if a == 1:
        return 'counterwise'
    return 'clockwise'



N = 18
d = 15
Prt = Protein(d,N)
Grid = Grid(d,N,Prt)
#Grid.update()
#print(Grid.getGrid())
#print(Prt.posArray)


Prt.pRotate(8,'counterwise')
"""
Prt.pRotate(8,'counterwise')
Prt.pRotate(2,'clockwise')
Prt.pRotate(5,'counterwise')
#2 bonds
Prt.pRotate(9,'counterwise')
#doom
Prt.pRotate(4,'clockwise')
"""
#current workbench


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
#plt.show()

#plott
arr = Grid.getGrid()
plt.plot()
X, Y = arr.nonzero()  #Gir to matriser som sammen gir indeksene til non-zero elements
                      #X=[4 4 4 5 5 6 6 6] Y=[4 5 6 4 6 4 5 6] (origo oppe til venstre)

#Any or all of X, Y, s, and c may be masked arrays, in which case all masks will be combined and only unmasked points will be plotted.

plt.scatter(Y, X, c=arr[X, Y], s=100, marker="o", cmap=plt.cm.viridis) #Marker size is scaled by s and marker color is mapped to c.
plt.axis([0, arr.shape[1], 0, arr.shape[0]])
plt.gca().invert_yaxis()
#plt.grid()
#plt.rcParams['axes.axisbelow'] = True
#plt.show()
"""