# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 07:45:32 2018

@author: HOS
"""
from numpy import *

#1.1 visualize
N = 10
grid = zeros((N,N))
#reset polynomial
start = int(N/2-5)
end = start + 10
print(start,end)
grid[int(N/2),start:end] = arange(1,11)
print(grid)