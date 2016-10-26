# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 19:16:36 2016

@author: Laila
"""
from matplotlib.pylab import*

def fortune(F,p,q,I,N):
    index_set = range(N+1) 
    x = zeros(len(index_set))
    c = zeros(len(index_set))
    x[0] = F
    c[0] = (p*q)/(1.0e4)*F
    for n in index_set[1:]:
        x[n] = x[n-1] + (p/100.)*x[n-1] - c[n-1]
        c[n] = c[n-1] + (I/100.)*c[n-1]
    return x,c

F=100 #initial fortune
p=5; q=3; I=4 
N=100 #years

x1, c1 = fortune(F,p,q,I,N)

n = linspace(0,N,N+1)
plot(n,x1)
xlabel('n-years')
ylabel('fortune')
show()
      
'''
(C:\Users\Daniel\Anaconda2) C:\Users\Daniel\Desktop\Laila\INF1100\Exercis
es A>python fortune_and_inflation1.py

'''