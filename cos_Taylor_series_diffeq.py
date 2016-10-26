# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 13:35:24 2016

@author: Laila
"""

import numpy as np
import matplotlib.pyplot as plt


#a)

'''
a[j] = -x**2/((2*j)*(2-1))*a[j-1]
s[j] = s[j-1] + a[-1]
'''

def cos_Taylor(x,n):
    a = np.zeros(n+2)
    s = np.zeros(n+2)
    
    a[0] = 1
    s[0] = 1

    for j in range(1,n+2):
        a[j] = -x**2/((2*j)*(2*j-1))*a[j-1]
        s[j] = s[j-1] + a[j]

    return s[n+1], abs(a[n+1])

print "b:"
n = 10
result = cos_Taylor(np.pi/2,n)
print "With polynomial order %g, cos(pi/2) = %g, error = %g." %(n,result[0], result[1])

x = np.linspace(0,6.28,1000)
y  = np.zeros(1000)
for i in range(1000):
    y[i] =  cos_Taylor(x[i],20)[0]
plt.plot(x,y)
plt.show()


#c)
"""
n = 2
S = x - (x**2)/2! + (x**4)/4!
"""
from math import factorial

def test_cos_Taylor():
    n = 2
    x = 3*np.pi/2
    tol = 1e-6
    
    expected = x - (x**2)/factorial(2.) + (x**4)/factorial(4.)
    computed = cos_Taylor(x,n)[0]
    success = abs(expected-computed) < tol

    assert success

print "c:"
test_cos_Taylor()

#d)

print "d:"
x_values = [0.25*np.pi,0.5*np.pi,np.pi]
n_values = [1, 5,10]

print "%10s %10s %10s %10s" %('x-value', 'order', 'approx', 'exact')
for x in x_values:
    for n in n_values:
        print "%10f  %10d %10f %10f" %(x,n,cos_Taylor(x,n)[0], np.cos(x))
