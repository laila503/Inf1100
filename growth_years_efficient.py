# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 16:25:48 2016

@author: Laila
"""

x0 = 100; p = 5; N = 4

outfile = open('growth_eff','w')
outfile.write('year |  money  ' + '\n')
outfile.write('-----|---------' + '\n')

x_old = x0
n=0
while n <= N:
    x_new = x_old + (p/100.)*x_old
    outfile.write('  ' +str(n) + '  |  ' + str(x_old) + '  '+'\n')
    x_old = x_new
    n=n+1
    
outfile.close()  

#A short file to see that the .dat-file made has espected values:
infile = open('growth_eff', 'r')
infile.readline()
infile.readline()
for line in infile:
    word = line.split()
    a = float(word[0])
    b = float(word[-1])
    print a,b
  
'''
(C:\Users\Daniel\Anaconda2) C:\Users\Daniel\Desktop\Laila\INF1100\Exercises 
A>python growth_years_efficient.py
0.0 100.0
1.0 105.0
2.0 110.25
3.0 115.7625
4.0 121.550625
'''