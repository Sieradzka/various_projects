# -*- coding: utf-8 -*-
#A Pythagorean triplet is a set of three natural numbers, a < b < c, for which,
#a^2 + b^2 = c^2
#For example, 3^2 + 4^2 = 9 + 16 = 25 = 5^2.
#There exists exactly one Pythagorean triplet for which a + b + c = 1000.
#Find the product abc.
"""
def Special_Pythagorean_triplet():
#    c = 1000 - (a + b)
#    c**2 = a**2 + b**2
#    a*b*c = 10**3*c*( a + b - 50)
    for a in range (10, 999):
        for b in range (a + 1, 999):
            for c in range (b + 1, 999):
                    L = a*b*c
                    R = 10**3*c*( a + b - 50)
                    if (L == R):
                        break
                        return L
"""
for a in range(1, 1000):
    for b in range(a, 1000):
        c = 1000 - a - b
        if c > 0:
            if c*c == a*a + b*b:
                print(a*b*c)
                break
                 
#------------------------------------------------------------------------------
#import time
#start_time = time.time()   
import timeit
start = timeit.default_timer()
#------------------------------------------------------------------------------
#The statements here
#print(Special_Pythagorean_triplet())
#------------------------------------------------------------------------------
#print("--- %s seconds ---" % (time.time() - start_time))
stop = timeit.default_timer()
print(stop - start)