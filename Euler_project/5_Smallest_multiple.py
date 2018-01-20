# -*- coding: utf-8 -*-
#2520 is the smallest number that can be divided by each of the numbers from 1 to 10 without any remainder.
#What is the smallest positive number that is evenly divisible by all of the numbers from 1 to 20?
#Answer: 232792560

def gcd(a,b):
    return b and gcd(b, a % b) or a
def lcm(a,b): 
    return a * b / gcd(a,b)
    
def Smallest_multiple ():   
    n = 1
    for i in range(1, 21):
        n = lcm(n, i) 
    print(n)
     
 #------------------------------------------------------------------------------
#import time
#start_time = time.time()   
import timeit
start = timeit.default_timer()
#------------------------------------------------------------------------------
#The statements here
Smallest_multiple()
#print(Smallest_multiple3())
#------------------------------------------------------------------------------
#print("--- %s seconds ---" % (time.time() - start_time))
stop = timeit.default_timer()
print(stop - start)               