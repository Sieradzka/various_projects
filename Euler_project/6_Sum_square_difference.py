# -*- coding: utf-8 -*-
#The sum of the squares of the first ten natural numbers is,
#1^2 + 2^2 + ... + 10^2 = 385
#The square of the sum of the first ten natural numbers is,
#(1 + 2 + ... + 10)2 = 55^2 = 3025
#Hence the difference between the sum of the squares of the first ten natural numbers and the square of the sum is 3025 − 385 = 2640.
#Find the difference between the sum of the squares of the first one hundred natural numbers and the square of the sum.
#Answer: 25164150

def Squre_sum(n):
    return (n*(n + 1))**2 / 4

def Sum_squre(n):
    return n/6*(n + 1)*(2*n + 1)
    
def Sum_square_difference(n):
    return Squre_sum(n) - Sum_squre(n)

#iterative approach:   
def Sum_square_difference1(n):
    Sum_squre, Squre_sum = 0, 0
    for i in range(1, n + 1):
        Sum_squre += i**2
        Squre_sum += i
    Squre_sum = Squre_sum**2
    return Squre_sum - Sum_squre
    
def Sum_square_difference1a(n):    
    r = range(1, n + 1)
    sum_n = sum(r)
    return sum_n ** 2 - sum(i**2 for i in r)
#------------------------------------------------------------------------------
#import time
#start_time = time.time()   
import timeit
start = timeit.default_timer()
#------------------------------------------------------------------------------
#The statements here
print(Sum_square_difference1a(100))
#------------------------------------------------------------------------------
#print("--- %s seconds ---" % (time.time() - start_time))
stop = timeit.default_timer()
print(stop - start)  