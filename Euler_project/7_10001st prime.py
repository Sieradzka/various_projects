# -*- coding: utf-8 -*-
#By listing the first six prime numbers: 2, 3, 5, 7, 11, and 13, we can see that the 6th prime is 13.
#What is the 10 001st prime number?
#Answer = 104743

def prime(number):
    prime = [2, 3]
    i = 3
    while len(prime) < number:
        for j in range (0, len(prime)): 
            if not i % prime[j]: #odd number modulo prime number            
                break 
            if (j == len(prime)-1):
                prime.append(i) 
        i += 2
    return prime
    
#------------------------------------------------------------------------------
#import time
#start_time = time.time()   
import timeit
start = timeit.default_timer()
#------------------------------------------------------------------------------
#The statements here

print_prime = prime(10001)
print(print_prime[-1])
#------------------------------------------------------------------------------
#print("--- %s seconds ---" % (time.time() - start_time))
stop = timeit.default_timer()
print(stop - start)

