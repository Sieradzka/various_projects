# -*- coding: utf-8 -*-
#The sum of the primes below 10 is 2 + 3 + 5 + 7 = 17.
#Find the sum of all the primes below two million.
#Answer = 142913828922
	
def find_prime(number):
    prime = [2]
    not_prime = [] 
    for i in range (3, number + 1, 2):
        prime.append(i)
        for j in range (3, i + 1, 2):
            ij = i*j
            if ij > number:
                break
            not_prime.append(ij)
    prime[:] =  set(prime) - set(not_prime)
    prime.sort()
    return prime
    
def Summation_of_primes(number):
    result = 0
    for i in find_prime(number):
        result += i
    return result
    
#------------------------------------------------------------------------------
#import time
#start_time = time.time()   
import timeit
start = timeit.default_timer()
#------------------------------------------------------------------------------
#The statements here
print (sum(find_prime(2000000)))
#or
print(Summation_of_primes(2000000))
#------------------------------------------------------------------------------
#print("--- %s seconds ---" % (time.time() - start_time))
stop = timeit.default_timer()
print(stop - start)
