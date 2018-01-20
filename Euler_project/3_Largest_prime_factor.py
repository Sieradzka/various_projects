# -*- coding: utf-8 -*-
#The prime factors of 13195 are 5, 7, 13 and 29.
#What is the largest prime factor of the number 600851475143 ?
#Answer: 6857

def largest_prime_factor1(number):
    i = 2
    while i * i < number:
        while number % i == 0:
            number = number / i
        i = i + 1
    return number

#not efficient way of finding primery numbers
def find_prime1(number):
    prime = [2, 3]
    for i in range (3, number + 1, 2):
        for j in range (0, len(prime)): 
            if not i % prime[j]: #odd number modulo prime number            
                break 
            if (j == len(prime)-1):
                prime.append(i) 
    return prime
    
def find_prime1a(number):
    prime = [2, 3]
    for i in range (3, number + 1, 2):
        for j in range (3, i + 1, 2): 
            if not i % j:
                break 
            if (j == i-2):
                prime.append(i) 
    return prime
    
#a better way of finding primery numbers  
def find_prime2(number):
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
    
def largest_prime_factor2(number):
    for i in find_prime2(number):
        if not number % i:
            number = number / i
        if number == 1:
            return i 
            
def largest_prime_factor3(number):
    largest_prime = 1
    prime = find_prime2(number)
    for i in range(1, len(prime)):
        if not number % prime[-i]:
            largest_prime = prime[-i]
            break
    return largest_prime   

#------------------------------------------------------------------------------
#import time
#start_time = time.time()   
import timeit
start = timeit.default_timer()
#------------------------------------------------------------------------------
#The statements here
#print(find_prime2(1000))
print(largest_prime_factor2(13195))
#------------------------------------------------------------------------------
#print("--- %s seconds ---" % (time.time() - start_time))
stop = timeit.default_timer()
print(stop - start)