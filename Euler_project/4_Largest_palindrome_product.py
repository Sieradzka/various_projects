# -*- coding: utf-8 -*-
#A palindromic number reads the same both ways. The largest palindrome made from the product of two 2-digit numbers is 9009 = 91 × 99.
#Find the largest palindrome made from the product of two 3-digit numbers.
#Answer: 906609 

def palindrome ():
    largest_palindrome = 0
    for i in range(999, 99, -1):
        for j in range(i, 99, -1):
            multiplication = i*j
            if multiplication > largest_palindrome:
                if str(multiplication) == str(multiplication)[::-1]:
                    largest_palindrome = multiplication
                    print("{} * {} = {}".format(i, j, largest_palindrome))
    return largest_palindrome
                 
                 
#------------------------------------------------------------------------------
#import time
#start_time = time.time()   
import timeit
start = timeit.default_timer()
#------------------------------------------------------------------------------
#The statements here
print(palindrome())
#------------------------------------------------------------------------------
#print("--- %s seconds ---" % (time.time() - start_time))
stop = timeit.default_timer()
print(stop - start)