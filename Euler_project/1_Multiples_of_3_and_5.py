#Multiples of 3 and 5
#If we list all the natural numbers below 10 that are multiples of 3 or 5,
#we get 3, 5, 6 and 9. The sum of these multiples is 23.
#Find the sum of all the multiples of 3 or 5 below 1000.
#Answer: 233168

#The simplest solution is to go from 1 (or actually from 3) 
#up to n (here 1000) and divide each number by 3. 
#If its modulo equals 0 then the number is added to the sum 
#if not the number is divided by 5; if its modulo equals 0 
#then the number is added to the sum.

def sum_of_multiples1(n):
    result = 0
    for i in range(n):
        if not (i % 3 and i % 5):
            result += i
    return result
 
#Different approched is to start from the higest number m (here 999).
#In the first step the number n is divided by 3.
#The floor of this integer divission, q, is the quotient.
#Calculating the sum of the first q natural numbers: q(q+1)/2
#and multiplying it by 3, gives the sum of all multiples of 3.
#The same step is to be repeated for number 5.
#Next the sum of multiples of 3 and 5 is added.
#The last step is to subtract from the sum these multiples
#which are the same for 3 and 5 (as they were added to the sum twice).
 
def sum_of_multiples2(m):
    n = m - 1 #this line is to be consistent with function: sum_of_multiples1(n)
    result = 0
    for i in (3, 5): #the loop to process step 1 and 2
        q = n // i #Calculating the floor of the argument divided by 3 and 5,
        result += i * q * (q + 1) // 2       
    q = n // 15
    result -= 15 * q * (q + 1) // 2  
    return result 

 #This is more compact version of function sum_of_multiples2(m)

def sum_of_multiples2a(m):
    n = m - 1 #this line is to be consistent with function: sum_of_multiples1(n)
    result = 0
    for i in (3, 5, -3*5): #the loop to process step 1 and 2
        q = n // abs(i) #the absolute value is needed as the floor would be different for negative value
        result += i * q * (q + 1) // 2       
    return result 
    
#------------------------------------------------------------------------------
#import time
#start_time = time.time()   
import timeit
start = timeit.default_timer()
#------------------------------------------------------------------------------
#The statements here
print(sum_of_multiples2(1000))
#------------------------------------------------------------------------------
#print("--- %s seconds ---" % (time.time() - start_time))
stop = timeit.default_timer()
print(stop - start)
