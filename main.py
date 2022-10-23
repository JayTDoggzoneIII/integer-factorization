from gc import disable
from time import perf_counter
from collections import Counter
from math import factorial
from random import randrange
from Tools import is_prime, primegen 
from Pollards2 import rho_Brent, rho, P1
from mpECM import mp_ECM
from Quadratic_sieve_algorithm2 import MPQS
from colorama import init
init()

from sys import stdin, stdout, setrecursionlimit
#stdin = open("input.txt","r") 
#stdout = open("output.txt","w")
      
#setrecursionlimit((1<<31)-1)

gets = input
puts = print
input = stdin.readline
print = stdout.write

def fact(n:int, ans:Counter, verbose:bool) -> None:
    if (n == 1):
        return;
    elif (is_prime(n)):
        ans[n] += 1
    else:
        if (verbose):
            stdout.write("Factor the number: %d\n\n"%n)
            stdout.flush()
        div = rho(n,verbose)
        if (div == 1): div = rho_Brent(n,verbose)
        if (div == 1): div = P1(n,verbose)
        #if (div == 1): div = ECM(n,verbose) if (len(str(n)) < 55) else mp_ECM(n,verbose)
        if (div == 1): div = mp_ECM(n,verbose)
        if (div == 1): div = MPQS(n,verbose)
        fact(div, ans, verbose)
        fact(n//div, ans, verbose)

def main() -> int:
    #disable()
    n = eval(gets("Input your number (0 to finish):\n"))
    if (n): verbose = True if (gets("Do you want displaying details? (y/n):\n") == 'y') else False
    while (n):
        start = perf_counter()
        if (type(n) != int):
            print("n must be integer\n")
            print("Working time: %.3fsec\n"%(perf_counter()- start))
        elif (n <= 1):
            print("n must be larger then 1\n")
            print("Working time: %.3fsec\n"%(perf_counter() - start))
        elif (is_prime(n)):
            print("Number of digits: %d\n\n"%(len(str(n))))
            print("This number is prime\n")
            print("Working time: %.3fsec\n"%(perf_counter() - start))
        else:
            print("Number of digits: %d\n\n"%(len(str(n))))
            ans = Counter()
            pg = primegen()
            p = next(pg)
            if (verbose):
                stdout.write("Starting the trial division method\n\n")
                stdout.flush()                
            while (p < 65535):
                if (verbose and not n%p):
                    puts('\033[92m' + "Found factor: %d"%p + '\033[0m')
                    stdout.write("Time elapsed: %.3f\n\n"%(perf_counter() - start))
                    stdout.flush()                    
                while (not n%p):
                    ans[p] += 1
                    n //= p
                p = next(pg)
            if (n != 1): fact(n,ans,verbose)
            print("Prime divisors:\n")
            for i in sorted(ans):
                if (ans[i] != 1):
                    print("%i^%i  "%(i,ans[i]))
                else:
                    print("%i  "%i)
            print("\n")
            print('\033[42m'+"Working time: %.3fsec"%(perf_counter() - start) + '\033[0m\n')
        print("\n\n")
        n = eval(gets("Input your number (0 to finish):\n"))
        if (not n): break;
        verbose = True if (gets("Do you want displaying details? (y/n):\n") == 'y') else False
    gets("Press Enter to exit...")
    return 0;

if (__name__ == "__main__"):   
    main()
