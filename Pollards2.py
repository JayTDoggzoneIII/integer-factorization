from math import gcd
from random import randrange, seed as rseed
from Tools import pow_log, isqrt, primegen
from time import perf_counter
from sys import stdout
'''
import ctypes
kernel32 = ctypes.windll.kernel32
kernel32.SetConsoleMode(kernel32.GetStdHandle(-11), 7)'''

rseed(perf_counter())

def Brent(n:int) -> int:
    if (not n&1): return 2
    y,c,m = randrange(1,n), randrange(1,n), randrange(1,n)
    g = r = q = 1
    count = 0
    while (g == 1 and count < 17):
        x = y
        for i in range(r):
            y = (y*y%n + c)%n
        k = 0
        while (k < r and g == 1):
            ys = y
            for i in range(m if (m < r-k) else r-k):
                y = (y*y%n + c)%n
                q = q*abs(x-y)%n
            g = gcd(q,n)
            k += m
        r <<= 1
        count += 1
    if (g == 1):
        while (count < 10_000):
            ys = (ys*ys%n+c)%n
            g = gcd(abs(x-ys),n)
            if (g > 1): return g;
            count += 1
    return g;

def rho_Brent(n:int, verbose:bool = False) -> int:
    if (verbose):
        stdout.write("Starting the Rho Brent method\n")
        stdout.flush()
    start = perf_counter()
    y,c,m = randrange(1,n), randrange(1,n), randrange(1,n)
    g = r = q = 1
    while (g == 1 and r < 600_000):
        x = y
        for i in range(r):
            y = (y*y%n + c)%n
        k = 0
        while (k < r and g == 1):
            ys = y
            for i in range(m if (m < r-k) else r-k):
                y = (y*y%n + c)%n
                q = q*abs(x-y)%n
            g = gcd(q,n)
            k += m
        r <<= 1
    if (g == n):
        ys = (ys*ys%n+c)%n
        g = gcd(abs(x-ys),n)        
        while (g == 1):
            ys = (ys*ys%n+c)%n
            g = gcd(abs(x-ys),n)
    if (verbose):
        if (1 < g < n):
            print('\033[92m' + "Found factor: %d"%g + '\033[0m')
            stdout.write("Time elapsed: %.3fsec\n\n"%(perf_counter() - start))
            stdout.flush()
        else:
            print('\033[31m' + "Failed with Rho Brent, trying Pollard's p - 1" + '\033[0m')
            stdout.write("Time elapsed: %.3fsec\n\n"%(perf_counter() - start))
            stdout.flush()               
    return g;

def P1(n:int, verbose:bool = False, B1:int = 100, B2:int = 300) -> int:
    
    if (verbose):
        stdout.write("Starting the Pollard's p - 1 method\n")
        stdout.flush()    
    start = perf_counter()
    while (B2 <= 3_000_000):
        it = primegen()
        q = 2
        p = next(it)
        while (p <= B1): q,p = pow(q,pow_log(B1,p),n), it.next_prime()
        g = gcd(q-1,n)
        if (1 < g < n): 
            if (verbose):
                print('\033[92m' + "Found factor: %d"%g + '\033[0m')
                stdout.write("Time elapsed: %.3fsec\n\n"%(perf_counter() - start))
                stdout.flush()
                
            return g
        while (p <= B2): q,p = pow(q,p,n), it.next_prime()
        g = gcd(q-1,n)
        if (1 < g < n): 
            if (verbose):
                print('\033[92m' + "Found factor: %d"%g + '\033[0m')
                stdout.write("Time elapsed: %.3fsec\n\n"%(perf_counter() - start))
                stdout.flush()

            return g
        B1 = B2
        B2 *= 10
    if (verbose):
        print('\033[31m'+"Failed with Pollard's p - 1, trying ECM"+'\033[0m')
        stdout.write("Time elapsed: %.3fsec\n\n"%(perf_counter() - start))
        stdout.flush()
    return 1;

def Floyd(n:int, c:int = 1) -> int:
    x,y,g,k = randrange(n), randrange(n), 1, 0
    while (g == 1 and k < 100_000):
        x = (x*x%n + c)%n;
        y = (y*y%n + c)%n;
        y = (y*y%n + c)%n;
        g = gcd(abs(x-y),n)
        k += 1
    return g;

def rho(n:int, verbose:bool = False) -> int:
    root = isqrt(n)
    if (n == root * root): return root;
    if (verbose):
        stdout.write("Starting the Modern Rho method\n")
        stdout.flush()
    start = perf_counter()    
    x = randrange(1,n)
    y = 1
    i = 0
    stage = 2
    while (gcd(n, abs(x - y)) == 1 and i < 360_000):
        if (i == stage):
            y = x
            stage <<= 1
        x = (x*x - 1)%n
        i += 1
    g = gcd(n,abs(x - y))
    if (verbose):
        if (1 < g < n):
            print('\033[92m' + "Found factor: %d"%g + '\033[0m')
            stdout.write("Time elapsed: %.3fsec\n\n"%(perf_counter() - start))
            stdout.flush()
        else:
            print('\033[31m' + "Failed with Modern Rho, trying Rho Brent" + '\033[0m')
            stdout.write("Time elapsed: %.3fsec\n\n"%(perf_counter() - start))
            stdout.flush()     
    return g;