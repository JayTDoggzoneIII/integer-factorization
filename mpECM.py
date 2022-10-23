from sys import stdin, stdout, setrecursionlimit
from gc import disable
from time import perf_counter
from multiprocessing import Queue as mpQueue, Process, Value, cpu_count
from math import gcd, exp, log
from random import randrange, seed as rseed
from Tools import isqrt, primegen

'''
import ctypes
kernel32 = ctypes.windll.kernel32
kernel32.SetConsoleMode(kernel32.GetStdHandle(-11), 7)'''


def add(p1:list, p2:list, p0:list, n:int) -> list:
    x1,z1 = p1; x2,z2 = p2; x0,z0 = p0
    t1, t2 = (x1-z1)*(x2+z2), (x1+z1)*(x2-z2)
    return (z0*(t1+t2)*(t1+t2) % n, x0*(t1-t2)*(t1-t2) % n)
 
def double(p:list, A:list, n:int) -> list: 
    x, z = p; An, Ad = A
    t1, t2 = (x+z)*(x+z)%n, (x-z)*(x-z)%n
    t = t1 - t2
    return (((t1*t2)<<2)*Ad % n, (((Ad*t2)<<2) + t*An)*t % n)
 
def multiply(m:int, p:list, A:list, n:int) -> int: 
    if (m == 0): return (0, 0)
    elif (m == 1): return p
    else:
        q = double(p, A, n)
        if (m == 2): return q
        b = m >> 1
        while (b > (b & -b)): b ^= b & -b
        r = p
        while (b):
            if m&b: q, r = double(q, A, n), add(q, r, p, n)
            else:   q, r = add(r, q, p, n), double(r, A, n)
            b >>= 1
        return r

def pow_log(x:int, b:int) -> int:
    ans = b
    while(x > ans): ans *= b
    return ans//b;    

def try_curve(n:int, B1:int, B2:int) -> int:
    seed = randrange(6,n)
    u, v = (seed*seed - 5) % n, (seed<<2) % n
    p = pow(u, 3, n)
    Q, C = (pow(v-u,3,n)*(3*u+v) % n, ((p*v)<<2) % n), (p, pow(v,3,n))
        
    pg = primegen()
    p = next(pg)
    while (p <= B1): Q, p = multiply(pow_log(B1, p), Q, C, n),pg.next_prime()
    g = gcd(Q[1], n)
    if (1 < g < n): 
        return g
    while (p <= B2):
        Q = multiply(p, Q, C, n)
        g = g*Q[1]%n
        p = next(pg)
    g = gcd(g, n)
    if (1 < g < n): 
        return g
    return 1

def _ECM(n:int, out:mpQueue, bup:mpQueue, verbose:bool = False, B1:int = 1000, B2:int = 3100) -> int:
        
    #if (n > 5e28): return 1;
    iters, size = 1, len(str(n))
    rseed(perf_counter())
    for k in range(3 if (size < 50) else 4):

        for _ in range(iters):
            if (not out.empty()): 
                out.put(1);
                return
            seed = randrange(6,n)
            u, v = (seed*seed - 5) % n, (seed<<2) % n
            p = pow(u, 3, n)
            Q, C = (pow(v-u,3,n)*(3*u+v) % n, ((p*v)<<2) % n), (p, pow(v,3,n))
            if (verbose and out.empty()):
                bup.put(0)
                #stdout.write("Trying curve #%d with %d and %d boundaries:\n    y^2 = x^3 + (%d/%d - 2)x^2 + x (mod %d)\n"%(bup.qsize(),B1,B2,Q[0],Q[1],n))
                print("Trying curve #%d with bounds B1 = %d B2 = %d\r"%(bup.qsize(),B1,B2),end='')
                stdout.flush()
                
            pg = primegen()
            p = next(pg)
            while (p <= B1): Q, p = multiply(pow_log(B1, p), Q, C, n),pg.next_prime()
            g = gcd(Q[1], n)
            if (1 < g < n): 
                if (verbose and out.empty()):
                    print('\n\033[K\033[92m' + "Found factor: %d"%g + '\033[0m')
                out.put(g)
                return
            while (p <= B2):
                Q = multiply(p, Q, C, n)
                g = g*Q[1]%n
                p = next(pg)
            g = gcd(g, n)
            if (1 < g < n): 
                if (verbose and out.empty()):
                    print('\n\033[K\033[92m' + "Found factor: %d"%g + '\033[0m')             
                out.put(g)
                return

        B1 = B2
        B2 <<= 2
        iters += k+2
    
    return


def _ECM2(n:int, out:mpQueue, bup:mpQueue, verbose:bool = False, B1:int = 1000, B2:int = 3000) -> int:
        
    #if (n > 5e28): return 1;
    iters, size = 1, len(str(n))
    rseed(perf_counter())
    for k in range(7 if (size < 50) else 8):

        for _ in range(2+k if (k < 5) else 6):
            if (not out.empty()): 
                out.put(1);
                return
            seed = randrange(6,n)
            u, v = (seed*seed - 5) % n, (seed<<2) % n
            p = pow(u, 3, n)
            Q, C = (pow(v-u,3,n)*(3*u+v) % n, ((p*v)<<2) % n), (p, pow(v,3,n))
            if (verbose and out.empty()):
                bup.put(0)
                #stdout.write("Trying curve #%d with %d and %d boundaries:\n    y^2 = x^3 + (%d/%d - 2)x^2 + x (mod %d)\n"%(bup.qsize(),B1,B2,Q[0],Q[1],n))
                print("Trying curve #%d with bounds B1 = %d B2 = %d\r"%(bup.qsize(),B1,B2),end='')
                stdout.flush()
                
            pg = primegen()
            p = next(pg)
            while (p <= B1): Q, p = multiply(pow_log(B1, p), Q, C, n),pg.next_prime()
            g = gcd(Q[1], n)
            if (1 < g < n): 
                if (verbose and out.empty()):
                    print('\n\033[K\033[92m' + "Found factor: %d"%g + '\033[0m')
                out.put(g)
                return
            while (p <= B2):
                Q = multiply(p, Q, C, n)
                g = g*Q[1]%n
                p = next(pg)
            g = gcd(g, n)
            if (1 < g < n): 
                if (verbose and out.empty()):
                    print('\n\033[K\033[92m' + "Found factor: %d"%g + '\033[0m')             
                out.put(g)
                return

        B1 <<= 1
        B2 <<= 1
        #iters += k+2
    
    return

def get_ret(f, n, out, bup, verbose):
    ans = f(n,out, bup, verbose)
    if (ans != 1):
        out.put(ans)

def mp_ECM(n,verbose):
    if (verbose):
        stdout.write("Starting the Lenstra elliptic-curve method\n")
        stdout.flush()
    start = perf_counter()
    
    out, bup = mpQueue(), mpQueue()
    procs = []
    #i = 0
    for p in range(cpu_count()):
        p = Process(target = _ECM2, args = (n, out, bup, verbose))
        p.start()
        procs.append(p)
    for p in procs:
        p.join()
            
       
    if (not out.empty()): 
        if (verbose):
            print("\033[KTime elapsed: %.3fsec\n"%(perf_counter() - start))
            stdout.flush()             
        return out.get()
    
    if (verbose):
        print('\n\033[31m'+"Failed with ECM, trying MPQS"+'\033[0m')
        stdout.write("Time elapsed: %.3fsec\n\n"%(perf_counter() - start))
        stdout.flush()
    
    return 1;