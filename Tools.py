from math import log2
from random import randrange
from collections import Counter
from itertools import count

def primegen():
    for p in [2,3,5,7]: yield p                 # base wheel primes
    gaps1 = [ 2,4,2,4,6,2,6,4,2,4,6,6,2,6,4,2,6,4,6,8,4,2,4,2,4,8 ]
    gaps = gaps1 + [ 6,4,6,2,4,6,2,6,6,4,2,4,6,2,6,4,2,4,2,10,2,10 ] # wheel2357
    def wheel_prime_pairs():
        yield (11,0); bps = wheel_prime_pairs() # additional primes supply
        p, pi = next(bps); q = p * p            # adv to get 11 sqr'd is 121 as next square to put
        sieve = {}; n = 13; ni = 1              #   into sieve dict; init cndidate, wheel ndx
        while True:
            if n not in sieve:                  # is not a multiple of previously recorded primes
                if n < q: yield (n, ni)         # n is prime with wheel modulo index
                else:
                    npi = pi + 1                # advance wheel index
                    if npi > 47: npi = 0
                    sieve[q + p * gaps[pi]] = (p, npi) # n == p * p: put next cull position on wheel
                    p, pi = next(bps); q = p * p  # advance next prime and prime square to put
            else:
                s, si = sieve.pop(n)
                nxt = n + s * gaps[si]          # move current cull position up the wheel
                si = si + 1                     # advance wheel index
                if si > 47: si = 0
                while nxt in sieve:             # ensure each entry is unique by wheel
                    nxt += s * gaps[si]
                    si = si + 1                 # advance wheel index
                    if si > 47: si = 0
                sieve[nxt] = (s, si)            # next non-marked multiple of a prime
            nni = ni + 1                        # advance wheel index
            if nni > 47: nni = 0
            n += gaps[ni]; ni = nni             # advance on the wheel
    for p, pi in wheel_prime_pairs(): yield p   # strip out indexes

def primes235(limit):
    yield 2; yield 3; yield 5
    if (limit < 7): return
    modPrms = [7,11,13,17,19,23,29,31]
    gaps = [4,2,4,2,4,6,2,6,4,2,4,2,4,6,2,6] # 2 loops for overflow
    ndxs = [0,0,0,0,1,1,2,2,2,2,3,3,4,4,4,4,5,5,5,5,5,5,6,6,7,7,7,7,7,7]
    lmtbf = (limit + 23) // 30 * 8 - 1 # integral number of wheels rounded up
    lmtsqrt = isqrt(limit) - 7
    lmtsqrt = lmtsqrt // 30 * 8 + ndxs[lmtsqrt % 30] # round down on the wheel
    buf = [True] * (lmtbf + 1)
    for i in range(lmtsqrt + 1):
        if (buf[i]):
            ci = i & 7; p = 30 * (i >> 3) + modPrms[ci]
            s = p * p - 7; p8 = p << 3
            for ci in range(ci,ci+8):
                c = s // 30 * 8 + ndxs[s % 30]
                buf[c::p8] = [False] * ((lmtbf - c) // p8 + 1)
                s += p * gaps[ci]
    for i in range(lmtbf - 6 + (ndxs[(limit - 7) % 30])): # adjust for extras
        if (buf[i]): yield (30 * (i >> 3) + modPrms[i & 7])

def pow_log(x:int, b:int) -> int:
    ans = b
    while(x > ans): ans *= b
    return ans//b;

def nextprime(n:int) -> int:
    if n < 2: return 2
    if n == 2: return 3
    n = (n + 1) | 1    # first odd larger than n
    m = n % 6
    if m == 3:
        if is_prime(n+2): return n+2
        n += 4
    elif (m == 5):
        if is_prime(n): return n
        n += 2
    for m in count(n, 6):
        if is_prime(m): return m
        if is_prime(m+4): return m+4

def isqrt(n:int) -> int:
    c = (n << 2)//3

    d = c.bit_length()

    a = d>>1
    if d&1:
        x = 1 << a
        y = (x + (n >> a)) >> 1
    else:
        x = (3 << a) >> 2
        y = (x + (c >> a)) >> 1

    if x != y:
        x = y
        y = (x + n//x) >> 1
        while y < x:
            x = y
            y = (x + n//x) >> 1
    return x

def _check(a:int, s:int, d:int, n:int) -> bool:
    x = pow(a,d,n)
    if (x == 1): return True;
    for i in range(s-1):
        if (x == n - 1):
            return True;
        x = x*x%n
    return x == n-1;

def is_prime(n:int, k:int = 5) -> bool:
    if (n == 2 or n == 3 or n == 5 or n == 7 or n == 11 or n == 13): return True;
    
    if (n < 2 or not n&1): return False;
    
    if (k is None): k = int(log2(n)) + 1
    
    s = 0
    d = n - 1
    while (not d&1):
        d >>= 1
        s += 1
    for i in range(k):
        a = randrange(2,n-1)
        if (not _check(a,s,d,n)):
            return False;
    return True;

def slow_fact(n:int, ans:Counter) -> list:
    while (not(n&1)):
        ans[2] += 1
        n >>= 1
        
    while (not n%3):
        ans[3] += 1
        n //= 3
        
    d = 5
    while (d <= n//d):
        while (not n%d):
            ans[d] += 1
            n //= d
        while (not n%(d+2)):
            ans[d+2] += 1
            n //= (d+2)
        d += 6
    if (n != 1): ans[n] += 1
