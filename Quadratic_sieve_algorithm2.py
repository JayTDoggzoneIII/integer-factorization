from math import log, gcd
from bisect import bisect_left
from itertools import count
from random import randrange
import time

from Tools import is_prime, isqrt, nextprime, primegen
'''
import ctypes
kernel32 = ctypes.windll.kernel32
kernel32.SetConsoleMode(kernel32.GetStdHandle(-11), 7)'''

def prod(a:list) -> int:
    [*p] = a
    jump = 1

    while (jump < len(p)):
        for i in range(0, len(p) - jump, jump << 1):
            p[i] *= p[i + jump]
            p[i + jump] = None

        jump <<= 1

    return p[0]

def legendre(a:int, p:int) -> int: 
    return ((pow(a, (p-1) >> 1, p) + 1) % p) - 1

def mod_sqrt(n:int, p:int) -> int:
    a = n%p
    if (p%4 == 3): return pow(a, (p+1) >> 2, p)
    elif (p%8 == 5):
        v = pow(a << 1, (p-5) >> 3, p)
        i = ((a*v*v << 1) % p) - 1
        return (a*v*i)%p
    elif (p%8 == 1): # Shank's method
        q, e = p-1, 0
        while (not q&1):
            e += 1
            q >>= 1
        n = 2
        while (legendre(n, p) != -1): n += 1
        w, x, y, r = pow(a, q, p), pow(a, (q+1) >> 1, p), pow(n, q, p), e
        while (True):
            if (w == 1): return x
            v, k = w, 0
            while (v != 1 and k+1 < r):
                v = (v*v)%p
                k += 1
            if k == 0: return x
            d = pow(y, 1 << (r-k-1), p)
            x, y = (x*d)%p, (d*d)%p
            w, r = (w*y)%p, k
    else: return a # p == 2

def mod_inv(a:int, m:int) -> int:
    a, x, u = a%m, 0, 1
    while (a): x, u, m, a = u, x - (m//a)*u, a, m%a
    return x

# Multiple Polynomial Quadratic Sieve
# assumes n is composite
def MPQS(n:int, verbose:bool = False) -> int:
    if (verbose): print("Starting the MPQS method")
    
    if verbose:
        time1 = time.time()

    root_n = isqrt(n)
    root_2n = isqrt(n<<1)

    # formula chosen by experimentation
    # seems to be close to optimal for n < 10^50
    #bound = int(7.5 * (len(str(n))-1)**2)
    bound = int(7.5 * log(n,10)**2)
    prime, mod_root, log_p, num_prime = [], [], [], 0

    # size of the sieve
    x_max = bound * 5
    x_max_2 = x_max<<1

    # maximum value on the sieved range
    m_val = (x_max * root_2n) >> 1

    # find a number of small primes for which n is a quadratic residue
    pg = primegen()
    p = next(pg)
    leg = n&1
    while p < bound or num_prime < 3:
        
        # legendre (n|p) is only defined for odd p
        if leg == 1:
            prime.append(p)
            log_p.append(log(p,10))
            r = mod_sqrt(n, p)
            roots = [r]
            q = p
            while q < x_max:
                # find all square roots mod p^n via Hensel Lifting
                r = (r + (n - r*r)*mod_inv(r<<1, q))%q
                #assert r*r%q == n%q
                roots.append(r)
                q *= p
            mod_root.append(roots)
            num_prime += 1
        elif not leg:
            if verbose:
                print('trial division found factors:')
                print(p, 'x', n//p)
            return p
        p = next(pg)
        leg = legendre(n%p, p)

        

    # fudging the threshold down a bit makes it easier to find powers of partial-partial
    # relationships, but it also makes the smoothness check slower. reducing by twice the log
    # of the largest prime in the factor base results in cofactors less than that value squared
    thresh = log(m_val,10) - (log_p[-1]*2)

    # skip small primes. they contribute very little to the log sum
    # and add a lot of unnecessary entries to the table
    # instead, fudge the threshold down a bit, according to expected number of factors
    min_prime = int(thresh*2)
    sp_idx = bisect_left(prime, min_prime)
    sieve_primes = prime[sp_idx:]

    fudge = sum(log_p[i]/(prime[i]-1) for i in range(sp_idx))

    sums = [fudge]*x_max_2

    if verbose:
        print('smoothness bound:', bound)
        print('sieve size:', x_max)
        print('log threshold:', thresh)
        print('skipping primes less than:', min_prime)

    smooth = []
    used_prime = set()
    partial = {}
    num_smooth = prev_num_smooth = num_used_prime = num_poly = num_partial = 0
    num_poly = 0
    root_A = isqrt(root_2n // x_max)

    if verbose:
        print('sieving for smooths...')
    prev_p = 0
    while (True):
    # find an integer value A such that:
    # A is =~ sqrt(2*n) / x_max
    # A is a perfect square
    # sqrt(A) is prime, and n is a quadratic residue mod sqrt(A)
        
        if (root_A < 18400_00000_00000_00000): pg.skipto(root_A)
        while (True):
            root_A = pg.next_prime() if (root_A < 18400_00000_00000_00000) else nextprime(root_A)
            leg = legendre(n, root_A)
            if leg == 1:
                break
            elif leg == 0:
                if verbose:
                    print('dumb luck found factors:')
                    print(root_A, 'x', n//root_A)
                return root_A

        A = root_A * root_A

        # solve for an adequate B
        # B*B is a quadratic residue mod n, such that B*B-A*C = n
        # this is unsolvable if n is not a quadratic residue mod sqrt(A)
        b = mod_sqrt(n, root_A)
        B = (b + (n - b*b) * mod_inv(b<<1, root_A))%A

        # B*B-A*C = n <=> C = (B*B-n)/A
        C = (B*B - n) // A

        num_poly += 1

        # sieve for prime factors
        i = sp_idx
        for p in sieve_primes:
            logp = log_p[i]

            e = 0
            q = p
            while q < x_max:
                inv_A = mod_inv(A, q)
                # modular root of the quadratic
                
                mrie = mod_root[i][e]
                a = ((mrie - B) * inv_A)%q
                b = ((q - mrie - B) * inv_A)%q

                amx = a+x_max
                bmx = b+x_max

                #apx = amx-q
                #bpx = bmx-q

                #k = q
                for k in range(q,x_max,q):
                    sums[amx-q+k] += logp
                    sums[bmx-q+k] += logp
                    sums[amx-k] += logp
                    sums[bmx-k] += logp
                    #k += q

                q *= p
                e += 1

            i += 1

        # check for smooths
        x = -x_max
        #i = 0
        
        for i in range(x_max_2):
            v = sums[i]
            if v > thresh:
                vec = set()
                sqr = []
                # because B*B-n = A*C
                # (A*x+B)^2 - n = A*A*x*x+2*A*B*x + B*B - n
                #               = A*(A*x*x+2*B*x+C)
                # gives the congruency
                # (A*x+B)^2 = A*(A*x*x+2*B*x+C) (mod n)
                # because A is chosen to be square, it doesn't need to be sieved
                sieve_val = (A*x + (B<<1))*x + C

                if sieve_val < 0:
                    vec = {-1}
                    sieve_val = -sieve_val

                for p in prime:
                    while not(sieve_val%p):
                        if p in vec:
                            # keep track of perfect square factors
                            # to avoid taking the sqrt of a gigantic number at the end
                            sqr.append(p)
                        vec ^= {p}
                        sieve_val //= p

                if sieve_val == 1:
                    # smooth
                    smooth.append((vec, (sqr, (A*x+B), root_A)))
                    used_prime |= vec
                elif sieve_val in partial:
                    # combine two partials to make a (xor) smooth
                    # that is, every prime factor with an odd power is in our factor base
                    pair_vec, pair_vals = partial[sieve_val]
                    sqr.extend(list(vec & pair_vec))
                    sqr.append(sieve_val)
                    vec ^= pair_vec
                    smooth.append((vec, (sqr + pair_vals[0], (A*x+B)*pair_vals[1], root_A*pair_vals[2])))
                    used_prime |= vec
                    num_partial += 1
                else:
                    # save partial for later pairing
                    partial[sieve_val] = (vec, (sqr, A*x+B, root_A))
            x += 1

            # reset the value for the next go
            sums[i] = fudge
            #i += 1
      
        prev_num_smooth = num_smooth
        num_smooth = len(smooth)
        num_used_prime = len(used_prime)
        if verbose:
            now_p = 100 * num_smooth // num_prime
            if (now_p > prev_p):
                #print("\r", end = "")
                print('\033[34m',min(100,100 * num_smooth // num_prime), '% complete\033[0m\r',sep = '', end="")
                prev_p = now_p

        if num_smooth > num_used_prime and num_smooth > prev_num_smooth:
            if verbose:
                print('%d polynomials sieved (%d values)'%(num_poly, num_poly*x_max_2))
                print('found %d smooths (%d from partials) in %.3f seconds'%(num_smooth, num_partial, time.time()-time1))
                print('solving for non-trivial congruencies...')

            # set up bit fields for gaussian elimination
            masks = []
            mask = 1
            bit_fields = [0]*num_used_prime
            for vec, vals in smooth:
                masks.append(mask)
                #i = 0
                for i,p in enumerate(used_prime):
                    if p in vec: bit_fields[i] |= mask
                    #i += 1
                mask <<= 1

            # row echelon form
            col_offset = 0
            null_cols = []
            for col in range(num_smooth):
                pivot = col-col_offset == num_used_prime or bit_fields[col-col_offset] & masks[col] == 0
                for row in range(col+1-col_offset, num_used_prime):
                    if bit_fields[row] & masks[col]:
                        if pivot:
                            bit_fields[col-col_offset], bit_fields[row] = bit_fields[row], bit_fields[col-col_offset]
                            pivot = False
                        else:
                            bit_fields[row] ^= bit_fields[col-col_offset]
                if pivot:
                    null_cols.append(col)
                    col_offset += 1

            # reduced row echelon form
            for row in range(num_used_prime):
                # lowest set bit
                mask = bit_fields[row] & -bit_fields[row]
                for up_row in range(row):
                    if bit_fields[up_row] & mask:
                        bit_fields[up_row] ^= bit_fields[row]

            # check for non-trivial congruencies
            for col in null_cols:
                all_vec, (lh, rh, rA) = smooth[col]
                lhs = lh   # sieved values (left hand side)
                rhs = [rh] # sieved values - n (right hand side)
                rAs = [rA] # root_As (cofactor of lhs)
                #i = 0
                for i,field in enumerate(bit_fields):
                    if field & masks[col]:
                        vec, (lh, rh, rA) = smooth[i]
                        lhs.extend(list(all_vec & vec))
                        lhs.extend(lh)
                        all_vec ^= vec
                        rhs.append(rh)
                        rAs.append(rA)
                    #i += 1

                factor = gcd(prod(rAs)*prod(lhs) - prod(rhs), n)
                if 1 < factor < n:
                    break
            else:
                if verbose:
                    print('none found.')
                continue
            break

    if verbose:
        print('\033[92m'+ 'factor found: %d'%factor + '\033[0m')
        print('Time elapsed: %.3f seconds'%(time.time()-time1))
    return factor
