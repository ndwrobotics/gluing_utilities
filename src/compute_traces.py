from sage.all_cmdline import * 

SAGE_EXTCODE = SAGE_ENV['SAGE_EXTCODE'] 
script_dir = os.path.dirname(os.path.abspath(__file__))
magma.chdir(script_dir)
magma.attach_spec("../../CHIMP/CHIMP.spec")
magma.attach("helpers.m")


## using Magma implementation, compute frobenius polynomial of curve C at prime p
def frob_poly(C, p):
    C1 = C.change_ring(GF(p))
    return magma.function_call('FrobeniusPoly', [C1]).sage()

## given a genus two curve C, conductor N, and the prime bound
## return a finite superset of gluable primes given a curve and conductor,
## where we filter up to primes of a given bound.
def gluable_primes(C, N, bound=Integer(200)):
    bad_primes = ZZ(N).prime_divisors()
    sqfree_cond = squarefree_part(N)
    d = sqrt(N / sqfree_cond)
    Zd = Integers(d)
    M = Integer(0)
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    last_witness = Integer(0)
    Q = {}
    for p in prime_range(bound):
        if p in bad_primes:
            continue
        try:
            Q[p] = C.change_ring(GF(p)).frobenius_polynomial()
        except:
            continue
        f = Zd(p).multiplicative_order()
        r = Q[p].resultant(x**f-Integer(1))
        if(M == Integer(0) or r % M != Integer(0)):
            last_witness = p
            M = gcd(M, r)
        if(M == Integer(1)):
            break
    gluable_bad_primes = set(bad_primes)
    for l in bad_primes:
        for p in prime_range(bound):
            if p in bad_primes:
                continue
            if p not in Q:
                try:
                    Q[p] = C.change_ring(GF(p)).frobenius_polynomial()
                except:
                    continue
            f = Integers(l*d)(p).multiplicative_order()
            if Q[p].resultant(x**f-Integer(1)) % l != 0:
                gluable_bad_primes.discard(l)
                break
    return set(ZZ(M).prime_divisors()).union(gluable_bad_primes), last_witness

#compute the conductor bound of the one-dimensional representation according to Proposition 3.3 of the paper 
def conductor_bound(N,l):
    N1 = N
    while(N1 % l == Integer(0)):
        N1 = N1 // l
    d = sqrt(N1 / squarefree_part(N1))
    return d*l if N%l == Integer(0) else d

#input: curve C with conductor N and prime l
#output: list of all characters corresponding to H
def all_characters(C, N, l, bound=Integer(300)):
    d = conductor_bound(N,l)
    X = set(DirichletGroup(d, GF(l)).list())
    for p in prime_range(bound):
        if N%p == Integer(0) or p==l:
            continue
        try:
            Q = C.change_ring(GF(p)).frobenius_polynomial().change_ring(GF(l))
        except:
            continue
        for chi in list(X):
            if Q(chi(p)) != Integer(0):
                X.remove(chi)
    return X

#Given a curve C, a prime l, and a character chi
#output a trace function that correspond to chi.
def trace_function_from_character(C, l, chi):
    def trace(p):
        Q = C.change_ring(GF(p)).frobenius_polynomial().change_ring(GF(l))
        ap = -Q.list()[3]
        return (ap - chi(p) - p/chi(p)).lift()
    return trace
    

#input: curve C with conductor N and prime l
#output: list of functions that takes prime p not dividing N to frob trace at p
def frob_trace_functions(C,N,l,bound=Integer(300)):
    X = all_characters(C,N,l,bound)
    result = []
    for chi in X:
        def trace(p):
            Q = C.change_ring(GF(p)).frobenius_polynomial().change_ring(GF(l))
            ap = -Q.list()[3]
            return (ap - chi(p) - p/chi(p)).lift()
        result.append(trace)
    return result
    
## Given a curve C, a prime l, and a prime p,
## check whether the matrix of Frob_p (in C) is diagonal mod l or not
## should work for both genus one and two
def is_frob_matrix_diagonal(C, l, p):
    if isinstance(C, sage.schemes.elliptic_curves.ell_generic.EllipticCurve_generic):
        f = C.change_ring(GF(p)).division_polynomial(l)
        x = f.variables()[0]
        result.append(1  if power_mod(x, p**(l-1) , f) == x  else 0 )
    else:
        return magma.function_call('IsFrobeniusMatrixDiagonal', [C,l,p]).sage()
## Given a curve C, a prime l, and list `primes` of prime numbers,
## check whether the matrix of Frob_p (in C) is diagonal mod l or not for each prime p < bound
## should work for both genus one and two
## return a list of -1 (bad reduction), 0 (not diagonal), 1 (diagonal)
def is_frob_matrix_diagonal_batch(C, l, primes):
    if isinstance(C, sage.schemes.elliptic_curves.ell_generic.EllipticCurve_generic):
        result = []
        for p in primes:
            if C.discriminant() % p == 0 :
                result.append(-1 )
            else:
                f = C.change_ring(GF(p)).division_polynomial(l)
                x = f.variables()[0]
                result.append(1  if power_mod(x, p**(l-1) , f) == x  else 0 )
        return result
    else:
        return magma.function_call('IsFrobeniusMatrixDiagonalBatch', [C,l,primes]).sage()

# Given a genus 2 curve C, a prime l, and a `bound`,
# find all primes up to `bound` for wihch Frob_p acts trivially on l-torsion of Jac(C).
def find_identity_prime(C, l, bound):
    return magma.function_call('FindIdentityPrime', [C,l,bound]).sage()

# Given a genus 2 curve C, its conductor N, and a prime l,
# always return True if the representation corresponding to Hperp/H is reducible.
# returns False if the representation is irreducible, and `bound` is sufficiently high.
def is_HperpH_reducible(C, N, l, bound=300):
    f = frob_trace_functions(C,N,l,bound)[0]
    d = conductor_bound(N,l)
    X = set(DirichletGroup(d, GF(l)).list())
    for p in prime_range(bound):
        if N%p == Integer(0) or p==l:
            continue
        try:
            for chi in list(X):
                if (f(p) % l) != chi(p) + p/chi(p):
                    X.remove(chi)
        except:
            continue
    return (len(X)!=0)
    