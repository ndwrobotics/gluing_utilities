from sage.all_cmdline import * 
sys.path.append(os.environ['HOME'] + "/lmfdb")
from lmfdb import db
from compute_traces import is_frob_matrix_diagonal_batch, find_identity_prime

#OLD
def search_mf_by_traces(trace_function, conductor, prime_bound, modulus, level_bound, high_weight=False):
    """Returns an iterator of all modular forms with dimension 1, weight 2,
    level <= level_bound, and a_p \\equiv trace_dict[a_p] \mod modulus for all
    prime p not dividing conductor and less than prime_bound, except those for which p is a prime factor of the 
    level of the modular form."""
    query = {'level': {'$lte': level_bound}, 'dim': 1}
    if high_weight:
        query['weight'] = modulus + 1
    else: query['weight'] = 2
    sub_query = []
    for p in prime_range(prime_bound):
        if(conductor % p == 0):
            continue
        #query['traces.' + str(index)] = {'$mod': [value, modulus]}
        sub_query.append({'$or': [
            {'traces.' + str(p) : {'$mod': [trace_function(p), modulus]}},
            {'level_primes': {'$contains': p}}
        ]})
    query['$and'] = sub_query
    return db.mf_newforms.search(query, projection=['label'])
    
def search_ec_by_traces(trace_function, conductor, prime_bound, modulus):
    """Returns an iterator of all ellilptic curves with 
    a_p \\equiv trace_dict[a_p] \mod modulus for all
    prime p not dividing conductor and less than prime_bound, except those for which p is a prime factor of the 
    level of the modular form."""
    query = []
    for i, p in enumerate(prime_range(prime_bound)):
        if(conductor % p == 0):
            continue
        try:
            trace = trace_function(p)
        except:
            continue
        sub_query = [{'aplist.' + str(i+1) : {'$mod' : [trace, modulus]}}]
        if(modulus == 3 or modulus == p):
            #additive reduction
            sub_query.append({'conductor' : {'$mod' : [0, p**2]}})
        x = (p+1) % modulus
        if(trace == x):
            #split multiplicative reduction
            sub_query.append({'$and' : 
                [{'conductor' : {'$mod' : [0,p]}}, {'aplist.' + str(i+1) : 1}]})
        if(trace == (-x) % modulus):
            #non-split multiplicative reduction
            sub_query.append({'$and' : 
                [{'conductor' : {'$mod' : [0,p]}}, {'aplist.' + str(i+1) : -1}]})
        query.append(sub_query[0] if len(sub_query)==1 else {'$or': sub_query})
    curve_labels = [x['lmfdb_iso'] for x in db.ec_classdata.search({'$and' : query}, projection=['lmfdb_iso'])]
    return [EllipticCurve(ec['ainvs'])
        for ec in db.ec_curvedata.search({'lmfdb_iso': {'$or': curve_labels}}, projection=['ainvs'])
    ]

def filter_reducible_case(C, N, l, elliptic_curves, bound=500):
    if len(elliptic_curves) == 0:
        return elliptic_curves
    if elliptic_curves[0].isogenies_prime_degree(l) == []:
        return elliptic_curves
    #Step 1: filter by discriminant (i.e., ramification of Gal rep)
    filtered = []
    for E in elliptic_curves:
        D = E.discriminant()
        good = true
        for p,e in factor(D):
            if e % l != 0 and (not E.has_additive_reduction(p)) and N % p != 0:
                good = false
        if good:
            filtered.append(E)
    elliptic_curves = filtered
    #Step 2: filter by diagonal matrix
    #generate a list of interesting primes
    E = elliptic_curves[0]
    primes = []
    for p in prime_range(bound):
        if (E.ap(p)**2 - 4*p) % l == 0:
            primes.append(p)
    print(primes)
    is_diagonal_C = is_frob_matrix_diagonal_batch(C, l, primes)
    print(is_diagonal_C)
    result = []
    for E in elliptic_curves:
        is_diagonal_E = is_frob_matrix_diagonal_batch(E, l, primes)
        good = True
        for p,q in zip(is_diagonal_C, is_diagonal_E):
            if(p == 1 and q == 0):
                good = False
                break
        if(good):
            result.append(E)
    if len(result) > 100:
        result = filter_large_prime(C, l, result, bound=40000)
    return result
    
def filter_large_prime(C, l, elliptic_curves, bound=40000):
    result = []
    primes = find_identity_prime(C,l,bound)
    for E in elliptic_curves:
        ok = True
        for p in primes:
            if E.discriminant() % p == 0:
                continue
            f = E.change_ring(GF(p)).division_polynomial(l)
            x = f.variables()[0]
            if power_mod(x, p , f) != x:
                ok = False
                break
        if ok:
            result.append(E)
    return result
