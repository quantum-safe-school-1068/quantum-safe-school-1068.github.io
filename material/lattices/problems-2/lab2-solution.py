from tools import *

DEG = 32
DEPTH = 5

# The following variables are set when calling init_roots(mod)
MOD = None
RTREE = None 
INVTREE = None
RLIST = None
VAN = None
IVAN = None

def pos_sqrt_mod(a, mod):
    """Returns the positive square root of 'a' modulo 'mod'."""
    return next(b for b in range(mod) if pow(b, 2, mod) == a % mod)

def inv_mod(a, mod):
    """Returns the inverse of 'a' modulo 'mod'."""
    return pow(a, -1, mod)

def get_roots_tree(depth):
    """Returns the tree of roots as a 2D array."""
    assert(depth > 0)

    rtree = [[pos_sqrt_mod(-1, MOD)]]
    for i in range(1, depth):
        temp = []
        for j in range(2**(i-1)):
            r = rtree[i-1][j]
            temp += [pos_sqrt_mod(r, MOD), pos_sqrt_mod(MOD-r, MOD)]
        rtree.append(temp)

    return rtree

def get_roots_inv_tree(rtree):
    """Returns the tree of inverse-roots as a 2D array."""
    return [[inv_mod(r, MOD) for r in rs] for rs in rtree]

def get_roots_list(rtree):
    """Returns the list of roots."""
    rlist = []
    for r in rtree[len(rtree)-1]:
        rlist.append(r)
        rlist.append(MOD-r)
    return rlist

def get_vander(rlist):
    """Returns the Vandermonde matrix of roots."""
    van = [[1 for _ in range(DEG)] for _ in range(DEG)]
    for i in range(DEG):
        for j in range(1, DEG):
            van[i][j] = (van[i][j-1]*rlist[i]) % MOD
    return van

def get_inv_vander(van):
    """Returns the inverse Vandermonde matrix of roots."""
    ivan = transpose(van)
    ivan[0] = [-x for x in ivan[0]]
    for i in range(1, DEG//2):
        ivan[i], ivan[DEG-i] = ivan[DEG-i], ivan[i]
    inv = pow(-DEG, MOD-2, MOD)
    ivan = [[(inv * x) % MOD for x in xs] for xs in ivan]
    return ivan

def init_roots(mod):
    """Initializes all the global variables for the roots."""
    assert(mod > 0 and mod % 64 == 1) # and isprime(mod)
    global MOD, RTREE, INVTREE, RLIST, VAN, IVAN

    MOD = mod
    RTREE = get_roots_tree(DEPTH)
    INVTREE = get_roots_inv_tree(RTREE)
    RLIST = get_roots_list(RTREE)
    VAN = get_vander(RLIST)
    IVAN = get_inv_vander(VAN)


def ntt(a):
    """Transforms a polynomial in coefficient form to CRT form."""
    polnum = 1 # number of polynomials in a level
    coefnum = len(a) # number of coefs per polynomial in a level

    for i in range(DEPTH):
        for j in range(polnum):
            polstart = coefnum*j
            for k in range(coefnum//2):
                temp1 = a[polstart + k]
                temp2 = a[polstart + k + coefnum//2] * RTREE[i][j]
                a[polstart + k] = (temp1 + temp2) % MOD
                a[polstart + k + coefnum//2] = (temp1 - temp2) % MOD
        
        polnum *= 2 
        coefnum //= 2

    return a

def intt(a):
    """Transforms a polynomial in CRT form to coefficient form."""
    polnum = len(a) # number of polynomials in a level
    coefnum = 1 # number of coefs per polynomial in a level

    for i in reversed(range(DEPTH)):
        for j in range(polnum//2):
            polstart = coefnum*j*2
            for k in range(coefnum):
                temp1 = (a[polstart + k] + a[polstart + k + coefnum])
                temp2 = (a[polstart + k] - a[polstart + k + coefnum])*INVTREE[i][j]
                a[polstart + k] = temp1 % MOD
                a[polstart + k + coefnum] = temp2 % MOD
        
        polnum //= 2 
        coefnum *= 2

    factor = inv_mod(DEG, MOD)
    a = [(factor * coef) % MOD for coef in a]
    return a

def pol_mul_point(a, b, mod):
    """Point-wise product of polynomials in CRT form: a ⊙ b  % mod."""
    assert(len(a) == len(b))

    return [(x * y) % mod for (x, y) in zip(a, b)]

def pol_mul_vander(a, b):
    """Multiplication of polynomials using Vandermonde matrices."""
    a_crt = mv_mul_mod(VAN, a, MOD)
    b_crt = mv_mul_mod(VAN, b, MOD)
    ab_crt = pol_mul_point(a_crt, b_crt, MOD)
    ab = mv_mul_mod(IVAN, ab_crt, MOD)
    return ab

def pol_mul_ntt(a, b):
    """Multiplication of polynomials using NTTs."""
    a_ntt = ntt(a)
    b_ntt = ntt(b)
    ab_ntt = pol_mul_point(a_ntt, b_ntt, MOD)
    ab = intt(ab_ntt)
    return ab

def test_ntt():
    pols = [[sample_uniform(MOD) for _ in range(DEG)] for _ in range(100)]
    
    test_passed = True
    for pol in pols:
        pol_van = mv_mul_mod(VAN, pol, MOD)
        pol_ntt = ntt(pol)
        if pol_ntt != pol_van:
            test_passed = False
    return test_passed

def test_ntt_intt():
    pols = [[sample_uniform(MOD) for _ in range(DEG)] for _ in range(100)]
    
    test_passed = True
    for pol in pols:
        pol_ntt_intt = intt(ntt(pol.copy()))
        if pol_ntt_intt != pol:
            test_passed = False
    return test_passed

def test_mul_ntt():
    apols = [[sample_uniform(MOD) for _ in range(DEG)] for _ in range(100)]
    bpols = [[sample_uniform(MOD) for _ in range(DEG)] for _ in range(100)]
    
    test_passed = True
    for (a, b) in zip(apols, bpols):
        ab_van = pol_mul_vander(a, b)
        ab_ntt = pol_mul_ntt(a, b)
        if ab_ntt != ab_van:
            test_passed = False
    return test_passed


if __name__ == "__main__":
    init_roots(257)

    print("Test NTT:", test_ntt())

    print("Test NTT-INTT:", test_ntt_intt())

    print("Test Mul-NTT:", test_mul_ntt())