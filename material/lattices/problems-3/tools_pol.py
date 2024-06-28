from tools import *
from lab2 import DEG, pol_mul_point

def sample_uniform_pol(mod):
    """Returns a uniformly random polynomial in coeff/CRT form."""
    return [sample_uniform(mod) for _ in range(DEG)]

def sample_small_pol(beta, mod):
    """Returns a random polynomial with coefficients in [-beta, beta].
    
        The coefficients are reduced modulo mod.
    """
    return [sample_small(beta, mod) for _ in range(DEG)]

def pol_add(a, b, mod):
    """Polynomial addition for coeff/CRT form."""
    return vec_add_mod(a, b, mod)

def pol_sub(a, b, mod):
    """Polynomial subtraction for coeff/CRT form."""
    return vec_sub_mod(a, b, mod)

def pol_mul_const(a, c, mod):
    """Polynomial multiplication by a constant for coeff/CRT form."""
    return [(c*x) % mod for x in a]

def pol_vec_add(a, b, mod):
    """Polynomial vector addition for coeff/CRT form."""
    return [pol_add(x, y, mod) for (x, y) in zip(a, b)]

def pol_vec_sub(a, b, mod):
    """Polynomial vector subtraction for both coeff and CRT form."""    
    return [pol_sub(x, y, mod) for (x, y) in zip(a, b)]

def pol_vec_dot(a, b, mod):
    """Polynomial vector dot-product in CRT form."""
    vec = [pol_mul_point(x, y, mod) for (x, y) in zip(a, b)]
    return [sum([pol[i] for pol in vec]) % mod for i in range(DEG)]

def pol_mv_mul(A, b, mod):
    """Polynomial matrix-vector multiplication in CRT form."""
    return [pol_vec_dot(row, b, mod) for row in A]
