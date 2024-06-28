from tools_pol import *
from lab2 import *

POLNUM = 1
MOD = 32257
COEFF_S = 1
COEFF_Y = 2048
MAX_SC = DEG*COEFF_S
Z_BOUND = COEFF_Y-MAX_SC

def hash(w, m: str):
    """Returns the hash of the polynomial (vector) w and the message m.
    
        The hash consists of DEG bits.
    """
    h_bytes = hashlib.shake_128(bytes(str(w) + m, 'utf-8')).digest(DEG//8)
    h_bits = []
    for x in h_bytes:
        h_bits += byte_to_bits(x)
    return h_bits

def sample_mask():
    """Returns a uniformly random mask for the signing procedure.

        The mask is a polynomial with DEG coefficients 
        in [-COEFF_Y, COEFF_Y]. The coefficients are then reduced 
        modulo MOD. The polynomial is returned in CRT form.
    """
    return ntt([(sample_uniform(2*COEFF_Y+1)-COEFF_Y) % MOD for _ in range(DEG)])             

def keygen():
    # A is sampled directly in NTT form.
    A = [[sample_uniform_pol(MOD) for _ in range(POLNUM)] for _ in range(POLNUM)]

    s1 = [ntt(sample_small_pol(COEFF_S, MOD)) for _ in range(POLNUM)]
    s2 = [ntt(sample_small_pol(COEFF_S, MOD)) for _ in range(POLNUM)]

    As1 = pol_mv_mul(A, s1, MOD)
    t = pol_vec_add(As1, s2, MOD)

    sk = (A, s1, s2)
    vk = (A, t)
    return (sk, vk)

def sign(sk, m):
    A, s1, s2 = sk
    
    valid = False
    while not valid:
        y1 = [sample_mask() for _ in range(POLNUM)]
        y2 = [sample_mask() for _ in range(POLNUM)]
        Ay1 = pol_mv_mul(A, y1, MOD)
        w = pol_vec_add(Ay1, y2, MOD)

        h = hash(w, m)
        c = ntt(h.copy())

        cs1 = [pol_mul_point(c, pol, MOD) for pol in s1]
        cs2 = [pol_mul_point(c, pol, MOD) for pol in s2]
        z1 = pol_vec_add(cs1, y1, MOD)
        z2 = pol_vec_add(cs2, y2, MOD)
        z1 = [intt(x) for x in z1]
        z2 = [intt(x) for x in z2]

        valid = True
        for pol in z1:
            for coef in pol:
                if abs_mod(coef, MOD) > Z_BOUND:
                    valid = False
        for pol in z2:
           for coef in pol:
                if abs_mod(coef, MOD) > Z_BOUND:
                    valid = False
    sig = (h, z1, z2)
    return sig

def verify(vk, m, sig):
    A, t = vk
    h, z1, z2 = sig
    c = ntt(h.copy())

    valid = True
    for pol in z1:
        for coef in pol:
            if abs_mod(coef, MOD) > Z_BOUND:
                valid = False
    for pol in z2:
        for coef in pol:
            if abs_mod(coef, MOD) > Z_BOUND:
                valid = False

    z1 = [ntt(x) for x in z1]
    z2 = [ntt(x) for x in z2]
    Az1 = pol_mv_mul(A, z1, MOD)
    ct = [pol_mul_point(c, pol, MOD) for pol in t]
    v = pol_vec_sub(pol_vec_add(Az1, z2, MOD), ct, MOD)

    valid = (valid and hash(v, m) == h)
    return valid

def test_keygen_sign_verify():
    init_roots(MOD)

    sk, vk = keygen()
    ms = [str(i) for i in range(100)]

    test_passed = True
    for m in ms:
        h, z1, z2 = sign(sk, m)
        if not verify(vk, m, sig=(h.copy(), z1.copy(), z2.copy())):
            test_passed = False
        if verify(vk, "-1", sig=(h.copy(), z1.copy(), z2.copy())):
            test_passed = False
            
    return test_passed

if __name__ == "__main__":
    init_roots(MOD)

    print("Test keygen-sign-verify:", test_keygen_sign_verify())