from tools_pol import *
from lab2 import *

POLNUM = 2
MOD = 257
MODDIV2 = MOD//2
MODDIV4 = MOD//4

def keygen():
    # A is sampled directly in NTT form.
    A = [[sample_uniform_pol(MOD) for _ in range(POLNUM)] for _ in range(POLNUM)]
    
    s = [ntt(sample_small_pol(2, MOD)) for _ in range(POLNUM)]
    e = [ntt(sample_small_pol(2, MOD)) for _ in range(POLNUM)]

    As = pol_mv_mul(A, s, MOD)
    t = pol_vec_add(As, e, MOD)

    sk = s
    pk = (A, t)
    return (sk, pk)

def enc(pk, m):
    A, t = pk

    r = [ntt(sample_small_pol(2, MOD)) for _ in range(POLNUM)]
    e1 = [ntt(sample_small_pol(2, MOD)) for _ in range(POLNUM)]
    e2 = ntt(sample_small_pol(2, MOD))

    rA = pol_mv_mul(transpose(A), r, MOD)
    u = pol_vec_add(rA, e1, MOD)

    rt = pol_vec_dot(r, t, MOD)
    m_ntt = ntt(pol_mul_const(m, MODDIV2, MOD))
    v = pol_add(pol_add(rt, e2, MOD), m_ntt, MOD)
    return (u, v)

def dec(sk, ct):
    s = sk
    u, v = ct

    res = pol_sub(v, pol_vec_dot(u, s, MOD), MOD)
    res_coef = intt(res)

    m = [0 if abs_mod(coef, MOD) <= MODDIV4 else 1 for coef in res_coef]
    return m

def test_kgen_enc_dec():
    init_roots(MOD)

    sk, pk = keygen()
    ms = [[0, 1]*(DEG//2)]*10

    test_passed = True
    for m in ms:
        ct = enc(pk, m.copy())
        mdec = dec(sk, ct)
        if m != mdec:
            test_passed = False
    return test_passed

if __name__ == "__main__":
    print("Test keygen-enc-dec:", test_kgen_enc_dec())