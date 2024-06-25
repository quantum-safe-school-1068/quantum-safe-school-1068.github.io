from tools import *

DIM = 64
MOD = 257
MODDIV2 = 128
MODDIV4 = 64

def keygen():
    A = [[sample_uniform(MOD) for _ in range(DIM)] for _ in range(DIM)]
    s = [sample_small(2, MOD) for _ in range(DIM)]
    e = [sample_small(2, MOD) for _ in range(DIM)]

    As = mv_mul_mod(A, s, MOD)
    t = vec_add_mod(As, e, MOD)

    sk = s
    pk = (A, t)
    return (sk, pk)

def enc(pk, m):
    A, t = pk

    r = [sample_small(2, MOD) for _ in range(DIM)]
    e1 = [sample_small(2, MOD) for _ in range(DIM)]
    e2 = sample_small(2, MOD)

    rA = mv_mul_mod(transpose(A), r, MOD)
    u = vec_add_mod(rA, e1, MOD)

    rt = vec_dot_mod(r, t, MOD)
    v = (rt + e2 + m*MODDIV2) % MOD
    return (u, v)

def dec(sk, ct):
    s = sk
    u, v = ct

    res = (v - vec_dot_mod(u, s, MOD)) % MOD
    m = 0 if abs_mod(res, MOD) <= MODDIV4 else 1
    return m

def test_kgen_enc_dec():
    sk, pk = keygen()
    ms = [0, 1]*100

    test_passed = True
    for m in ms:
        ct = enc(pk, m)
        mdec = dec(sk, ct)
        if m != mdec:
            test_passed = False
    return test_passed

if __name__ == "__main__":
    print("Test keygen-enc-dec:", test_kgen_enc_dec())