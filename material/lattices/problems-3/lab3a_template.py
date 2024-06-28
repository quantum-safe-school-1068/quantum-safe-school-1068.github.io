from tools_pol import *
from lab2 import *

POLNUM = 2
MOD = 257
MODDIV2 = MOD//2
MODDIV4 = MOD//4

def keygen():
    raise NotImplementedError

def enc(pk, m):
    raise NotImplementedError

def dec(sk, ct):
    raise NotImplementedError

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
