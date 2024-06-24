from tools import *

DIM = 64
MOD = 257

def keygen():
    raise NotImplementedError

def enc(pk, m):
    raise NotImplementedError

def dec(sk, ct):
    raise NotImplementedError

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