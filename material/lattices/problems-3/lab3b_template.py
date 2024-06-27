from tools_pol import *
from lab2 import *

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
    raise NotImplementedError

def sign(sk, m):
    raise NotImplementedError

def verify(vk, m, sig):
    raise NotImplementedError

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