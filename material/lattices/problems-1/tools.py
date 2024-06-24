import hashlib
from math import ceil, log2

GLOBAL_COUNTER = 0 # counter used as seed for the randomness generation.

def randombytes(num: int) -> bytes:
    """Returns 'num' uniformly random bytes."""
    assert(num >= 0)

    global GLOBAL_COUNTER
    seed = GLOBAL_COUNTER
    blist = []
    for i in range(0, 8):
        blist.append(seed % 256)
        seed = seed // 256
    GLOBAL_COUNTER += 1
    return hashlib.shake_128(bytes(blist)).digest(num)

def byte_to_bits(a: int) -> list[int]:
    """Decomposes an integer in the range [0,255] into a list of 8 bits."""
    assert(0 <= a and a <= 255)

    return [int(x) for x in format(a, '08b')]

def sample_uniform(mod: int) -> int:
    """Returns a uniformly random integer in the range [0, mod)."""
    bits_count = ceil(log2(mod))
    bytes_count = ceil(bits_count/8)
    valid = False
    while not valid:
        res_bytes = randombytes(bytes_count)
        res = int.from_bytes(res_bytes, "big")
        res = res % (2**bits_count)
        if res < mod:
            valid = True
    return res

def sample_small(beta: int, mod: int) -> int:
    """Returns a random integer in the range [-beta, beta].
    
        The returned integer is of the form
            a_{1} + ⋅⋅⋅ + a_{beta} - b_{1} - ⋅⋅⋅ - b_{beta}
        where the a_i's, b_i's are independent uniformly random bits.

        The result is then reduced modulo mod.

        The argument beta has to satisfy 1 <= beta <= 4.
    """
    assert(1 <= beta and beta <= 4)

    random_byte = randombytes(1)[0]
    random_bits = byte_to_bits(random_byte)
    return (sum(random_bits[:beta]) - sum(random_bits[-beta:])) % mod

def abs_mod(a: int, mod: int) -> int:
    """Returns the absolute value of 'a' modulo 'mod'.
    
        The absolute value of 'a' modulo 'mod' is the absolute value of
        the integer congruent to 'a' in the range [-mod/2, mod/2).
    """
    assert(0 <= a and a < mod)

    return a if a <= mod//2 else mod-a

def vec_add_mod(a: list[int], b: list[int], mod: int) -> list[int]:
    """Modular addition of vectors (a + b) % mod"""
    assert(len(a) == len(b))

    return [(x + y) % mod for (x, y) in zip(a, b)]

def vec_sub_mod(a: list[int], b: list[int], mod: int) -> list[int]:
    """Modular subtraction of vectors (a - b) % mod"""
    assert(len(a) == len(b))

    return [(x - y) % mod for (x, y) in zip(a, b)]

def vec_dot_mod(a: list[int], b: list[int], mod: int) -> int:
    """Modular dot product of vectors <a, b> % mod"""
    assert(len(a) == len(b))

    return sum([x * y for (x, y) in zip(a, b)]) % mod

def vec_mul_mod(a: list[int], b: list[int], mod: int) -> list[int]:
    """Modular element-wise product of vectors a ⊙ b  % mod"""
    assert(len(a) == len(b))

    return [(x * y) % mod for (x, y) in zip(a, b)]

def mv_mul_mod(A: list[list[int]], b: list[int], mod: int) -> int:
    """Modular matrix-vector multiplication (A*b) % mod"""
    assert(len(A[0]) == len(b))

    return [vec_dot_mod(row, b, mod) for row in A]

def transpose(A):
    """Returns the transpose of a square matrix"""
    n = len(A)
    return [[A[i][j] for i in range(n)] for j in range(n)]