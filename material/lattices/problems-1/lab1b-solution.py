from numpy import log2
from numpy.polynomial import Polynomial

# Polynomial that represents the probability distribution of the random
# variable s = a0 + a1 - a2 - a3, where the ai's are uniformly random 
# bits sampled independently. It is shifted by x**2, i.e., 
#   Pr[s = -2] = single.coef[0], 
#   Pr[s = -1] = single.coef[1], 
#   Pr[s =  0] = single.coef[2], etc.
#
# single = 1/16 + 1/4 x + 3/8 x**2 + 1/4 x**3 + 1/16 x**4
single = Polynomial((1/16, 1/4, 3/8, 1/4, 1/16))

# Polynomial that represents the probability distribution of the random
# variable z = y*w, where y and w are independently distributed
# according to the distribution of s above. It is shifted by x**4.
#  
# product = 1/128 + 1/16 x**2 + 1/8 x**3 + 39/64 x**4 + 1/8 x**5
#           + 1/16 x**6 + 1/128 x**8
product = Polynomial((1/128, 0, 1/16, 1/8, 39/64, 1/8, 1/16, 0, 1/128))

# Polynomial that represents the probability distribution of the random
# variable e = z1 + ... + z128 + s, where the zi's are independently
# distributed according to the distribution of z above, and s is
# independent of the zi's and distributed according to the distribution
# of the single variable from the start.
# It is shifted by (x**4)**128 * x**2 = x**514
error = (product**64)*(product**64)*single

# Decryption will be correct as long as the error e is in [-63, 63],
# and we need to take into account the x**514 shift.
prob_correct_dec = sum([error.coef[i+514] for i in range(-63, 64)])
prob_dec_fail = 1 - prob_correct_dec

print("Probability of correct decryption =", prob_correct_dec)
print("Probability of decryption failure = {} ≈ 2**({})"
      .format(prob_dec_fail, log2(prob_dec_fail)))
