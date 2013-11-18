#!/usr/bin/env python
import binascii
import mathlib
import random
import sys


class AsmuthBloom(object):
    def __init__(self, threshold):
        # threshold is (shares to recombine, all shares)
        self.threshold = threshold
        self.shares = None
        self._m_0 = 0
        self._y = 0
        self._secret = 0
        
    def _find_group_for_secret(self, k):
        """Generate group Z/Zm_0 for secret, where m_0 is prime and m_0 > secret."""
        while True:
            m_0 = mathlib.get_prime(k)
            if (mathlib.primality_test(m_0)):
                return m_0

    def _check_base_condition(self, d):
        """Check if d satisfy the Asmuth-Bloom base condition.

        """
        recomb_count, all_count = self.threshold

        left = 1
        for i in xrange(1, recomb_count + 1):
            left = left * d[i]

        right = d[0]
        for i in xrange(0, recomb_count - 1):
            right = right * d[all_count - i]
        return left > right
    
    def _get_pairwise_primes(self, k, h):
        """Generate d = n+1 primes for Asmuth-Bloom threshold scheme and secret 
        such that d_0 is k-bit prime and d_1 is h-bit prime.
        (d_1...d_n should be pairwise coprimes)
        """
        if (h < k):
            raise Exception('Not enought bits for m_1')
        _, all_count = self.threshold
        # p is picked randomly simple number
        p = self._find_group_for_secret(k)
        while True:
            d = [p]
            # all_count consecutive primes starting from h-bit prime
            for prime in mathlib.get_consecutive_primes(all_count, h):
                d.append(prime)
            if (self._check_base_condition(d)):
                return d
            
    def _prod(self, coprimes):
        """Calculate M=m_1*m_2*...*m_t."""
        M = 1
        t, _ = self.threshold
        for i in xrange(0,t):
            M = M * coprimes[i]
        return M
    
    def _get_modulo_base(self, secret, coprimes):
        """Calculate M' = secret + some_number * taken_prime
        that should be less that coprimes prod.
        Modulos from this number will be used as shares.
        """
        prod = self._prod(coprimes)
        while True:
            A = mathlib.get_random_range(1, (prod - secret) / self._m_0)
            y = secret + A * self._m_0
            if (0 <= y < prod):
                break
        return y
    
    # k is m_0_bits and h is m_1_bits
    def generate_shares(self, secret, k, h):
        if (mathlib.bit_len(secret) > k):
            raise ValueError("Secret is too long")

        m = self._get_pairwise_primes(k, h)
        self._m_0 = m.pop(0)
        
        self._y = self._get_modulo_base(secret, m)
        
        self.shares = []
        for m_i in m:
            self.shares.append((self._y % m_i, m_i))
        # shares item format: (ki, di) ki - mods, di - coprimes
        return self.shares
            
    def combine_shares(self, shares):
        y_i = [x for x, _ in shares] # remainders
        m_i = [x for _, x in shares] # coprimes
        y = mathlib.garner_algorithm(y_i, m_i)
        d = y % self._m_0
        return d


def stringToLong(s):
    return long(binascii.hexlify(s), 16)

if len(sys.argv) < 4:
    print 'Usage: ./bloom.py (--random <bits> | <path>) <M> <N>'
    print ' --random <bits>     - generate random secret'
    print ' <path>              - read secret from file'
    print ' <M>                 - number of shares'
    print ' <N>                 - number of shares needed for recovery'
    sys.exit(1)

source = sys.argv[1]
if source == '--random':
    random = random.SystemRandom()
    secret = random.getrandbits(int(sys.argv.pop(2)))
else:
    try:
        secret = stringToLong(open(source).read())
    except:
        print 'Could not read the source file'
        sys.exit(1)

try:
    m = int(sys.argv[2])
    print 'Got M = %d' % m
except:
    print 'Invalid M'

try:
    n = int(sys.argv[3])
    print 'Got N = %d' % n
except:
    print 'Invalid N'

if n > m:
    print 'N should be less or equal than M'
    sys.exit(1)

threshold = (n, m)
m_0_bits = 500
m_1_bits = 800
    
print '--------------------------------------'
print "Secret: %s" % secret
    
ab = AsmuthBloom(threshold)

try:
    shares = ab.generate_shares(secret, m_0_bits, m_1_bits)
except ValueError, e:
    print 'Cannot generate shares: ' + str(e)
    sys.exit(1)
    
print "Secret shares:"
for i in xrange(0,m):
    print "%s: %s\n" % (i+1, shares[i])
    
print '--------------------------------------'

print 'Checking result'
d = ab.combine_shares(shares[0:n])
print "Recombined secret: %s" % d
print "Test %s" % ('successful' if d == secret else 'failed')
print '--------------------------------------'
