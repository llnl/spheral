#ATS:test(SELF, label="Unit tests of uniform_random random number generator")
from math import *
from Spheral import *
from SpheralTestUtilities import *

import unittest

# Create a global random number generator.
import random
random.seed(4599281940)

ntests = 1000

class TestRandom01(unittest.TestCase):

    #===========================================================================
    # Various ways of constructing
    #===========================================================================
    def testConstructors(self):
        seed1 = random.randint(1, 2**64)
        seed3 = random.randint(1, 2**64)
        while seed3 == seed1:
            seed3 = random.randint(1, 2**64)
        gen1 = uniform_random(seed1)
        gen2 = uniform_random(gen1)
        gen3 = uniform_random(seed3)
        assert gen1 == gen2
        assert gen1 != gen3
        for i in range(ntests):
            assert gen1() == gen2()

    #===========================================================================
    # seed
    #===========================================================================
    def testSeed(self):
        seed = random.randint(1, 2**64)
        gen1 = uniform_random(seed)
        assert gen1.seed == seed
        gen2 = uniform_random()
        assert gen1 != gen2
        gen2.seed = seed
        assert gen2.seed == seed
        for i in range(ntests):
            assert gen1() == gen2()

    #===========================================================================
    # Comparisons
    #===========================================================================
    def testComparisons(self):
        seed = random.randint(1, 2**64)
        gen1 = uniform_random(seed)
        gen2 = uniform_random(seed + 1)
        assert gen1 != gen2
        gen2.seed = seed
        assert gen1 == gen2
        gen3 = uniform_random(seed, 2.0, 3.0)
        assert gen3 != gen1
        gen1.range(2.0, 3.0)
        assert gen3 == gen1

    #===========================================================================
    # advance
    #===========================================================================
    def testAdvance(self):
        seed = random.randint(1, 2**64)
        gen1 = uniform_random(seed)
        throwaway = [gen1() for i in range(ntests)]
        vals1 = [gen1() for i in range(ntests)]
        gen2 = uniform_random(seed)
        gen2.advance(ntests)
        vals2 = [gen2() for i in range(ntests)]
        assert vals1 == vals2

    #===========================================================================
    # range
    #===========================================================================
    def testRange(self):
        seed = random.randint(1, 2**64)
        gen1 = uniform_random(seed)
        assert gen1.min == 0.0
        assert gen1.max == 1.0
        gen1.range(5.0, 10.0)
        assert gen1.min == 5.0
        assert gen1.max == 10.0

    #===========================================================================
    # Serialization
    #===========================================================================
    def testSerialize(self):
        seed = random.randint(1, 2**64)
        gen1 = uniform_random(seed)
        throwaway = [gen1() for i in range(ntests)]
        buf = vector_of_char()
        gen1.serialize(buf)
        gen2 = uniform_random()
        i = gen2.deserialize(buf, 0)
        assert i == len(buf)
        vals1 = [gen1() for i in range(ntests)]
        vals2 = [gen2() for i in range(ntests)]
        assert vals1 == vals2

if __name__ == "__main__":
    unittest.main()

