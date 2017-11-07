## @package: DNASeq.DNARand.Random
#  Generates random bitarrays
from Utils.MathHelpers import MathHelpers
from numpy import random
from bitarray import bitarray

## Class for generating random bit array sequences
class RandBits():

    mathh = MathHelpers()

    ## The constructor
    def __init__(self):
        pass

    ## Generate random bits
    #  @type self: RandBits
    #  @param self: The object
    #
    #  @param amt: int
    #  @param amt: The number of bits to generate
    #
    #  @rtype: bitarray
    #  @return: Bit array of random bits
    def getRandBits(self, amt = 2):
        bits = bitarray()
        if(amt >= 8): # 8 bit bytes
            for i in range(self.mathh.floatToInt(float(amt) * 0.125)):
                bits.frombytes(random.bytes(1))
        # and bits
        sz = amt % 8
        abyte = bitarray()
        abyte.frombytes(random.bytes(1))
        for i in range(sz):
            bits.append(abyte[i])
        return bits