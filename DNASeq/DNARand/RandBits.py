## @package: DNASeq.DNARand.Random
#  Methods for generating random DNASequence.
from DNASeq.DNASequence import DNASequence
from Utils.MathHelpers import MathHelpers
from numpy.random import bytes as randbytes
from bitarray import bitarray

## Class for generating random bit array sequences
class RandBits():

    mathh = MathHelpers()

    ## The constructor
    def __init__(self):
        pass

    ## Generates a random DNASequence
    #  @type self: RandBits
    #  @param self: The object
    #
    #  @type length: int
    #  @param length: The number of base pairs
    #
    #  @rtype: DNASequence
    #  @return: The DNASeuence generated
    def getRandDNASequence(self, length = 1):
        rseq = DNASequence()
        length = length * 2
        rseq.setDNABits(self.getRandBits(length))
        return rseq

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
                bits.append(randbytes(1))
        # and bits
        sz = amt % 8
        abyte = bitarray(randbytes(1))
        for i in range(sz):
            bits.append(abyte[i])
        return bits