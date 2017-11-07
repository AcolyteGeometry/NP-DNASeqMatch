## @package: DNASeq.DNARand.RandDNA
#  Package for randomizing DNASequences for testing
from DNASeq.DNASequence import DNASequence
from Utils.MathHelpers import MathHelpers
from .RandBits import RandBits
from numpy.random import randint
import numpy as np


## Class for constructing random main DNASequence and DNASequences with overlapping segments
class RandDNA():

    rbits = RandBits()
    mathh = MathHelpers()

    ## The constructor
    def __init__(self):
        pass

    ## Generates a series of partially overlapping DNASequence objects
    #  @type self: RandDNA
    #  @param: The object
    #
    #  @type length: int
    #  @param length: The langth of the overall sequence bp
    #
    #  @type fraglen: int
    #  @param fraglen: The length of the sequence fragments bp
    #
    #  @type fragvar: float
    #  @param fragvar: Normal variation of the fragment length (0.0 - 1.0 [0-100%])
    #
    #  @type fragovmin: float
    #  @param fragovmin: The minimum normal overlap of fragment
    #
    #  @type fragovmax: float
    #  @param fragovmax: The maximum normal overlap of fragment
    #
    #  @rtype: np.array(bitarray, np.array(bitarray, bitarray, ...))
    #  @return: Array containing DNASequence complete sequence and array of DNASequence objects
    #  of overlapping subseuences.
    def getDNASeqSeries(self, length = 100, fraglen = 20, fragvar = 0.05, fragovmin = 0.1, fragovmax = 0.2):
        farrl = self.mathh.floatToIntUp((length / fraglen) / fragovmax) # Array size for output with overlap rounded up
        ret = np.array([None, np.array([None * farrl])]) # Initialize array
        ret.flags.writeable = True # Make arrays writable
        ret[1].flags.writeable = True
        ret[0] = self.getRandDNASequence(length) # Generate full sequence
        flen = randint(self.mathh.floatToInt(fraglen - (fraglen * fragvar)), \
                       self.mathh.floatToInt(fraglen + (fraglen * fragvar)), farrl) # Fragment lengthts
        folen = randint(self.mathh.floatToInt(fraglen * fragovmin),\
                        self.mathh.floatToInt(fraglen * fragovmax), farrl - 1) # Fragment overlap for each length -> next length
        ret[1][0] = ret[0].getDNASubseqBits(0, flen[0]) # Get initial fragment
        idx = flen[0]
        for i in range(1, farrl, 1): # Get fragments sequence with overlap
            ret[1][i] = ret[0].getDNASubseqBits((idx - folen[i]), flen[i])
            idx += flen[i]
        return ret

    ## Generates a random DNASequence
    #  @type self: RandDNA
    #  @param self: The object
    #
    #  @type length: int
    #  @param length: The number of base pairs
    #
    #  @rtype: DNASequence
    #  @return: The DNASequence generated
    def getRandDNASequence(self, length = 1):
        rseq = DNASequence()
        length = length * 2
        rseq.setDNABits(self.rbits.getRandBits(length))
        return rseq