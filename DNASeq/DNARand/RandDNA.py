## @package: DNASeq.DNARand.RandDNA
#  Package for randomizing DNASequences for testing
from DNASeq.DNASequence import DNASequence
from Utils.MathHelpers import MathHelpers
from DNASeq.DNARand.RandBits import RandBits
from numpy.random import randint


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
    def getDNASeqSeries(self, length = 100, fraglen = 20, fragvar = 0.05, fragovmin = 0.3, fragovmax = 0.5):
        farrl = self.mathh.floatToInt((length / fraglen) * (1 + fragovmax))# Array size for output with overlap rounded up
        ret = [None, [None for i in range(farrl)]] # Initialize array
        dna = self.getRandDNASequence(length) # Generate full sequence
        ret[0] = dna
        fllow = self.mathh.floatToInt(fraglen - self.mathh.floatToIntUp(fraglen * fragvar))
        flhigh = self.mathh.floatToIntUp(fraglen + self.mathh.floatToIntUp(fraglen * fragvar))
        if(fllow == flhigh):
            fllow -=1
            flhigh += 1
        flen = randint(fllow, flhigh, farrl) # Fragment lengthts
        fllow = self.mathh.floatToInt(fraglen * fragovmin)
        flhigh = self.mathh.floatToInt(fraglen * fragovmax)
        if(fllow == flhigh):
            fllow -= 1
            flhigh += 1
        folen = randint(fllow, flhigh, farrl - 1) # Fragment overlap for each length -> next length
        idx = flen[0]
        dna2 = DNASequence()
        dna2.setDNABits(dna.getDNASubseqBits(0, idx)) # Get initial fragment
        ret[1][0] = dna2
        for i in range(1, farrl, 1): # Get fragments sequence with overlap
            if(idx < dna.size):
                dna2 = DNASequence()
                dna2.setDNABits(ret[0].getDNASubseqBits(idx, flen[i]))
                ret[1][i] = dna2
                idx += (flen[i] - folen[i-1])
        arrt = []
        for i in range(len(ret[1])):
            if(ret[1][i] != None):
                arrt.append(ret[1][i])
        ret[1] = arrt
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

'''
rdna = RandDNA()
dnas = rdna.getDNASeqSeries(3000, 100, 0.1, 0.05, 0.3)
print(str(dnas))
print(dnas[0].getDNABits().to01())
for i in dnas[1]:
    print(i.getDNABits().to01())
'''