## @package: DNASeq.DNARand.RandDNA
#  Package for randomizing DNASequences for testing
from DNASeq.DNASequence import DNASequence
from .RandBits import RandBits


## Class for constructing random main DNASequence and DNASequences with overlapping segments
class RandDNA():

    rbits = RandBits()

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
    #  @return: The DNASequence generated
    def getRandDNASequence(self, length = 1):
        rseq = DNASequence()
        length = length * 2
        rseq.setDNABits(self.rbits.getRandBits(length))
        return rseq