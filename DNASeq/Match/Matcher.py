## @package DNASeq.Match.Matcher
#  This package contains tools for DNA Sequence matching
from DNASeq.DNASequence import DNASequence
from bitarray import bitarray
import numpy as np


## DNA Sequence Matching functions
class Matcher():

    ## The constructor
    def __init__(self):
        pass

    ## Calculate base pair matches
    #  @type self: Matcher
    #  @param self: The matcher
    #
    #  @type seqa: DNASequence
    #  @param seqa: The DNA sequence to match
    #
    #  @type seqb: DNASequence
    #  @param seqb: The DNA sequence to match against
    #
    #  @rtype: np.array(int, int, bitarray)
    #  @return: Array containing number of matches for best match, offset, and bits (bool) whether the base pair matches.
    def matchDNASequence(self, seqa = DNASequence(), seqb = DNASequence()):
        matches = 0
        bestoffset = 0
        matcharr = bitarray()
        if(seqa.getBpCount() > 0 & seqb.getBpCount() > 0): # If there is something to match
            pass
        return np.array([matches, bestoffset, matcharr])

    ## Check matches for bitarrays
    #  @type self: Matcher
    #  @param self: The matcher
    #
    #  @type seqa: bitarray
    #  @param seqa: The array of bits to match
    #
    #  @type seqb: bitarray
    #  @param seqb: The array of bits to match against
    #
    #  @rtype:  np.array(int, bitarray)
    #  @return: Array containing the munber of matches and an array of bits for bp matches.
    def matchBitArray(self, seqa, seqb):
        matches = 0
        matcharr = bitarray()
        la = len(seqa)
        lb = len(seqb)
        if(la >= 2 & lb >= 2): # If there is something to match
            h = 0
            for i in range(0, la, 2):
                if(la > i):
                    h = i+1
                    if(seqa[i] == seqb[i] & seqa[h] == seqb[h]):
                        matches += 1
                        matcharr.append(True)
                    else:
                        matcharr.append(False)
        return np.array([matches, matcharr])
    


