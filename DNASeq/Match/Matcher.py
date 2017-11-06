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

    '''
    TODO: Invert matching
    TODO: Reverse sequence matching
    '''

    ## Calculate base pair matches, checks regular, inverted, reversed and inverted, and reversed.
    #  @type self: Matcher
    #  @param self: The matcher
    #
    #  @type seqa: DNASequence
    #  @param seqa: The DNA sequence to match
    #
    #  @type seqb: DNASequence
    #  @param seqb: The DNA sequence to match against
    #
    #  @type min: float
    #  @param min: Normal percentage for minimum match overlap (default 0.01 == 1%)
    #
    #  @rtype: np.array(float, int, int, bitarray)
    #  @return: Array containing best match normal percentage, number of matches, offset, and bits (bool) whether index base pair matches.
    def matchDNASequence(self, seqa = DNASequence(), seqb = DNASequence(), min = 0.01):
        matches = 0
        bestoffset = 0
        matcharr = bitarray()
        matchnormal = 0.0
        if(seqa.getBpCount() > 0 & seqb.getBpCount() > 0): # If there is something to match
            # Regular matching
            mbest = self.matchBitArray(seqa.getDNABits(), seqb.getDNABits()) # Match with no offset
            matches = mbest[0]
            matcharr = mbest[1]
            matchnormal = float(matches) / float(seqa.getBpCount())
            mtest = None
            for i in range(1, seqa.getBpCount() - seqa.getBpCount() * min): # for bp offset in seqa (Direction 1)
                mtest = self.matchBitArray(seqa.getDNASubseqBits(i), seqb.getDNABits())
                if(float(mtest[0]) / float(seqa.getBpCount() - i) > matchnormal): # If this is a better match
                    matches = mtest[0]
                    matcharr = mtest[1]
                    matchnormal = float(matches) / float(seqa.getBpCount() - i)
            for i in range(1, seqb.getBpCount() - seqb.getBpCount() * min): # for bp offset in seqb (Direction 2)
                mtest = self.matchBitArray(seqb.getDNASubseqBits(i), seqa.getDNABits())
                if(float(mtest[0]) / float(seqb.getBpCount() - i) > matchnormal): # If this is a better match
                    matches = mtest[0]
                    matcharr = mtest[1]
                    matchnormal = float(matches) / float(seqb.getBpCount() - i)

            # Inverted matching
            seqa.invertDNA()
            for i in range(0, seqa.getBpCount() - seqa.getBpCount() * min): # for bp offset in seqa (Direction 1)
                mtest = self.matchBitArray(seqa.getDNASubseqBits(i), seqb.getDNABits())
                if(float(mtest[0]) / float(seqa.getBpCount() - i) > matchnormal): # If this is a better match
                    matches = mtest[0]
                    matcharr = mtest[1]
                    matchnormal = float(matches) / float(seqa.getBpCount() - i)
            for i in range(0, seqb.getBpCount() - seqb.getBpCount() * min): # for bp offset in seqb (Direction 2)
                mtest = self.matchBitArray(seqb.getDNASubseqBits(i), seqa.getDNABits())
                if(float(mtest[0]) / float(seqb.getBpCount() - i) > matchnormal): # If this is a better match
                    matches = mtest[0]
                    matcharr = mtest[1]
                    matchnormal = float(matches) / float(seqb.getBpCount() - i)

            # Inverted and Reversed matching
            seqa.reverseDNA()
            for i in range(0, seqa.getBpCount() - seqa.getBpCount() * min): # for bp offset in seqa (Direction 1)
                mtest = self.matchBitArray(seqa.getDNASubseqBits(i), seqb.getDNABits())
                if(float(mtest[0]) / float(seqa.getBpCount() - i) > matchnormal): # If this is a better match
                    matches = mtest[0]
                    matcharr = mtest[1]
                    matchnormal = float(matches) / float(seqa.getBpCount() - i)
            for i in range(0, seqb.getBpCount() - seqb.getBpCount() * min): # for bp offset in seqb (Direction 2)
                mtest = self.matchBitArray(seqb.getDNASubseqBits(i), seqa.getDNABits())
                if(float(mtest[0]) / float(seqb.getBpCount() - i) > matchnormal): # If this is a better match
                    matches = mtest[0]
                    matcharr = mtest[1]
                    matchnormal = float(matches) / float(seqb.getBpCount() - i)

            # Reversed matching
            seqa.invertDNA() # Revert inversion
            for i in range(1, seqa.getBpCount() - seqa.getBpCount() * min): # for bp offset in seqa (Direction 1)
                mtest = self.matchBitArray(seqa.getDNASubseqBits(i), seqb.getDNABits())
                if(float(mtest[0]) / float(seqa.getBpCount() - i) > matchnormal): # If this is a better match
                    matches = mtest[0]
                    matcharr = mtest[1]
                    matchnormal = float(matches) / float(seqa.getBpCount() - i)
            for i in range(1, seqb.getBpCount() - seqb.getBpCount() * min): # for bp offset in seqb (Direction 2)
                mtest = self.matchBitArray(seqb.getDNASubseqBits(i), seqa.getDNABits())
                if(float(mtest[0]) / float(seqb.getBpCount() - i) > matchnormal): # If this is a better match
                    matches = mtest[0]
                    matcharr = mtest[1]
                    matchnormal = float(matches) / float(seqb.getBpCount() - i)
        return np.array([matchnormal, matches, bestoffset, matcharr])

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



