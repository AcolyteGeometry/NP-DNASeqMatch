## @package DNASeq.Seuence
#  A class and methods for holding and querying a DNA sequence
from bitarray import bitarray
from .BasePair import BasePair
from Utils.MathHelpers import MathHelpers

## Class representing a sequence of DNA
class Sequence():

    mathh = MathHelpers()
    dna = bitarray()
    # This will hold a bit array of a R/DNA sequence.
    # A => 00
    # T => 01
    # U => 01
    # C => 10
    # G => 11

    ## The constructor
    def __init__(self):
        pass

    ## Gets the dna bitarray
    #  @type self: Sequence
    #  @param self: The class
    #
    #  @rtype: bitarray
    #  @return: The complete DNA sequence
    def getDNA(self):
        return self.dna

    ## Gets the length in base pairs
    #  @type self: Sequence
    #  @param self: The object
    #
    #  @rtype: float
    #  @return: The number of base pairs in sequence
    def getBpCount(self):
        l = (len(self.dna) * 1.0)
        if(l % 2.0 == 0.0):
            return l * 0.5
        else:
            return (l * 0.5) - 0.5

    ## Gets a number of nucleotides from a DNA sequence
    #  @type self: Sequence
    #  @param self: The class
    #
    #  @type start: float
    #  @param start: The starting index (0 - inf)
    #
    #  @type length: float
    #  @param length: The number of nucleotides to return
    #
    #  @rtype: bitarray
    #  @return: The subsequence of DNA, up to the length of the sequence from index.
    def getDNASubseq(self, start = 0.0, length = float('inf')):
        odna = bitarray() # Array of bits to return
        tba = bitarray(2) # Temp array of 2 bits for appending
        tba.setall(False) # Initialize false
        start = start * 2.0 # Since each base pair is 2 bits
        length = length * 2.0 # We double both values
        if(start < len(self.dna)): # Make sure we aren't starting out past the end of seuence.
            for i in range(self.mathh.floatToInt2(start),\
                           self.mathh.floatToInt2(start + length), 2): # For every even bit
                tba[0] = self.dna[i] # set the temp bit and the next (odd) bit
                tba[1] = self.dna[i+1]
                odna.append(tba) # Append bits to bitarray
        return odna

    ## Append the R/DNA Sequence
    #  @type self: Sequence
    #  @param self: The object
    #
    #  @type seq: bitarray
    #  @param seq: bitarray containing the DNA bp sequence to append. If the array is not an even number, the last bit will be discarded.
    #
    #  @rtype: float
    #  @return: The new length in base pairs
    def appendSeq(self, seq = bitarray()):
        l = (len(self.dna) * 2.0)
        if(len(seq) >= 2 ): # If there is something worth adding
            if(len(seq) % 2 != 0): # If not an even number of bits
                del seq[-1] # Remove last bit to align to bp
            self.dna.append(seq) # actually append the sequence
            l += (len(seq) * 2.0) # new length
        return l


