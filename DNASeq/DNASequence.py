## @package DNASeq.DNASequence
#  A class and methods for holding and querying a DNA sequence
from Utils.MathHelpers import MathHelpers
from .BasePair import BasePair
from bitarray import bitarray



## Class representing a sequence of DNA
class DNASequence():

    mathh = MathHelpers()
    bp = BasePair() # Base pair helper class object
    dna = bitarray()
    # This will hold a bit array of a R/DNA sequence.
    # A => 00
    # T => 01
    # U => 01
    # C => 10
    # G => 11
    size = 0 # Size of sequence

    ## The constructor
    def __init__(self):
        pass

    ## Gets the dna bitarray
    #  @type self: DNASequence
    #  @param self: The class
    #
    #  @rtype: bitarray
    #  @return: The complete DNA sequence
    def getDNABits(self):
        return self.dna

    ## Sets the dna bitarray with base pair alignment
    #  @type self: DNASequence
    #  @param self: The class
    #
    #  @type seq: bitarray
    #  @param seq: The DNA sequence as bitarray
    #
    #  @rtype: int
    #  @return: The size of DNASeuence in base pairs
    def setDNABits(self, seq = bitarray()):
        self.dna = seq
        return self.setBpCount()

    ## Gets the length in base pairs
    #  @type self: DNASequence
    #  @param self: The object
    #
    #  @rtype: float
    #  @return: The number of base pairs in sequence
    def getBpCount(self):
        return self.size

    ## Gets a number of base pairs from a DNA sequence
    #  @type self: DNASequence
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
    def getDNASubseqBits(self, start = 0, length = float('inf')):
        odna = bitarray() # Array of bits to return
        tba = bitarray(2) # Temp array of 2 bits for appending
        tba.setall(False) # Initialize false
        start = self.mathh.floatToInt2(start * 2.0) # Since each base pair is 2 bits
        length = self.mathh.floatToInt2(length * 2.0) # We double both values
        if(start < self.size): # Make sure we aren't starting out past the end of seuence.
            for i in range(start,(start + length), 2): # For every even bit
                tba[0] = self.dna[i] # set the temp bit and the next (odd) bit
                tba[1] = self.dna[i+1]
                odna.append(tba) # Append bits to bitarray
        return odna

    ## Append the R/DNA Sequence
    #  @type self: DNASequence
    #  @param self: The object
    #
    #  @type seq: bitarray
    #  @param seq: bitarray containing the DNA bp sequence to append. If the array is not an even number, the last bit will be discarded.
    #
    #  @rtype: int
    #  @return: The new length in base pairs
    def appendSeqBits(self, seq = bitarray()):
        if(len(seq) >= 2 ): # If there is something worth adding
            self.dna.append(seq) # actually append the sequence
        self.setBpCount()
        return self.size

    ## Append the R/DNA Sequence
    #  @type self: DNASequence
    #  @param self: The object
    #
    #  @type seq: str
    #  @param seq: String containing the DNA bp sequence to append.
    #
    #  @rtype: int
    #  @return: The new length in base pairs
    def appendSeqChars(self, seq=''):
        if (len(seq) >= 1 & isinstance(seq, str)):  # If there is something worth adding
            bseq = bitarray()
            for i in seq:
                bseq.append(self.bp.getBpBits(i))
            self.dna.append(bseq)  # actually append the sequence
        self.setBpCount()
        return self.size

    ## Fix the alignment of base pairs (removes last bit in dna bitarray if number of bits is uneven)
    #  @type self: DNASequence
    #  @param self: The object
    #
    #  @rtype: int
    #  @return: The number of base pairs in DNASequence
    def fixBpAlignment(self):
        l = self.mathh.floatToInt2(len(self.dna))
        if(l % 2 != 0):
            del self.dna[-1]
            l -= 1
            self.size = self.mathh.floatToInt2(l * 0.5)
        return self.size

    ## Sets the length in base pairs, fixes alignment and length if off.
    #  @type self: DNASequence
    #  @param self: The object
    #
    #  @rtype: int
    #  @return: The number of base pairs in sequence
    def setBpCount(self):
        l = self.mathh.floatToInt2(len(self.dna))
        if(l != self.mathh.floatToInt(len(self.dna))):
            return self.fixBpAlignment()
        return self.size




