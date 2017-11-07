## @package DNASeq.DNASequence
#  A class and methods for holding and querying a DNA sequence
from Utils.MathHelpers import MathHelpers
from .BasePair import BasePair
from bitarray import bitarray



## Class representing a sequence of DNA
class DNASequence():

    mathh = MathHelpers()
    bp = BasePair() # Base pair helper class object
    # This will hold a bit array of a R/DNA sequence.
    # A => 00
    # T => 01
    # U => 01
    # C => 10
    # G => 11
    dna = bitarray()
    size = 0 # Size of sequence
    inverted = False # If the sequence is inverted from input
    reversed = False # If the sequence is reversed from input

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
        self.inverted = False
        self.reversed = False
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
        if(length == float('inf') | length > self.size):
            length = self.size - start
        start = self.mathh.floatToIntm2(start * 2) # Since each base pair is 2 bits
        length = self.mathh.floatToIntm2(length * 2) # We double both values
        if(start < self.size): # Make sure we aren't starting out past the end of seuence.
            for i in range(start, length, 2): # For every even bit
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
        l = self.mathh.floatToInt(len(self.dna))
        if(l % 2 != 0):
            del self.dna[-1]
            l -= 1
            self.size = self.mathh.floatToIntm2(l * 0.5)
        return self.size

    ## Sets the length in base pairs, fixes alignment and length if off.
    #  @type self: DNASequence
    #  @param self: The object
    #
    #  @rtype: int
    #  @return: The number of base pairs in sequence
    def setBpCount(self):
        l = self.mathh.floatToIntm2(len(self.dna))
        if(l != self.mathh.floatToInt(len(self.dna))):
            return self.fixBpAlignment()
        return self.size

    ## Inverts the DNA sequence: A => T, T => A, C => G and G => C
    #  @type self: DNASequence
    #  @param self: The DNA Sequence
    #
    #  @rtype: bool
    #  @return: Whether the DNA sequence is inverted from original input.
    def invertDNA(self):
        self.setBpCount() # Make sure count is correct and check/repair offset
        for i in range(1, len(self.dna), 2): # invert every odd bit
            self.dna[i] = self.dna[i] ^ True
        self.inverted = self.inverted ^ True
        return self.inverted

    ## Reverses the DNA sequence: ATGC => CGTA
    #  @type self: DNASequence
    #  @param self: The DNA Sequence
    #
    #  @rtype: bool
    #  @return: Whether the DNA seuence is reversed from original input.
    def reverseDNA(self):
        self.setBpCount() # Make sure count is correct and set/repair offset
        rev = bitarray() # Temp hold reversed array
        h = 0;
        for i in range(len(self.dna), 0, -2): # For each base pair in reverse
            h = i - 1
            rev.append([self.dna[h], self.dna[i]]) # Add bp bits in reverse
        self.dna = rev
        self.reversed = self.reversed ^ True
        return self.reversed

    ## Returns whether the sequence is inverted
    #  @type self: DNASequence
    #  @param self: The DNA Sequence
    #
    #  @rtype: bool
    #  @return: If the DNA sequence is inverted from original input.
    def isInverted(self):
        return self.inverted

    ## Returns whether the sequence is reversed
    #  @type self: DNASequence
    #  @param self: The DNA Sequence
    #
    #  @rtype: bool
    #  @return: If the DNA sequence is reversed from original input.
    def isReversed(self):
        return self.reversed



