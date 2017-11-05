## @package DNASeq.BassePair
#  A class and methods that define a base pair
from bitarray import bitarray
import copy

## Defines an R/DNA base pair.
class BasePair():

    bases = [['A', bitarray('00')], \
             ['T', bitarray('01')], \
             ['U', bitarray('01')], \
             ['C', bitarray('10')], \
             ['G', bitarray('11')]]

    ## The constructor
    def __init__(self):
        pass

    ## Returns a bitarray containing the binary representation of a base pair
    #  @type self: BasePair
    #  @param self: The object
    #
    #  @type base: str
    #  @param base: The letter representing the base pair A, T, U, C, G
    #
    #  @rtype: bitarray
    #  @return: The bitarray containing the bits representing the base bair, or an empty bitarray on err
    def getBpBits(self, base):
        if(isinstance(base, str)):
            for i in self.bases:
                if(i[0] == base):
                    return bitarray(i[1])
        return bitarray()

    ## Returns a string containing the character representation of a base pair
    #  @type self: BasePair
    #  @param self: The object
    #
    #  @type base: bitarray
    #  @param base: The bitarray containing the bits representing the base bair
    #
    #  @rtype: str
    #  @return: The letter representing the base pair A, T, U, C, G , or an empty string on err
    def getBpChar(self, base):
        if (isinstance(base, bitarray) & len(base) == 2):
            for i in self.bases:
                if (i[1][0] == base[0] & i[1][1] == base[1]):
                    return i[0]
        return ''

