## @package: Utils.MathHelpers
#  Math helpers for calculating various parts of the TSP


## Class contining methods related to number format conversion
#  for example changing a Float into rounded Integer.
class MathHelpers():

    ## The Constructor.
    def __init__(self):
        pass # Do nothing.

    ## String to integer
    #  @type self: MathHelpers
    #  @param self: The object
    #
    #  @type str: string
    #  @param str: the integer string to convert
    #
    #  @rtype: int
    #  @return: The integer from string or 0 on error.
    def stoi(self, str):
        if (str.isdigit()):
            str = int(str)
            return str
        return 0

    ## Converts a float to a rounded integer
    # @type self: MathHelpers
    # @param self: The object.
    #
    # @type fl: Float, String
    # @param fl: Float to convert to integer
    #
    # @rtype: Integer
    # @return: An integer rounded to one place, or 0 if input is invalid.
    def floatToInt(self, fl):
        if(isinstance(fl, float)): # If this is a float
            if(fl.is_integer()): # If it is a whole number return as integer;
                return int(fl)
            else: # Otherwise split at the decimal.
                intarr = str(fl).split(".")
                if(int(intarr[1]) >=5): # Read 1/10ths place, check if >= 5.
                    return int(intarr[0]) + 1 # Return integer rounded up.
                else:
                    return int(intarr[0]) # Otherwise return integer rounded down.
        elif(isinstance(fl, str)): # If float is a string, convert to float and recurse.
            return self.floatToInt(float(fl))
        elif(isinstance(fl, int)): # If float is already an integer just return it.
            return fl
        return 0 # Return 0 if the float is invalid.

    ## Converts a float rounded up to the nearest even integer.
    # @type self: MathHelpers
    # @param self: The object.
    #
    # @type fl: Float, String
    # @param fl: Float to convert to integer
    #
    # @rtype: Integer
    # @return: An even integer rounded to one place, or 0 if input is invalid.
    def floatToInt2(self, fl):
        ret = self.floatToInt(fl)
        if((ret % 2) != 0):
            ret += 1
        return ret


    ## Converts a float rounded down to the nearest even integer.
    # @type self: MathHelpers
    # @param self: The object.
    #
    # @type fl: Float, String
    # @param fl: Float to convert to integer
    #
    # @rtype: Integer
    # @return: An even integer rounded to one place, or 0 if input is invalid.
    def floatToIntm2(self, fl):
        ret = self.floatToInt(fl)
        if((ret % 2) != 0):
            ret -= 1
        return ret

    ## Converts a float to an int rounded up
    #  @type self: MathHelpers
    #  @param self: The object
    #
    #  @type fl: float, string
    #  @param fl: Float to be converted to integer
    #
    #  @rtype: int
    #  @return: An integer rounded up from to the nearest whole number
    def floatToIntUp(self, fl):
        if(isinstance(fl, float)): # If this is a float
            if(fl.is_integer()): # If it is a whole number return as integer;
                return int(fl)
            else: # Otherwise split at the decimal.
                intarr = str(fl).split(".")
                if(len(intarr) == 1 & int(intarr[1]) != 0): # Read 1/10ths place, check if >= 5.
                    return int(intarr[0]) + 1 # Return integer rounded up.
                else:
                    return int(intarr[0]) # Otherwise return integer rounded down.
        elif(isinstance(fl, str)): # If float is a string, convert to float and recurse.
            return self.floatToInt(float(fl))
        elif(isinstance(fl, int)): # If float is already an integer just return it.
            return fl
        return 0 # Return 0 if the float is invalid.
