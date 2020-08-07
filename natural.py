def natural(formula):
# Inputs
# formula: a string representing the formula whose natural MID is being calculated

    #Import necessary packages and object classes
    from PolyMID import Formula
    from pdb import set_trace

    #Create a formula object
    formula = Formula(formula)

    #Determine the NaturalMID attribute of the formula object using the calc_natural_mid() method
    formula.calc_natural_mid()

    #Return the natural MID
    return(formula.NaturalMID)
