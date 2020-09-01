def natural(formula,AtomLabeled='C',HighRes=False):
    # Inputs:
    # formula: a string representing the formula whose natural MID is being calculated
    # AtomLabeled:
    # HighRes:

    # Outputs:
    # formula.NaturalMID

    #Import necessary packages and object classes
    from PolyMID import Formula
    from pdb import set_trace

    #Create a formula object
    formula = Formula(formula=formula,AtomLabeled=AtomLabeled,HighRes=HighRes)
    #     The Formula object designation takes AtomLabeled and HighRes inputs
    #         If HighRes is False, AtomLabeled does not matter (which is default behaviour)
    #         High resolution natural MIDs can be calculated by specifying HighRes and AtomLabeled, where AtomLabeled is the atom contributing heavy isotopes

    #Determine the NaturalMID attribute of the formula object using the calc_natural_mid() method
    formula.calc_natural_mid()

    #Return the natural MID
    return(formula.NaturalMID)
