def main(fragment=None,TextFile=None):

    import pdb
    import numpy as np
    import pdb
    from Pesciolini import Fragment
    from Pesciolini import Input

    Inputs = Input(fragment,TextFile)

    Inputs.check()
    if not Inputs.FragmentOrText:
        print('Either a Fragment object or a TextFile should be passed to CorrectMID.main(), but not both.')
        return

    #create a Fragment object if one was not passed to this function
    if (fragment is None):
        #initialize fragment object attributes not taken as function inputs
        MIDc=None
        PeakArea=None

        #creates a fragment object
        fragment = Fragment(formula=Inputs.formula, MetaboliteAtoms=Inputs.MetaboliteAtoms, MIDu=Inputs.MIDu, AtomLabeled=Inputs.AtomLabeled, FragmentName=Inputs.FragmentName, CM=Inputs.CM, MIDc=MIDc, PeakArea=PeakArea)


    #create a correction matrix if one is not provided in the Fragment object
    if fragment.CM is None:
        fragment.create_correction_matrix()

    #calculate the corrected MID
    fragment.calc_corrected_mid()

    #return the corrected MID
    return(fragment.MIDc)
