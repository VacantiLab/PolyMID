def correct(fragment=None,TextFile=None,AtomLabeled='C'):

    import numpy as np
    import pdb
    from PolyMID import Fragment
    from PolyMID import Input
    from pdb import set_trace

    #Initialize the Inputs variable as an Input object
    Inputs = Input(fragment,TextFile,AtomLabeled)

    #Check if imputs were passed correctly and import the attributes to the Inputs variable from the fragment variable or the text file
    Inputs.check_and_import()
    if not Inputs.FragmentOrText:
        print('Either a Fragment object or a TextFile should be passed to CorrectMID.main(), but not both.')
        return

    #Create a Fragment object
    #    If a fragment object was passed to this function, this will be equivalent
    fragment = Inputs.fragment

    #Create a correction matrix if it was not an input
    #    the CM attribute of a Fragment object is a dictionary with keys corresponding to different atom identities
    if fragment.CM[AtomLabeled] is None:
        fragment.create_correction_matrix(Inputs.AtomLabeled)

    #If there is already a correction matrix, calculate its inverse
    #    CMi is not a dictionary and is calculated every time a correction is performed for an atom
    if fragment.CM[AtomLabeled] is not None:
        fragment.assign('CMi',np.linalg.pinv(fragment.CM[AtomLabeled]))

    #Calculate the corrected MID
    fragment.calc_corrected_mid()

    #Return the corrected MID
    print('\nThe corrected MID has been successfully calculated and returned as a numpy array.')
    return(fragment.MIDc)
