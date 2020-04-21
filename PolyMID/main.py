def main(fragment=None,TextFile=None):

    import pdb
    import numpy as np
    import pdb
    from Pesciolini import Fragment
    from Pesciolini import Input

    #Initialize the Inputs variable as an Input object
    Inputs = Input(fragment,TextFile)

    #Check if imputs were passed correctly and import the attributes to the Inputs variable from the fragment variable or the text file
    Inputs.check_and_import()
    if not Inputs.FragmentOrText:
        print('Either a Fragment object or a TextFile should be passed to CorrectMID.main(), but not both.')
        return

    #create a Fragment object if one was not passed to this function
    fragment = Inputs.fragment

    #create a correction matrix if one is not provided in the Fragment object
    if fragment.CM is None:
        fragment.create_correction_matrix()

    #calculate the corrected MID
    fragment.calc_corrected_mid()

    #return the corrected MID
    return(fragment.MIDc)
