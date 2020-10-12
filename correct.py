def Correct(fragment=None,TextFile=None):

    # Inputs:
    #     fragment: A fragment object containing information of the fragment whose MID is being corrected for natural isotopic abundances
    #     TextFile: A character that is the path to a .txt file containing information of the fragment whose MID is being corrected for natural isotopic abundances

    #     Defined within the fragment object:
    #         LabeledElement: A character that is the chemical symbol of the atom which is assumed to be labeled in the fragment whose MID is being corrected for natural isotopic abundances
    #         TracerEnrichment: The percent of the molecule assumed to be the tracer that is actually the tracer (e.g. 50% of glucose is [U-13C6]glucose)
    #         LabelEnrichment: The percent of the atom that is assumed to be labeled that is actually labeled (e.g. 99% of the atoms said to be 13C in [U-13C6]glucose are actually 13C)
    #         HighRes: A  list of strings indicating which elements have heavy isotopes whose mass differences are resolved from the mass differences due to heavy isotopes of the tracer element
    #             can also be 'all' or 'none'
    #             i.e. whether M1, M2, M3 are aggregate measurements of heavy isotopes of all atoms (low resolution, indicated by HighRes='none') or just the atom that is labeled and those not indicated in HighRes
    #             Note: Correcting high resolution data is achieved by setting all atom MIDs, for the elements indicated in HighRes, in the atom objects to [1 0 0] except that of the LabeledElement
    #                 This is accomplished in the definition of an Atom object where atom MIDs for high resolution data are taken from a separate file
    #             Note: this is now entered as a list separated by spaces in the input text file


    # Outputs:
    #     fragment.MIDc: The MID corrected for natural isotopic abundances

    import numpy as np
    import pandas as pd
    import pdb
    from PolyMID import Fragment
    from PolyMID import Input
    from PolyMID import Tracer
    from pdb import set_trace

    #Initialize the Inputs variable as an Input object
    Inputs = Input(fragment=fragment,TextFile=TextFile)

    #Check if imputs were passed correctly and import the attributes to the Inputs variable from the fragment variable or the text file
    #    One of the arguemtns fragment or TextFile to correct() should be None. Both cannot have values.
    #    If fragment has a value other than None, it will be used to define Input.fragment
    #    If TextFile has a value, it will be read as a directory to a .txt file to reference to define Input.fragment
    #        The format of this .txt file is provided in the References folder
    #    If neither has a value, the user will be prompted to provide a .txt file with the fragment definition
    #        The format of this .txt file is provided in the References folder
    Inputs.check_and_import()
    if not Inputs.FragmentOrText:
        print('Either a Fragment object or a TextFile should be passed to CorrectMID.main(), but not both.')
        return

    print('\nCalculating corrected MIDs...')

    #Create a Fragment object
    #    If a fragment object was passed to this function, this will be equivalent
    fragment = Inputs.fragment

    #Create a correction matrix if it was not an input
    #    the CM attribute of a Fragment object is a dictionary with keys corresponding to different atom identities
    LabeledElement = fragment.Tracer.LabeledElement
    if fragment.CM[LabeledElement] is None:
        fragment.create_correction_matrix()

    #If there is already a correction matrix, calculate its inverse
    #    CMi is not a dictionary and is calculated every time a correction is performed for an atom
    if fragment.CM[LabeledElement] is not None:
        fragment.assign('CMi',np.linalg.pinv(fragment.CM[LabeledElement]))

    #Calculate the corrected MID
    fragment.calc_corrected_mid()

    #Return the corrected MID
    print('The corrected MIDs are as follows:\n')
    MIDc = pd.DataFrame(fragment.MIDc)
    MIDc = MIDc.transpose()
    MIDc.columns=['M0','M1','M2','M3']
    MIDc.index=[fragment.name]
    print(MIDc)
    print('\n\n')

    return(fragment.MIDc)
