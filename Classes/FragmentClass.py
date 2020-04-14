class Fragment:

    # Initializer and Instance Attributes
    def __init__(self,formula,MetaboliteAtoms,MIDu,AtomLabeled,FragmentName,CM,MIDc,PeakArea):
        from Pesciolini import Formula
        import numpy as np

        self.name = FragmentName
        self.Formula = Formula(formula)
        self.MetaboliteAtoms = MetaboliteAtoms
        self.MIDu = MIDu
        self.MIDc = MIDc
        self.CM = CM
        self.CMi = None
        self.PeakArea = PeakArea
        self.AtomLabeled = AtomLabeled

        if CM != None:
            self.CMi = np.linalg.pinv(CM) #find the right inverse (pseudo-inverse in numpy jargon) of the correction matrix

    # instance method to assign new values to a Fragment object
    def assign(self,attribute,NewValue):
        if attribute == 'name':
            self.name = NewValue
        if attribute == 'formula':
            self.formula = NewValue
        if attribute == 'MetaboliteAtoms':
            self.MetaboliteAtoms = NewValue
        if attribute == 'MIDu':
            self.MIDu = NewValue
        if attribute == 'MIDc':
            self.MIDc = NewValue
        if attribute == 'AtomLabeled':
            self.AtomLabeled = NewValue
        if attribute == 'CM':
            self.CM = NewValue
        if attribute == 'PeakArea':
            self.PeakArea = NewValue

    # instance method to create correction matrix for a Fragment object
    def create_correction_matrix(self):

        import pdb
        import numpy as np
        import re
        import copy
        import pandas
        from Pesciolini import quantity_of_atom
        from Pesciolini import Formula

        pdb.set_trace
        #break the formula up so atoms and quantities are consecutive entries in a numpy array
        broken_formula = np.array(re.findall('[A-Z][a-z]?|[0-9]+', self.Formula.formula))
        #    '[A-Z][a-z]?|[0-9]+': A capital letter [A-Z], followed by a lowercase letter [a-z], which is optional '?', or '|' a number '[0-9]', and possibly more numbers '+'
        #    example: this command will take formula = C6H12N1O3Si1 and return broken_formula = array(['C','6','H','12','N','1','O','3','Si','1'])
        #        all components are strings
        n_formula_entries = len(broken_formula)

        #find the index of the number of the atom which can acquire a label in the formula for the full fragment
        #    it is used in creating the correction matrix below because successive quantities of this atom need to be subtracted and a heavy atom put in its place
        atom_index = np.where(broken_formula==self.AtomLabeled)[0][0]
        atom_quantity_index = atom_index+1 #refering to full fragment

        #the number of rows of the correction matrix is equal to the quantity of the atom being corrected for that are in the fragment and the original metabolite
        atom_quantity = quantity_of_atom(self.MetaboliteAtoms,self.AtomLabeled) #this does not refer to the full fragment!

        #add the "heavy atom to the end of the broken formula array", initially its quantity is 0
        broken_formula = np.append(broken_formula,np.array(['Hv','0']))
        n_formula_entries_heavy = len(broken_formula)

        #replace each atom of interest with a heavy atom and get the natural mid of the result
        #    these mids fill the rows of the correction matrix
        broken_formula_correct = copy.copy(broken_formula) #initialize the array to carry the formula with a heavy atom
        correction_matrix_dict = dict() #initialize a dictionary to hold the rows of the correction matrix
        n_theoretical_mid_entries = atom_quantity + 4 #the theoretical length is all possible atoms to label plus 4 (after this the relative abundances are assumed negligible)
        for i in range(0,atom_quantity+1):
            #subtract an atom of interest from the formula
            broken_formula_correct[atom_quantity_index] = broken_formula[atom_quantity_index].astype(np.int) - i
            broken_formula_correct[atom_quantity_index] = broken_formula_correct[atom_quantity_index].astype(np.str)

            #replace that atom with a heavy atom
            broken_formula_correct[n_formula_entries+1] = broken_formula[n_formula_entries+1].astype(np.int) + i
            broken_formula_correct[n_formula_entries+1] = broken_formula_correct[n_formula_entries+1].astype(np.str)

            #update the string version of the formula from the array version
            new_formula = ''
            for j in range(0,n_formula_entries_heavy):
                new_formula = new_formula + broken_formula_correct[j]

            #make a Formula object
            new_formula = Formula(new_formula)

            #get the mid due to natural abundances of the updated formula (with one or more heavy atoms)
            new_formula.calc_natural_mid()
            correction_matrix_dict[i] = new_formula.NaturalMID

            #shorten each theoretical MID with given quantities of heavy atoms to the specified length of the theoretical MIDs
            correction_matrix_dict[i] = new_formula.NaturalMID[0:n_theoretical_mid_entries]
            correction_row_normalization_factor = sum(correction_matrix_dict[i])
            correction_matrix_dict[i] = correction_matrix_dict[i]/correction_row_normalization_factor

        #make the correction matrix dictionary into a matrix
        CM = pandas.DataFrame(correction_matrix_dict)
        CM = pandas.DataFrame.transpose(CM)
        CM = np.asarray(CM)

        #find the right inverse (pseudo-inverse in numpy jargon) of the correction matrix
        CMi = np.linalg.pinv(CM)

        self.CM = CM
        self.CMi = CMi

    def calc_corrected_mid(self):

        import numpy as np

        #find the right inverse (pseudo-inverse in numpy jargon) of the correction matrix
        CMi = self.CMi

        #the theoretical MIDs are the rows of the correction matrix.
        #    Their length corresponds to the number of rows of the right inverse of the correcion matrix
        n_theoretical_mid_entries = CMi.shape[0]

        #the measured MID must have the same number of entries as each of the theoretical MIDs
        #    if it is short, add zeros to make up for the difference
        MIDu = self.MIDu
        if len(MIDu) < n_theoretical_mid_entries:
            mid_u_appendage = np.zeros(n_theoretical_mid_entries-len(MIDu))
            MIDu = np.append(MIDu,mid_u_appendage)

        #calculate corrected MID
        MIDc = np.dot(MIDu,CMi)
        MIDc = MIDc/sum(MIDc)
        self.assign('MIDc',MIDc)
