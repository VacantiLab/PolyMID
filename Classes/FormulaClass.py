class Formula:
    # Initializer a Formula Instance and its attributes
    def __init__(self,formula):
        self.FormulaInput = formula
        self.formula = None
        self.AtomArray = None
        self.NaturalMID = None
        self.FormatFormula()
        self.CreateAtomArray()

    def FormatFormula(self):
        # This method redefines the formula attribute of the Formula class
        #     This method adds 1s (ones) where they are implied to the input string representing the chemical formula (self.FormulaInput)

        import re
        from pdb import set_trace

        # Find the flanking characters where a 1 should be inserted (1s that are internal to the formula)
        FormulaToEdit = self.FormulaInput
        missing1s = re.findall('(?:[A-Z]|[a-z])[A-Z]',FormulaToEdit)
        #    '(?:)' in '(?:RegExpHere)' denotes a non-capturing group, i.e. it is used for specifying order of operations as opposed to just '()'

        # Insert the internal 1s
        for entry in missing1s:
            NewEntry = entry[0] + '1' + entry[1]
            FormulaToEdit = re.sub(entry,NewEntry,FormulaToEdit)

        # Append an external 1 if necessary
        #     If there is a letter at the end of the string, replace it with that letter and a 1
        LetterAtEnd = re.findall('(?:[A-Z]|[a-z])$',FormulaToEdit)
        #    '(?:)' in '(?:RegExpHere)' denotes a non-capturing group, i.e. it is used for specifying order of operations as opposed to just '()'
        if (len(LetterAtEnd) > 0):
            FormulaToEdit = re.sub(LetterAtEnd[0]+'$',LetterAtEnd[0]+'1',FormulaToEdit)

        # Redefine the formula attribute of the Formula obeject
        self.formula = FormulaToEdit

    def CreateAtomArray(self):
        # This method creates an array of atom objects corresponding to the formula

        # Import required modules
        from pdb import set_trace
        import numpy as np
        import re
        from PolyMID import Atom

        # Break the fragment formula up into its atomic symbol and letter components
        broken_formula = np.array(re.findall('[A-Z][a-z]?|[0-9]+', self.formula))
        #    '[A-Z][a-z]?|[0-9]+': A capital letter [A-Z], followed by a lowercase letter [a-z], which is optional '?', or '|' a number '[0-9]', and possibly more numbers '+'
        #    example: this command will take formula = C6H12N1O3Si1 and return broken_formula = array(['C','6','H','12','N','1','O','3','Si','1'])
        #        all components are strings

        # Get separate arrays containing the atomic symbols and corresponding formula numbers
        odd_indices = np.array(range(1,len(broken_formula),2)) #all of the odd indices of the array broken_formula
        even_indices = np.array(range(0,len(broken_formula),2)) #all of the even indices of the array broken_formula
        formula_numbers = broken_formula[odd_indices].astype(np.int) #the atomic numbers array
        formula_atoms = broken_formula[even_indices] #the atomic symbols array

        AtomArray = []
        AtomCounter = 0
        for atom in formula_atoms:
            AtomStoich = formula_numbers[AtomCounter]
            AtomArray.append(Atom(atom,AtomStoich))
            AtomCounter = AtomCounter + 1

        self.AtomArray = AtomArray


    # Method for calculating natural MID of a Formula object
    def calc_natural_mid(self):
        # This method defines the NatrualMID attribute of the Formula class
        # self.formula: the formula for the fragment whose natural isotopic abundance is being calculated
        #    each atomic symbol must start with a capital letter, if there is a second letter it must be lowercase
        #    if there is no number following an atomic symbol, that number is assumed to be 1

        # Import required modules
        from pdb import set_trace
        import numpy as np
        import re
        from PolyMID import expand_polynomial

        # Define the atomic isotopic abundances
        atom_abundances = dict()
        atom_abundances['C'] = np.array([0.9893,0.0107,0])
        atom_abundances['H'] = np.array([0.999885,0.000115,0])
        atom_abundances['N'] = np.array([0.99632,0.00368,0])
        atom_abundances['O'] = np.array([0.99757,0.00038,0.00205])
        atom_abundances['Si'] = np.array([0.922297,0.046832,0.030872])
        atom_abundances['S'] = np.array([0.9500,0.0075,0.0425])
        atom_abundances['P'] = np.array([1,0,0])
        atom_abundances['Hv'] = np.array([0,1,0])

        # Break the fragment formula up into its atomic symbol and letter components
        broken_formula = np.array(re.findall('[A-Z][a-z]?|[0-9]+', self.formula))
        #    '[A-Z][a-z]?|[0-9]+': A capital letter [A-Z], followed by a lowercase letter [a-z], which is optional '?', or '|' a number '[0-9]', and possibly more numbers '+'
        #    example: this command will take formula = C6H12N1O3Si1 and return broken_formula = array(['C','6','H','12','N','1','O','3','Si','1'])
        #        all components are strings

        # Get separate arrays containing the atomic symbols and corresponding formula numbers
        odd_indices = np.array(range(1,len(broken_formula),2)) #all of the odd indices of the array broken_formula
        even_indices = np.array(range(0,len(broken_formula),2)) #all of the even indices of the array broken_formula
        formula_numbers = broken_formula[odd_indices].astype(np.int) #the atomic numbers array
        formula_atoms = broken_formula[even_indices] #the atomic symbols array

        # Create a dictionary containing all of the the atom mids
        #     Each key corresponds to an atom in the fragment
        #     If an atom appears more than once in the fragment, the dictionary will contain idenical entries for they keys
        #     The mass isotopomer distributions are stored as the entries for each keys
        #     Example:
        #                 carbon1:   M0 M1 M2 M3 M4...
        #                 carbon2:   M0 M1 M2 M3 M4...
        #                                    .
        #                                    .
        #                 hydrogen1: M0 M1 M2 M3 M4 M5...
        #                 hydrogen2: M0 M1 M2 M3 M4 M5...
        #                                        .
        #                                        .

        # Initialize the dictionary containing the atom MIDs
        atom_mids_dict = {}

        # Determine the number of atoms in the fragment (not unique atoms)
        n_atoms = np.sum(formula_numbers)

        # Initialize the atom counter, this tracks the identity of the atom
        #     If you go from the first to the second C, this will not changed
        #     If you go from the last C to the first H this will increase by one
        atom_index = 0

        # Populate the atom_mid_dict
        #     Iterate through the atoms in the fragment
        for i in range(0,n_atoms): #iterate through the atoms in the fragment

            # Find the next row that the identity of the atom will change
            #     '+1' is necessary because the end of the range is not included in python
            new_atom_id_index = sum(formula_numbers[0:atom_index+1])

            if i < new_atom_id_index: #if the atom identity has not changed
                current_atom = formula_atoms[atom_index]
                atom_mids_dict[current_atom + str(i)] = atom_abundances[current_atom]
            if i >= new_atom_id_index: #if the atom identity has changed
                current_atom = formula_atoms[atom_index+1]
                atom_mids_dict[current_atom + str(i)] = atom_abundances[current_atom]
                atom_index = atom_index + 1 #only iterate the atom_index if the atom identity has changed

        # Given the formula A2B3 where A and B have isotopic abundances of [A_M0 A_M1 A_M2] and [B_M0 B_M1 B_M2]
        #     The vector of mass isotopomeric abundace of A2B3 can be found by multiplying out the vecotors of atomic isotopic abundances as if they were polynomials, and combining like terms
        #     The M0, M1, ... are be treated as exponents; thus (A_M0)*(B_M2) gives an M2 term (0+2)

        # Initialize the MID vector of the fragemnt specified by the formula
        expanded_placeholder = np.array([1])

        # Perform the polynomial expansion; expanding and collecting like terms for two MIDs (atoms) at a time
        for key in atom_mids_dict.keys():
            expanded = expand_polynomial(expanded_placeholder,atom_mids_dict[key])
            expanded_placeholder = np.array(expanded['prob'])

        # Define the NaturalMID attribute of the Formula object
        natural_mid = np.array(expanded['prob'])
        self.NaturalMID = natural_mid
