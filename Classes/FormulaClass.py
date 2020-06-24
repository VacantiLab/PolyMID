class Formula:
    # Initializer a Formula Instance and its attributes
    def __init__(self,formula):
        self.FormulaInput = formula
        self.formula = None
        self.NaturalMID = None
        self.FormatFormula()

    def FormatFormula(self):
        # This method adds 1s (ones) where they are implied to the input string representing the chemical formula
        # It then redefines the formula attribute of the Formula object

        import re
        from pdb import set_trace

        #Find the flanking characters where a 1 should be inserted (1's that are internal to the formula)
        FormulaToEdit = self.FormulaInput
        missing1s = re.findall('(?:[A-Z]|[a-z])[A-Z]',FormulaToEdit)
        #    '(?:)' in '(?:RegExpHere)' denotes a non-capturing group, i.e. it is used for specifying order of operations as opposed to just '()'

        #Insert the 1s
        for entry in missing1s:
            NewEntry = entry[0] + '1' + entry[1]
            FormulaToEdit = re.sub(entry,NewEntry,FormulaToEdit)

        #If there is a letter at the end of the string, replace it with that letter and a 1
        LetterAtEnd = re.findall('(?:[A-Z]|[a-z])$',FormulaToEdit)
        #    '(?:)' in '(?:RegExpHere)' denotes a non-capturing group, i.e. it is used for specifying order of operations as opposed to just '()'
        if (len(LetterAtEnd) > 0):
            FormulaToEdit = re.sub(LetterAtEnd[0]+'$',LetterAtEnd[0]+'1',FormulaToEdit)

        self.formula = FormulaToEdit

    # Method for calculating natural MID of a Formula object
    def calc_natural_mid(self):
        #formula: the formula for the fragment whose natural isotopic abundance is being calculated
        #    each atomic symbol must start with a capital letter, if there is a second letter it must be lowercase
        #    a number needs to follow each atomic symbol, even if that number is 1

        #import required modules
        from pdb import set_trace
        import numpy as np #this is numpy
        import pandas
        import re
        from Pesciolini import expand_polynomial

        #define the atomic isotopic abundances, there must be an equal number for each atom
        atom_abundances = dict()
        atom_abundances['C'] = pandas.Series([0.9893,0.0107,0],index=['M0','M1','M2'])
        atom_abundances['H'] = pandas.Series([0.999885,0.000115,0],index=['M0','M1','M2'])
        atom_abundances['N'] = pandas.Series([0.99632,0.00368,0],index=['M0','M1','M2'])
        atom_abundances['O'] = pandas.Series([0.99757,0.00038,0.00205],index=['M0','M1','M2'])
        atom_abundances['Si'] = pandas.Series([0.922297,0.046832,0.030872],index=['M0','M1','M2'])
        atom_abundances['S'] = pandas.Series([0.9500,0.0075,0.0425],index=['M0','M1','M2'])
        atom_abundances['P'] = pandas.Series([1,0,0],index=['M0','M1','M2'])
        atom_abundances['Hv'] = pandas.Series([0,1,0],index=['M0','M1','M2'])

        #break the fragment formula up into its atomic symbol and letter components
        broken_formula = np.array(re.findall('[A-Z][a-z]?|[0-9]+', self.formula))
        #    '[A-Z][a-z]?|[0-9]+': A capital letter [A-Z], followed by a lowercase letter [a-z], which is optional '?', or '|' a number '[0-9]', and possibly more numbers '+'
        #    example: this command will take formula = C6H12N1O3Si1 and return broken_formula = array(['C','6','H','12','N','1','O','3','Si','1'])
        #        all components are strings

        #get separate arrays containing the atomic symbols and corresponding formula numbers
        odd_indices = np.array(range(1,len(broken_formula),2)) #all of the odd indices of the array broken_formula
        even_indices = np.array(range(0,len(broken_formula),2)) #all of the even indices of the array broken_formula
        formula_numbers = broken_formula[odd_indices].astype(np.int) #the atomic numbers array
        formula_atoms = broken_formula[even_indices] #the atomic symbols array

        #create a dictionary containing all of the the atom mids
        #    each key corresponds to an atom in the fragment
        #    if an atom appears more than once in the fragment, the dictionary will contain idenical entries for they keys
        #    the mass isotopomer distributions are stored as the entries for each keys
        #    example:
        #                 carbon1:   M0 M1 M2 M3 M4...
        #                 carbon2:   M0 M1 M2 M3 M4...
        #                                    .
        #                                    .
        #                 hydrogen1: M0 M1 M2 M3 M4 M5...
        #                 hydrogen2: M0 M1 M2 M3 M4 M5...
        #                                        .
        #                                        .

        #Initialize the dictionary containing the atom MIDs
        atom_mids_dict = {}

        #determine the number of atoms in the fragment (not unique atoms)
        n_atoms = np.sum(formula_numbers)

        #initialize the atom counter, this tracks the identity of the atom
        #    if you go from the first to the second C, this will not changed
        #    if you go from the last C to the first H this will increase by one
        atom_index = 0

        #Populate the atom_mid_dict
        #    iterate through the atoms in the fragment
        for i in range(0,n_atoms): #iterate through the atoms in the fragment

            #Find the next row that the identity of the atom will change
            #    '+1' is necessary because the end of the range is not included in python
            new_atom_id_index = sum(formula_numbers[0:atom_index+1])

            if i < new_atom_id_index: #if the atom identity has not changed
                current_atom = formula_atoms[atom_index]
                atom_mids_dict[current_atom + str(i)] = atom_abundances[current_atom]
            if i >= new_atom_id_index: #if the atom identity has changed
                current_atom = formula_atoms[atom_index+1]
                atom_mids_dict[current_atom + str(i)] = atom_abundances[current_atom]
                atom_index = atom_index + 1 #only iterate the atom_index if the atom identity has changed

        #given the formula A2B3 where A and B have isotopic abundances of [A_M0 A_M1 A_M2] and [B_M0 B_M1 B_M2]
        #the vector of mass isotopomeric abundace of A2B3 can be found by multiplying out the vecotors of atomic isotopic abundances as if they were polynomials, and combining like terms
        #the M0, M1, ... are be treated as exponents; thus (A_M0)*(B_M2) gives an M2 term (0+2)
        expanded_placeholder = [1] #initializes the MID vector of the fragemnt specified by the formula
        for key in atom_mids_dict.keys():
            expanded = expand_polynomial(expanded_placeholder,atom_mids_dict[key])
            expanded_placeholder = np.array(expanded['prob'])

        natural_mid = np.array(expanded.iloc[:,0])
        self.NaturalMID = natural_mid
