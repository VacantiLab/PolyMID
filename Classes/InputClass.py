class Input:
    # Initializer a Formula Instance and its attributes
    def __init__(self,fragment,TextFile):
        self.fragment = fragment
        self.TextFile = TextFile
        self.FragmentOrText = True #True if there are paratmers OR a text file provided, but not both
        self.AllRequired = True #True if all required parameters are provided directly or through the text file

    # Method for calculating natural MID of a Formula object
    def check(self):
        #Detemrine if the inputs were passed to the function properly

        IsFragment = not (self.fragment is None)
        IsTextFile = not (self.TextFile is None)
        PassedCorrectly1 = not (IsFragment & IsTextFile)
        PassedCorrectly2 = IsFragment | IsTextFile
        self.FragmentOrText = PassedCorrectly1 & PassedCorrectly2

        # # Get parameter values from text file if that is where they are provided
        # if self.ParamsOrText & (self.TextFile not None):
        #     with open(TextFile, 'r') as read_file:
        #         for line in read_file:
        #             line_split = strip(line.split(':'))
        #             if (line_split[0] == 'formula') | (line_split[0] == 'Formula'):
        #                 self.formula = line_split[1]
        #             if (line_split[0] == 'MetaboliteAtoms') | (line_split[0] == 'Metabolite Atoms'):
        #                 self.MetaboliteAtoms = line_split[1]
        #             if (line_split[0] == 'MIDu'):
        #                 self.MIDu = line_split[1]
        #             if (line_split[0] == 'AtomLabeled') | (line_split[0] == 'Atom Labeled'):
        #                 self.AtomLabeled = line_split[1]
        #             if (line_split[0] == 'FragmentName') | (line_split[0] == 'Fragment Name'):
        #                 self.FragmentName = line_split[1]
        #             if (line_split[0] == 'CM'):
        #                 self.CM = line_split[1]



        # Check that all required parameters are entered
        #if self.ParamsOrText:
