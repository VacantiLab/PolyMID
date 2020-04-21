class Input:
    # Initializer a Formula Instance and its attributes
    def __init__(self,fragment,TextFile):
        self.fragment = fragment
        self.TextFile = TextFile
        self.FragmentOrText = True #True if there are paratmers OR a text file provided, but not both
        self.AllRequired = True #True if all required parameters are provided directly or through the text file

    # Method for calculating natural MID of a Formula object
    def check_and_import(self):
        #Detemrine if the inputs were passed to the function properly
        #import fragment values if in a txt file

        from pdb import set_trace
        from Pesciolini import Fragment
        import numpy as np
        from Pesciolini import get_directory

        #If both inputs are None, then a GUI prompts for a .txt file input
        if ((self.fragment==None) & (self.TextFile==None)):
            self.TextFile = get_directory('gui_file')

        IsFragment = not (self.fragment is None)
        IsTextFile = not (self.TextFile is None)
        PassedCorrectly1 = not (IsFragment & IsTextFile)
        PassedCorrectly2 = IsFragment | IsTextFile
        self.FragmentOrText = PassedCorrectly1 & PassedCorrectly2

        # Get parameter values from text file if that is where they are provided
        if self.FragmentOrText & (not(self.TextFile is None)):
            # initialize variables
            FragmentName = None
            formula = None
            MetaboliteAtoms = None
            MIDu = None
            MIDc = None
            CM = None
            AtomLabeled = None

            #import values from text file
            with open(self.TextFile, 'r') as read_file:
                for line in read_file:
                    line_split = line.split(':')
                    line_split[0] = line_split[0].strip()
                    line_split[1] = line_split[1].strip()
                    if (line_split[0] == 'formula') | (line_split[0] == 'Formula'):
                        formula = line_split[1]
                    if (line_split[0] == 'MetaboliteAtoms') | (line_split[0] == 'Metabolite Atoms'):
                        MetaboliteAtoms = line_split[1]
                    if (line_split[0] == 'MIDu'):
                        MIDu = line_split[1]
                        MIDu = np.fromstring(MIDu,dtype=float,sep=' ')
                    if (line_split[0] == 'AtomLabeled') | (line_split[0] == 'Atom Labeled'):
                        AtomLabeled = line_split[1]
                    if (line_split[0] == 'FragmentName') | (line_split[0] == 'Fragment Name'):
                        FragmentName = line_split[1]
                    if (line_split[0] == 'CM'):
                        CM = line_split[1]
                        CM = None if CM == '' else CM

            self.fragment = Fragment(formula=formula, MetaboliteAtoms=MetaboliteAtoms, MIDu=MIDu, AtomLabeled=AtomLabeled, FragmentName=FragmentName, CM=CM, MIDc=None, PeakArea=None)
