class Input:
    #Initializer a Formula Instance and its attributes
    def __init__(self,fragment,TextFile,AtomLabeled):
        self.fragment = fragment
        self.TextFile = TextFile
        self.AtomLabeled = AtomLabeled
        self.FragmentOrText = True #True if there are paratmers OR a text file provided, but not both
        self.AllRequired = True #True if all required parameters are provided directly or through the text file

    #Method for checking input given correctly and importing values
    def check_and_import(self):
        #Detemrine if the inputs were passed to the function properly
        #    One of the arguemtns fragment or TextFile to correct() should be None. Both cannot have values.
        #import fragment values if in a txt file

        from pdb import set_trace
        from PolyMID import Fragment
        import numpy as np
        from PolyMID import get_directory
        from PolyMID import TextToCM

        #If both inputs are None, then a GUI prompts for a .txt file input
        if ((self.fragment==None) & (self.TextFile==None)):
            self.TextFile = get_directory('gui_file')

        IsFragment = not (self.fragment is None)
        IsTextFile = not (self.TextFile is None)
        PassedCorrectly1 = not (IsFragment & IsTextFile)
        PassedCorrectly2 = IsFragment | IsTextFile
        self.FragmentOrText = PassedCorrectly1 & PassedCorrectly2

        #Get parameter values from text file if that is where they are provided
        if self.FragmentOrText & IsTextFile:
            #Initialize variables
            FragmentName = None
            formula = None
            CanAcquireLabel = None
            MIDu = None
            MIDc = None
            #CM is initialized as a dictionary with None entries for atoms
            CM = {'C':None, 'H':None, 'N':None, 'O':None}

            #import values from text file
            with open(self.TextFile, 'r') as read_file:
                for line in read_file:
                    line_split = line.split(':')
                    line_split[0] = line_split[0].strip()
                    line_split[1] = line_split[1].strip()
                    if (line_split[0] == 'formula') | (line_split[0] == 'Formula'):
                        formula = line_split[1]
                    if (line_split[0] == 'CanAcquireLabel') | (line_split[0] == 'Metabolite Atoms'):
                        CanAcquireLabel = line_split[1]
                    if (line_split[0] == 'MIDu'):
                        MIDu = line_split[1]
                        MIDu = np.fromstring(MIDu,dtype=float,sep=' ')
                    if (line_split[0] == 'FragmentName') | (line_split[0] == 'Fragment Name'):
                        FragmentName = line_split[1]
                    if (line_split[0] == 'CM'):
                        if line_split[1] == 'None':
                            CM = CM
                        #If the text file contains information on CM
                        #    That information is added to the already initialized CM dictionary
                        if line_split[1] != 'None':
                            CM = TextToCM(line_split[1],CM)

            self.fragment = Fragment(formula=formula, CanAcquireLabel=CanAcquireLabel, MIDu=MIDu, FragmentName=FragmentName, CM=CM, MIDc=None, PeakArea=None)
