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

        # Get parameter values from text file if that is where they are provided
        #if self.ParamsOrText & (self.TextFile not None):


        # Check that all required parameters are entered
        #if self.ParamsOrText:
