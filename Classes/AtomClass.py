class Atom:

    # Initializer and Instance Attributes
    def __init__(self,AtomSymbol,AtomStoich):
        import numpy as np
        from pdb import set_trace

        self.symbol = AtomSymbol
        self.stoich = AtomStoich
        self.MID = None

        # Retrieve the Atom MID from the text file PolyMID/SupportingFiles/AtomMIDs.txt
        self.ReadMID()

    # instance method to assign new values to an Atom object
    def assign(self,attribute,NewValue):
        if attribute == 'symbol':
            self.symbol = NewValue
        if attribute == 'MID':
            self.MID = NewValue
        if attribute == 'name':
            self.name = NewValue

    def ReadMID(self):
        # This method opens the text file defining the AtomMIDs and imports that information into the attribute MID as a numpy array

        import os.path
        import numpy as np
        import PolyMID
        from pdb import set_trace

        PolyMID_Path = os.path.abspath(PolyMID.__file__)
        PolyMID_Path = PolyMID_Path.split(sep='/')
        PolyMID_Path = PolyMID_Path[:-1]
        AtomMIDs_txtPath = '/'.join(PolyMID_Path) + '/SupportingFiles/AtomMIDs.txt'

        with open(AtomMIDs_txtPath,'r') as AtomMIDsFile:
            for line in AtomMIDsFile:
                line_split = line.split(':')
                FileAtomSymbol = line_split[0].strip()
                FileAtomMID_String = line_split[1].strip()

                if FileAtomSymbol == self.symbol:
                    AtomMID_String = FileAtomMID_String
                    AtomMID_StringArray = AtomMID_String.split(sep=' ')
                    AtomMID_FloatArray = [float(i) for i in AtomMID_StringArray]
                    self.MID = np.asarray(AtomMID_FloatArray)
