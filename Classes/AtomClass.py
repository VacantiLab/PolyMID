class Atom:

    # Initializer and Instance Attributes
    def __init__(self,AtomSymbol,AtomMID,AtomName):
        import numpy as np

        self.symbol = AtomSymbol
        self.MID = AtomMID
        self.name = AtomName

    # instance method to assign new values to an Atom object
    def assign(self,attribute,NewValue):
        if attribute == 'symbol':
            self.symbol = NewValue
        if attribute == 'MID':
            self.MID = NewValue
        if attribute == 'name':
            self.name = NewValue
