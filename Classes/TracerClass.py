class Tracer:

    # Initializer and Instance Attributes
    def __init__(self,AtomLabeled,TracerEnrichment,LabelEnrichment):
        from pdb import set_trace

        self.AtomLabeled = AtomLabeled # A string indicating the atom that is considered to be labeled
        self.TracerEnrichment = TracerEnrichment
        self.LabelEnrichment = LabelEnrichment
