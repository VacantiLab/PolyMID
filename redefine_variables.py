pyruvate174 = Fragment('C6H12N1O3Si1','C3H3O2',np.array([0.854,0.104,0.0388,0.00280]),np.array([]),np.array([]),np.array([]))

from Pesciolini import PolyMID
from Pesciolini import Fragment
import numpy as np
formula='C6H12N1O3Si1'
MetaboliteAtoms = 'C3H3O2'
MIDu = np.array([0.854,0.104,0.0388,0.00280])
FragmentName = 'pyruvate174'
CM = None
AtomLabeled = 'C'
MIDc = None
PeakArea = None
fragment = Fragment(formula=formula,MetaboliteAtoms=MetaboliteAtoms,MIDu=MIDu,CM=CM,MIDc=MIDc,PeakArea=PeakArea,AtomLabeled=AtomLabeled,FragmentName=FragmentName)

#Run Command with text file input
#Run this from the directory containing the Pesciolini folder
from Pesciolini import PolyMID
TextFile = '/Users/nate/git_hub_projects/Pesciolini/PolyMID/ExampleInputWithoutCM.txt'
PolyMID.main(TextFile=TextFile)


fragment.create_correction_matrix()
