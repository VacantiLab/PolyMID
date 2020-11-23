#Run Command with text file input
#Run this from the directory containing the Pesciolini folder
import PolyMID
TextFile = '/Users/nate/git_hub_projects/PolyMID/references/ExampleInput.txt'
f = PolyMID.Correct(CorrectInput=TextFile)



FragmentName = 'TryptophanProtonated'
FragmentFormula = 'C11H13N2O2'
CanAcquireLabel = 'C11H13N2O2'
MIDm = [0.88885, 0.106829, 0.004322]
LabeledElement = 'C'
TracerEnrichment = 1
LabelEnrichment = 1
HighRes = ['N', 'O', 'H']
Input = PolyMID.Fragment(FragmentName, FragmentFormula, CanAcquireLabel, MIDm, LabeledElement, TracerEnrichment, LabelEnrichment, HighRes)
Output = PolyMID.Correct(Input)
