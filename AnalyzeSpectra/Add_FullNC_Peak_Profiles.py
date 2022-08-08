def Add_FullNC_Peak_Profiles(metabolite_dict,metabolite_list):
    from PolyMID import quantity_of_atom
    from pdb import set_trace

    FullNC_Metabolites = ['alanine_2tbdms','glycine_2tbdms','valine_2tbdms','leucine_2tbdms','isoleucine_2tbdms','serine_3tbdms', \
                          'threonine_3tbdms','methionine_2tbdms','phenylalanine_2tbdms','aspartate_3tbdms','proline_MOX_2tbdms', \
                          'glutamate_3tbdms','asparagine_3tbdms','lysine_3tbdms','glutamine_3tbdms','arginine_4tbdms','tyrosine_3tbdms', \
                          'histidine_3tbdms','tryptophan_3tbdms']

    for metabolite in metabolite_list:
        if metabolite in FullNC_Metabolites:
            Fragment_Keys = metabolite_dict[metabolite]['fragments'].keys()
            for key in Fragment_Keys:
                Metabolite_Atoms = metabolite_dict[metabolite]['fragments'][key]['metabolite_atoms']
                Fragment_Parent_MZ = metabolite_dict[metabolite]['fragments'][key]['mzs_to_integrate'][0]
                N_quantity = quantity_of_atom(Metabolite_Atoms,'N')
                C_quantity = quantity_of_atom(Metabolite_Atoms,'C')
                Adjustment_to_Peak_Profile = N_quantity + C_quantity
                Peak_Profile = metabolite_dict[metabolite]['peak_profile']
                Indices_to_Change = Peak_Profile - Fragment_Parent_MZ == 0
                Peak_Profile[Indices_to_Change] = Peak_Profile[Indices_to_Change] + Adjustment_to_Peak_Profile
                metabolite_dict[metabolite]['peak_profile'] = Peak_Profile

    return(metabolite_dict)
