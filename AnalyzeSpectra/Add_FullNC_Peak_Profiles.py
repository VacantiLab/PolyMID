def Add_FullNC_Peak_Profiles(metabolite_dict,metabolite_list):
    from PolyMID import quantity_of_atom
    import numpy as np
    from pdb import set_trace

    print('this opened')

    FullNC_Metabolites = ['alanine_2tbdms','glycine_2tbdms','valine_2tbdms','leucine_2tbdms','isoleucine_2tbdms','serine_3tbdms', \
                          'threonine_3tbdms','methionine_2tbdms','phenylalanine_2tbdms','aspartate_3tbdms','proline_MOX_2tbdms', \
                          'glutamate_3tbdms','asparagine_3tbdms','lysine_3tbdms','glutamine_3tbdms','arginine_4tbdms','tyrosine_3tbdms', \
                          'histidine_3tbdms','tryptophan_3tbdms','glucose_clean']

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
                
                # Lengthen the mzs to integrate to include the tail peaks from the fully labeled (N and C) internal standard
                mzs_to_integrate_array = metabolite_dict[metabolite]['fragments'][key]['mzs_to_integrate']
                mzs_to_append = np.arange((mzs_to_integrate_array[-1]+1),(mzs_to_integrate_array[-1]+N_quantity+1),1)
                new_mzs_to_integrate = np.append(mzs_to_integrate_array,mzs_to_append)
                metabolite_dict[metabolite]['fragments'][key]['mzs_to_integrate'] = new_mzs_to_integrate

    return(metabolite_dict)
