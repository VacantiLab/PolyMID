def integrate_peaks(ic_smooth_dict,peak_start_t_dict,peak_end_t_dict,peak_start_i_dict,peak_end_i_dict,x_data_numpy,metabolite_dict,metabolite_list,ri_array,mz_vals,coelut_dict,coelut_dict_val,sample_name,Labeled_Element='C',Assume_All_MZs_Present=False):
    import importlib #allows fresh importing of modules
    from pdb import set_trace #python debugger
    import numpy as np #this is numpy, allows for data frame and matrix handling
    import pandas #a module which allows for making data frames
    import scipy #contains simpsons integration function
    import copy
    import PolyMID
    from PolyMID.AnalyzeSpectra import fragment_library #a custom function
    from PolyMID.AnalyzeSpectra import ri_to_rt
    from PolyMID.AnalyzeSpectra import quantity_of_atom
    from PolyMID.AnalyzeSpectra import match_fingerprint

    #initialize a dictionary that will contain all of the quantified metabolite information
    metabolite_dict_complete = copy.deepcopy(metabolite_dict)

    #store the scan acquisition times in a variable with a more sensible name (rename the original variable in source code later)
    sat_array = x_data_numpy

    #iterate through the metabolite names so you can then iterate through the fragments of each metabolite for integration
    for metabolite_iter in metabolite_list:
        fragments_list = list(dict.keys(metabolite_dict[metabolite_iter]['fragments']))
        ri_window = 3
        met_present,ri = match_fingerprint.match_fingerprint(ri_array,coelut_dict,coelut_dict_val,metabolite_dict,mz_vals,ic_smooth_dict,metabolite_iter,sample_name,ri_window)
        if met_present:
            rt = ri_to_rt.ri_to_rt(sat_array,ri_array,ri)
        if not met_present:
            ri = 0
            rt = 0
        metabolite_dict_complete[metabolite_iter]['ri'] = ri
        metabolite_dict_complete[metabolite_iter]['rt'] = rt

        #iterate through the fragments of each metabolite and integrate
        for frag_iter in fragments_list:
            mzs_to_integrate = metabolite_dict_complete[metabolite_iter]['fragments'][frag_iter]['mzs_to_integrate']
            M0 = True # indicates the first mz indicated corresponds to M0
            for i in mzs_to_integrate:
                #find the index of the peak corresponding to the given rt in the mz curve (the first peak is 0, the second 1, the third 2, ...)
                peak_present = False
                if i in mz_vals: #it is possible there were no values above threshhold recorded in the ms scan so there would be no entry for that mz in the data dictionary
                    if met_present: #the presence of the metabolite is checked here instead of outside the fragment iteration loop because the fragments are iterated through to fill with zeros for printing even if the metabolite is not present.
                                    #this seems inefficient and can be fixed later - integration is not really expensive at this point
                        if Assume_All_MZs_Present:
                            peak_present = True
                        if i in coelut_dict[ri]: #recall coelut dict has keys of ri and for each key is an array of mz's with peaks eluting at that ri
                            peak_present = True

                #if there is no peak detected for that mz_to_integrate of the fragment, record zero as the peak area
                if not peak_present:
                    metabolite_dict_complete[metabolite_iter]['fragments'][frag_iter]['areas'] = np.append(metabolite_dict_complete[metabolite_iter]['fragments'][frag_iter]['areas'],0)

                #if there is a peak, then integrate it to find its area
                if peak_present:
                    possible_peak_starts = np.where(peak_start_t_dict[i] < rt)[0] #there must be at least one peak start before the rt because there is a peak
                    prosp_peak_start_nm = max(possible_peak_starts) #prospective peak start
                    x_start_i = peak_start_i_dict[i][prosp_peak_start_nm]
                    x_end_i = peak_end_i_dict[i][prosp_peak_start_nm]
                    x_range_i = np.arange(x_start_i,x_end_i)
                    x_range_t = x_data_numpy[x_range_i]
                    y_range_ic = ic_smooth_dict[i][x_range_i]

                    # can specify certain fragments will be integrated only as the left or right side of the peak
                    #     hard-coded for now
                    Decon_Left_List = ['glucose-QuantOnly_273']
                    Decon_Right_List = []

                    if frag_iter in Decon_Left_List:
                        # redefine the right border of the peak as the midway of the M0 peak
                        if M0:
                            domain_length = np.floor(0.5*len(x_range_t)).astype(int)
                            max_peak_time = x_range_t[domain_length-1]
                        x_range_t = x_range_t[x_range_t<max_peak_time]
                        y_range_ic = y_range_ic[0:len(x_range_t)]

                    if frag_iter in Decon_Right_List:
                        if M0:
                        # redefine the left border of the peak as the midway of the M0 peak
                            domain_length = np.floor(0.5*len(x_range_t)).astype(int)
                            min_peak_time = x_range_t[domain_length-1]
                        x_range_t = x_range_t[x_range_t>min_peak_time]
                        y_range_ic = y_range_ic[len(y_range_ic)-len(x_range_t):-1]

                    integrated_area = scipy.integrate.simps(y_range_ic,x_range_t)
                    metabolite_dict_complete[metabolite_iter]['fragments'][frag_iter]['areas'] = np.append(metabolite_dict_complete[metabolite_iter]['fragments'][frag_iter]['areas'],integrated_area)
                    M0 = False # indicates the next mz in the loop is not M0

            #record the total area for the fragment by summing the areas of all mz members of that fragment
            metabolite_dict_complete[metabolite_iter]['fragments'][frag_iter]['tot_area'] = np.sum(metabolite_dict_complete[metabolite_iter]['fragments'][frag_iter]['areas'])

            #if peaks are found, the mid must be normalized and corrected for natural abundances
            if metabolite_dict_complete[metabolite_iter]['fragments'][frag_iter]['tot_area'] > 0:
                metabolite_dict_complete[metabolite_iter]['fragments'][frag_iter]['mid'] = metabolite_dict_complete[metabolite_iter]['fragments'][frag_iter]['areas']/metabolite_dict_complete[metabolite_iter]['fragments'][frag_iter]['tot_area']
                #correct the mids, will work if the MID is all zeros
                #    in order to print, all fragments must have a corrected mid
                #    these corrected MIDs must be the same length for a fragment across all samples (are formula dependent so they will be)

                fragment = PolyMID.Fragment(FragmentName = frag_iter, \
                                            FragmentFormula = metabolite_dict_complete[metabolite_iter]['fragments'][frag_iter]['formula'], \
                                            CanAcquireLabel = metabolite_dict_complete[metabolite_iter]['fragments'][frag_iter]['metabolite_atoms'], \
                                            MIDm = metabolite_dict_complete[metabolite_iter]['fragments'][frag_iter]['mid'], \
                                            LabeledElement = Labeled_Element, \
                                            TracerEnrichment = 1, \
                                            LabelEnrichment = 1, \
                                            HighRes = 'none', \
                                            CM = metabolite_dict_complete[metabolite_iter]['fragments'][frag_iter]['CM'])
                fragment = PolyMID.Correct(fragment)
                metabolite_dict_complete[metabolite_iter]['fragments'][frag_iter]['mid_c'] = fragment.MIDc

            #if there are no peaks found, the metabolite is not present and all mid entries are 0
            if metabolite_dict_complete[metabolite_iter]['fragments'][frag_iter]['tot_area'] == 0:
                metabolite_atoms = metabolite_dict_complete[metabolite_iter]['fragments'][frag_iter]['metabolite_atoms']
                atom_labeled = 'C'
                atom_quantity = quantity_of_atom.quantity_of_atom(metabolite_atoms,atom_labeled)
                n_mid_entries = atom_quantity + 1 #accounts for M0
                metabolite_dict_complete[metabolite_iter]['fragments'][frag_iter]['mid'] = metabolite_dict_complete[metabolite_iter]['fragments'][frag_iter]['areas']
                metabolite_dict_complete[metabolite_iter]['fragments'][frag_iter]['mid_c'] = np.zeros(n_mid_entries)


            #clever way to iterate through dictionary - not used here anymore
            #if fragment_dict_complete[frag_iter]['tot_area'] == 0:
            #    fragment_dict_complete[frag_iter]['mid'] = {k: 0 for k, v in fragment_dict_complete[frag_iter]['areas'].items()}

    return(metabolite_dict_complete)
