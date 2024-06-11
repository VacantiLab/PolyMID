def print_integrated_peaks(file_directory,samples,metabolite_list,file_data,corrected,Full_NC):

    import numpy as np
    from pdb import set_trace

    CorrectedOptionsArray = [True,False]

    for corrected in CorrectedOptionsArray:

        #specify if the corrected or uncorrected (for natural isotopic abundances) MIDs are printed
        if corrected:
            mid_to_print = 'mid_c'
            text_file_name = 'IntegratedPeaks_CorrectedNaturalAbundances.txt'

        if not corrected:
            mid_to_print = 'mid'
            text_file_name = 'IntegratedPeaks_UncorrectedNaturalAbundances.txt'

        output_text_file = file_directory + text_file_name
        file_object_text = open(output_text_file,'w')

        samples_string = '\t'.join(samples)
        file_object_text.write(' \t' + samples_string)
        file_object_text.write('\n')

        file_object_text.write('Retention Indices\n')
        for metabolite_iter in metabolite_list:
            fragment_list = list(dict.keys(file_data[samples[0]]['metabolites'][metabolite_iter]['fragments']))
            for fragment in fragment_list:
                file_object_text.write(fragment)
                file_object_text.write('\t')
                for sample_name in samples:
                    ri = file_data[sample_name]['metabolites'][metabolite_iter]['ri']
                    ri = np.round(ri,decimals=1)
                    ri = str(ri)
                    file_object_text.write(ri)
                    file_object_text.write('\t')
                file_object_text.write('\n')

        file_object_text.write('\n')
        file_object_text.write('Peak Areas (ion counts * seconds)\n')
        for metabolite_iter in metabolite_list:
            fragment_list = list(dict.keys(file_data[samples[0]]['metabolites'][metabolite_iter]['fragments']))
            for fragment in fragment_list:
                file_object_text.write(fragment)
                file_object_text.write('\t')
                for sample_name in samples:
                    area_to_round = file_data[sample_name]['metabolites'][metabolite_iter]['fragments'][fragment]['tot_area']
                    area_to_print = np.round(area_to_round,decimals=0)
                    area = str(area_to_print)
                    file_object_text.write(area)
                    file_object_text.write('\t')
                file_object_text.write('\n')

        file_object_text.write('\n')
        file_object_text.write('MIDs\n')
        for metabolite_iter in metabolite_list:
            fragment_list = list(dict.keys(file_data[samples[0]]['metabolites'][metabolite_iter]['fragments']))
            for fragment in fragment_list:
                mid_length = len(file_data[samples[0]]['metabolites'][metabolite_iter]['fragments'][fragment][mid_to_print])
                mid_members = range(0,mid_length)
                for M in mid_members:
                    fragment_mi_name = fragment + ' ' + 'M' + str(mid_members[M])
                    file_object_text.write(fragment_mi_name)
                    file_object_text.write('\t')
                    for sample_name in samples:
                        mi_to_round = file_data[sample_name]['metabolites'][metabolite_iter]['fragments'][fragment][mid_to_print][M]
                        mi_to_print = np.round(mi_to_round,decimals=6)
                        current_mi = str(mi_to_print)
                        file_object_text.write(current_mi)
                        file_object_text.write('\t')
                    #print the MID for the fragment as it exists in nature
                    mi_to_round = file_data[sample_name]['metabolites'][metabolite_iter]['fragments'][fragment]['natural_mid'][M]
                    mi_to_print = np.round(mi_to_round,decimals=6)
                    current_mi = str(mi_to_print)
                    file_object_text.write(current_mi)
                    file_object_text.write('\t')

                    file_object_text.write('\n')
                file_object_text.write('\n')

        file_object_text.close()
    return()
