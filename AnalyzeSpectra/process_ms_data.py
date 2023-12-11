def process_ms_data(sat,ic_df,output_plot_directory,n_scns,mz_vals,low_sensitivity):

    import importlib #allows fresh importing of modules
    from pdb import set_trace #python debugger
    from PolyMID.AnalyzeSpectra import savitzky_golay
    import bokeh.plotting as bkp #allows for making interactive plots with bokeh
    import numpy as np #this is numpy, allows for data frame and matrix handling
    import pandas #a module which allows for making data frames
    import peakutils #allows for peak and baseline finding
    from PolyMID.AnalyzeSpectra import FindBorders
    from PolyMID.AnalyzeSpectra import ExtendBounds
    import copy
    from scipy import signal

    mz_to_plot = [174,175,176,177,178,179]
    mz_colors = ['red','green','blue','yellow','purple','orange']
    #bkp.output_file(output_plot_directory)
    p=bkp.figure(title='ic vs. time', x_axis_label='time',y_axis_label='ic',plot_width=1000)
    ic_smooth_dict = dict()
    ic_nsmooth_dict = dict()
    ic_dict = dict()
    peak_start_t_dict = dict()
    peak_end_t_dict = dict()
    peak_start_i_dict = dict()
    peak_end_i_dict = dict()
    peak_i_dict = dict()
    peak_max_dict = dict()
    peak_sat_dict = dict()

    j=0
    for plotted_mz in mz_vals:
        #arrange data to be plotted
        x_data = sat
        y_data = ic_df.loc[plotted_mz,:]

        #convert to numpy objects
        x_data_numpy = np.array(x_data)
        y_data_numpy = y_data.to_numpy()

        #smooth the data
        y_data_smooth = savitzky_golay.savitzky_golay(y_data_numpy, window_size=7, order=2, deriv=0, rate=1)

        #if values go negative because of smoothing, replace them with the original values
        smooth_negative = y_data_smooth < 0
        smooth_neg_indices = np.where(smooth_negative)[0] #np.where returns a tuple, the first item gives the locations of True
        n_y = len(y_data_numpy)

        smooth_neg_indices = ExtendBounds.ExtendBounds(smooth_neg_indices,n_y) #extend out two from each negative interval for better continuity
        y_data_smooth[smooth_neg_indices] = y_data_numpy[smooth_neg_indices]

        #find the indices of peaks
        # set the peak threshold parameter to pass to peakutils.indexes
        thres_numerator = 500
        if low_sensitivity:
            thres_numerator = 20
        thres = thres_numerator/np.amax(y_data_smooth)
        if thres > 1:
            thres = 1

        # Find the indices of peaks
        indexes = signal.find_peaks(x=y_data_smooth,height=thres_numerator)[0]
        #     Enforcing a prominence here causes missed peaks

        # Find the intensities associated with the peaks
        peak_values = y_data_smooth[indexes]

        # Find the peak prominences and relative prominences
        prominences = signal.peak_prominences(y_data_smooth, indexes, wlen=None)[0]
        relative_prominence = prominences/peak_values

        # enforce a relative peak prominence threshold
        indices_of_peak_indices_to_keep = relative_prominence = prominences/peak_values >= 0.5
        indexes = indexes[indices_of_peak_indices_to_keep]

        #a peak must be higher than the max radius points before and after it
        max_radius = 5
        delete_indices = np.array([], dtype=int)
        index_counter = 0
        for index in indexes:
            point_to_test = y_data_smooth[index]
            if (index >= max_radius) & (index <= n_scns-max_radius): #if there is a radius length around the point (it's not too close to the end)
                neighbor_indices = np.arange(index-max_radius,index+max_radius,1)
                keep_index = all(point_to_test >= y_data_smooth[entry] for entry in neighbor_indices)
                if not keep_index:
                    delete_indices = np.append(delete_indices, int(index_counter))
            index_counter = index_counter+1
        indexes = np.delete(indexes,delete_indices)

        #find the baseline
        base = peakutils.baseline(y_data_smooth,deg=5)

        #if the baseline is found to be negative in places, replace it with 0 in those places
        #    a negative baseline can raise peaks interfering with border detection and integration
        base_negative = base < 0
        base_neg_indices = np.where(base_negative)[0] #np.where returns a tuple, the first item gives the locations of True
        base[base_neg_indices] = 0

        #if the baseline is found to be higher than the measured value, replace it with the measured value
        #    artifacts of fitting a polynomial that cause the baseline to rise above the measurement are not taken seriously
        base_too_high = base > y_data_smooth
        base_too_high_indices = np.where(base_too_high)[0]
        base[base_too_high_indices] = y_data_smooth[base_too_high_indices]
        y_smooth_cor = y_data_smooth - base
        ic_smooth_dict[plotted_mz] = y_smooth_cor
        ic_dict[plotted_mz] = y_data_numpy
        y_nsmooth_cor = y_data_numpy - base
        ic_nsmooth_dict[plotted_mz] = y_nsmooth_cor

        #get the x and y peak locations and values resepctively.
        x_peak_locs = x_data_numpy[indexes]
        y_peak_vals_base_cor = y_smooth_cor[indexes]

        #Find the peak borders
        borders_prior_indices,borders_after_indices = FindBorders.FindBorders(indexes,y_smooth_cor,sat) #user defined function
        borders_prior_sat = x_data_numpy[borders_prior_indices]
        borders_prior_values = y_smooth_cor[borders_prior_indices]
        borders_after_sat = x_data_numpy[borders_after_indices]
        borders_after_values = y_smooth_cor[borders_after_indices]

        peak_start_t_dict[plotted_mz] = borders_prior_sat
        peak_end_t_dict[plotted_mz] = borders_after_sat
        peak_start_i_dict[plotted_mz] = borders_prior_indices
        peak_end_i_dict[plotted_mz] = borders_after_indices
        peak_i_dict[plotted_mz] = indexes
        peak_sat_dict = sat[indexes]
        peak_max_dict[plotted_mz] = y_peak_vals_base_cor

        #plot baseline corrected data and the peak locations
        #p.line(x_data_numpy,y_smooth_cor,line_width=1,color=mz_colors[j],legend=str(mz_to_plot[j]))
        #p.scatter(x_peak_locs,y_peak_vals_base_cor,marker='circle',color='red',alpha=0.5,size=10)
        #p.scatter(borders_prior_sat,borders_prior_values,marker='circle',color='blue',alpha=0.5,size=5)
        #p.scatter(borders_after_sat,borders_after_values,marker='square',color='green',alpha=0.5,size=5)

        j=j+1 #index of mz

    return(ic_smooth_dict,peak_start_t_dict,peak_end_t_dict,peak_start_i_dict,peak_end_i_dict,x_data_numpy,peak_i_dict,peak_max_dict,p,peak_sat_dict,ic_dict)
