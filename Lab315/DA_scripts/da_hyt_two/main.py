#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 11:08:44 2021

@author: pietro
"""
import load_dataset
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal._peak_finding import find_peaks
from scipy.ndimage import gaussian_filter1d
global_bool_for_photn_counting = True
plt.close('all')
new_sign_int_time_win = [0.3e-6, 5e-6]
new_bckgrnd_int_time_win = [0.3e-6, 5e-6]
# today
path_of_dataset = r'C:/Users/Gruppe Willitsch/Desktop/data/2022-01-11/2022-01-11-001.h5'
tt = load_dataset.SingleDataset(path_of_dataset, 18, photon_counting=global_bool_for_photn_counting,
                                # signal_integration_time_win=new_sign_int_time_win,
                                # background_integration_time_win=new_bckgrnd_int_time_win,
                                smoothing_value=.4)
tt.plot_signal_and_background()

tt2 = load_dataset.SingleDataset(path_of_dataset, 3, photon_counting=global_bool_for_photn_counting,
                                 # signal_integration_time_win=new_sign_int_time_win,
                                 # background_integration_time_win=new_bckgrnd_int_time_win,
                                 smoothing_value=.4)
tt2.plot_signal_and_background()
print('Mean background first file', np.mean(tt.analyzed_background))
print('Mean background second file', np.mean(tt2.analyzed_background))
# tt3 = load_dataset.SingleDataset(path_of_dataset, 19, photon_counting=global_bool_for_photn_counting, smoothing_value=.4)
# tt3.plot_signal_and_background(wavelength_scan=True)





# best signals evah
path_of_dataset_ref = r'C:/Users/Gruppe Willitsch/Desktop/data/2022-01-06/2022-01-06-002.h5'
# ref_90 = load_dataset.SingleDataset(path_of_dataset_ref, 56, photon_counting=global_bool_for_photn_counting, smoothing_value=.4)
ref_49 = load_dataset.SingleDataset(path_of_dataset_ref, 79, photon_counting=global_bool_for_photn_counting, smoothing_value=1.4)
# ref_30 = load_dataset.SingleDataset(path_of_dataset_ref, 82, photon_counting=global_bool_for_photn_counting, smoothing_value=1.4)




plt.figure()
plt.plot(tt.timebase_np*1e6, (tt.analyzed_signal-tt.analyzed_background), label = 'file n. ' + str(tt.file_numbers[0]) + ' of 11/01/2022 \n -265 us laser power')
plt.plot(tt2.timebase_np*1e6, (tt2.analyzed_signal-tt2.analyzed_background), label = 'file n. ' + str(tt2.file_numbers[0]) + ' of 11/01/2022 \n -20 us laser power')
# plt.plot(tt3.timebase_np*1e6, (tt3.analyzed_signal-tt3.analyzed_background), label = 'file n. ' + str(tt3.file_numbers[0]))
# plt.plot(ref_90.timebase_np*1e6, (ref_90.analyzed_signal-ref_90.analyzed_background), label = 'file n. ' + str(ref_90.file_numbers[0]) + ' of 28/06/2021, REF 90 m/s ')
plt.plot(ref_49.timebase_np*1e6, (ref_49.analyzed_signal-ref_49.analyzed_background), label = 'file n. ' + str(ref_49.file_numbers[0]) + ' of 28/06/2021, REF 49 m/s ')
# plt.plot(ref_30.timebase_np*1e6, (ref_30.analyzed_signal-ref_30.analyzed_background), label = 'file n. ' + str(ref_30.file_numbers[0]) + ' of 28/06/2021, REF 30 m/s ')

plt.xlabel(tt.xlabel)
plt.ylabel(tt.ylabel)
plt.title('30 m/s - comparison with old references')
# plt.title('92 m/s, 12.5 kV focusing mode with trap inside')
plt.legend()



# plt.close('all')
# scan = load_dataset.MultiplesCoherentDataset(path_of_dataset, 141, 159, [152], photon_counting=True)
# scan.plot_signal_and_background()
# scan.plot_repetitions()
#
# scan2 = load_dataset.MultiplesCoherentDataset(path_of_dataset, 141, 159, [152], photon_counting=True,
#                                              signal_integration_time_win=new_sign_int_time_win,
#                                              background_integration_time_win=new_bckgrnd_int_time_win)
# scan2.plot_signal_and_background()
# scan2.plot_repetitions()
#
# plt.figure()  # first standard plot
# # analyzed data
# plt.plot(scan.timebase_np * 1e6, scan.signal_avg_over_filenames - scan.background_avg_over_filenames,
#          lw=3, ms=0.3, label = 'bigger int window', color='red')
# # fill an area of +- yerr
# plt.fill_between(scan.timebase_np * 1e6, scan.signal_avg_over_filenames-scan.background_avg_over_filenames - scan.yerr,
#                  scan.signal_avg_over_filenames-scan.background_avg_over_filenames + scan.yerr, alpha=.4,
#                  color='red', lw=1)
# # plt.fill_between(scan.timebase_np * 1e6, 0., scan.signal_avg_over_filenames-scan.background_avg_over_filenames, alpha=.3, lw=0)
#
# plt.plot(scan2.timebase_np * 1e6, scan2.signal_avg_over_filenames - scan2.background_avg_over_filenames,
#          lw=3, ms=0.3, color='green', label ='smaller int window')
# # fill an area of +- yerr
# plt.fill_between(scan2.timebase_np * 1e6, scan2.signal_avg_over_filenames-scan2.background_avg_over_filenames - scan2.yerr,
#                  scan2.signal_avg_over_filenames-scan2.background_avg_over_filenames + scan2.yerr, alpha=.4,
#                  color='green', lw=1)
# plt.xlabel(scan2.xlabel), plt.ylabel(scan2.ylabel), plt.title(
#     'Signal averaged over all the files\n+- possion error ( = sqrt(rate) / (number of experiment) )')
#
# plt.legend()
#
#





path_of_dataset = r'C:/Users/Gruppe Willitsch/Desktop/data/2022-01-10/2022-01-10-001.h5'
five_hertz = load_dataset.MultiplesCoherentDataset(path_of_dataset, 90, 94, photon_counting=True)
one_hertz = load_dataset.MultiplesCoherentDataset(path_of_dataset, 98, 103, photon_counting=True)
half_hertz = load_dataset.MultiplesCoherentDataset(path_of_dataset, 104, 109, photon_counting=True)
point_two_hertz = load_dataset.MultiplesCoherentDataset(path_of_dataset, 118, 125, photon_counting=True)
point_one_hertz = load_dataset.MultiplesCoherentDataset(path_of_dataset, 126, 138, [127, 128],  photon_counting=True)
point_o_five_hertz = load_dataset.MultiplesCoherentDataset(path_of_dataset, 141, 159, [152], photon_counting=True)

tt.timebase_np = np.concatenate((five_hertz.timebase_np, one_hertz.timebase_np, half_hertz.timebase_np,
                                 point_two_hertz.timebase_np, point_one_hertz.timebase_np,
                                 point_o_five_hertz.timebase_np), axis=0)
tt.analyzed_signal = np.concatenate((five_hertz.signal_avg_over_filenames, one_hertz.signal_avg_over_filenames,
                                     half_hertz.signal_avg_over_filenames, point_two_hertz.signal_avg_over_filenames,
                                     point_one_hertz.signal_avg_over_filenames,
                                     point_o_five_hertz.signal_avg_over_filenames), axis=0)
tt.analyzed_background = np.concatenate((five_hertz.background_avg_over_filenames,
                                         one_hertz.background_avg_over_filenames,
                                         half_hertz.background_avg_over_filenames,
                                         point_two_hertz.background_avg_over_filenames,
                                         point_one_hertz.background_avg_over_filenames,
                                         point_o_five_hertz.background_avg_over_filenames), axis=0)

tt.signal_np = tt.analyzed_signal
tt.background_np = tt.analyzed_background
tt.yerr =  np.concatenate((five_hertz.yerr, one_hertz.yerr, half_hertz.yerr, point_two_hertz.yerr, point_one_hertz.yerr,
                           point_o_five_hertz.yerr), axis=0)

tt.fit_OH_lifetime(weights=True)
