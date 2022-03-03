# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 10:06:53 2021

@author: Gruppe Willitsch - Pietro
This script, similar to main.py above, is made to study the effect of the
integration time window onto the peak height and SNR of the guided and 
decelerated peak.
Dong and Dominik where using 4e-7, 2.5e-6 as time range, which is also the default on the HyT 0.1 of Claudio. 

Enough smoothing must be done, above 1.3

"""

import numpy as np
import matplotlib.pyplot as plt
# import matplotlib
# matplotlib.use('TkAgg')
import load_dataset

plt.close('all')

# guiding path, file is 005
# path_of_dataset = r'C:\Users\Gruppe Willitsch\Desktop\Lab315\Lab315\DA_scripts\da_hyt_two\int_win_analysis\2021-02-10-005.h5'
# path_of_dataset = r'/Users/pietro/Documents/PhD/coding/Lab315/Lab315/DA_scripts/da_hyt_two/int_win_analysis/2021-02-18-005.h5'
# path_of_dataset = r'/home/pietro/Documents/PhD/coding/Lab315/Lab315/DA_scripts/da_hyt_two/int_win_analysis/2021-02-18-005.h5'
# deceleration path high res, 40 degrees, file is 017
path_of_dataset = r'/Users/pietro/Documents/PhD/coding/Lab315/Lab315/DA_scripts/da_hyt_two/int_win_analysis/2021-02-18-007.h5'
# path_of_dataset = r'C:\Users\Gruppe Willitsch\Desktop\Lab315\Lab315\DA_scripts\da_hyt_two\int_win_analysis\2021-02-19-017.h5'
# deceleration path, file is 007
# path_of_dataset = r'C:\Users\Gruppe Willitsch\Desktop\Lab315\Lab315\DA_scripts\da_hyt_two\int_win_analysis\2021-02-18-007.h5'

smoothing = 5
file_num = 5
find_peak_pos = list(range(200, 212))

my2 = load_dataset.SingleDataset(path_of_dataset, file_num,
                                 bin_value=0, smoothing_value=smoothing)
#    [0.2e-6, 5e-6], [0., 5e-6]
# my2.plot_lif_signal(), my2.plot_signal_and_background()
# my2.find_peaks(plotting=False)

# generate appropriate time window range from eg 0e-6 to 2.2e-6, in 0.1 us steps.
start = np.arange(0.e-6, .5e-6, .25e-6)
end = np.arange(5e-6, 2e-6, -2.5e-6)
int_time_win = np.zeros([len(start), len(end), 2], np.float64)
for i in range(len(start)):
    for j in range(len(end)):
        int_time_win[i, j, :] = [start[i], end[j]]

# for every time window, get a new set of data
peak_height_signal = np.zeros([len(start), len(end)], np.float64)
peak_height_difference = np.zeros([len(start), len(end)], np.float64)
SNR = np.zeros([len(start), len(end)], np.float64)

for i in range(len(start)):
    for j in range(len(end)):
        print(i, j, int_time_win[i, j, :])
        my_dataset = load_dataset.SingleDataset(path_of_dataset, file_num, smoothing_value=smoothing,
                                                signal_integration_time_win=int_time_win[i, j],
                                                background_integration_time_win=int_time_win[i, j])
        # peak finder
        [signal_peaks, signal_peaks_heights, difference_peaks, difference_peaks_heights] = my_dataset.find_peaks(
            plotting=False, verbose=True)
        #        find_peak_pos = [109, 110, 111, 108]
        if np.any(np.isin(signal_peaks, find_peak_pos)) and np.any(
                np.isin(difference_peaks, find_peak_pos)):  # guiding peak is at 109
            print(np.any(np.isin(signal_peaks, find_peak_pos)), '\t', signal_peaks)
            peak_height_signal[i, j] = signal_peaks_heights[np.isin(signal_peaks, find_peak_pos)]
            peak_height_difference[i, j] = difference_peaks_heights[np.isin(difference_peaks, find_peak_pos)]
        else:
            print('Ahiahiahi quanti guai in Paraguay\nI did NOT FOUND THE PEAK in ', find_peak_pos)
            peak_height_signal[i, j] = 0.0005
            peak_height_difference[i, j] = 0.0005

        # SNR estimation
        mean_signal, mean_background = np.sum(my_dataset.analyzed_signal), np.sum(my_dataset.analyzed_background)
        SNR[i, j] = mean_signal / mean_background

my_max = np.max([np.max(peak_height_signal), np.max(peak_height_difference)])
my_min = np.min([np.min(peak_height_signal), np.min(peak_height_difference)])
plt.figure()
plt.subplot(1, 2, 1)
plt.imshow(peak_height_signal, origin='lower', aspect='auto', interpolation='none', vmin=my_min, vmax=my_max,
           extent=[end[0] * 1e6, end[-1] * 1e6, start[0] * 1e6, start[-1] * 1e6])
plt.title('Central peak height, signal')
plt.colorbar()
plt.xlabel('End of integration window (us)')
plt.ylabel('Start of integration window (us)')

plt.subplot(1, 2, 2)
plt.imshow(peak_height_difference, origin='lower', aspect='auto', interpolation='none', vmin=my_min,
           extent=[end[0] * 1e6, end[-1] * 1e6, start[0] * 1e6, start[-1] * 1e6])
plt.title('Central peak height, signal-background')
plt.colorbar()
plt.xlabel('End of integration window (us)')
plt.ylabel('Start of integration window (u s)')

plt.figure()
plt.imshow(SNR, origin='lower', aspect='auto', interpolation='none',
           extent=[end[0] * 1e6, end[-1] * 1e6, start[0] * 1e6, start[-1] * 1e6])
plt.colorbar()
plt.title('Signal to noise ratio')
plt.xlabel('End of integration window (us)')
plt.ylabel('Start of integration window (u s)')

# print maxima of peak_height_signal, peak_height_difference, SNR
# peak height signal
i, j = np.unravel_index(peak_height_signal.argmax(), peak_height_signal.shape)
print('Peak height signal has maximum value of ', np.max(peak_height_signal), '\t for time win of ', int_time_win[i, j])
i, j = np.unravel_index(peak_height_difference.argmax(), peak_height_difference.shape)
print('Peak height difference has maximum value of ', np.max(peak_height_difference), '\t for time win of ',
      int_time_win[i, j])
i, j = np.unravel_index(SNR.argmax(), SNR.shape)
print('SNR has maximum value of \t\t', np.max(SNR), '\t for time win of ', int_time_win[i, j])
