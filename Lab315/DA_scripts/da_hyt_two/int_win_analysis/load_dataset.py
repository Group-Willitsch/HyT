#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 11:14:59 2021

@author: pietro
Notes:
    -np.mean could be replaced by np.nanmean to avoid bat things were
-CHECK FROM WERE RANGE STARTS
-direct subtraction of lif signal does not work because of difference length, 
even with same integration window. Not good, they should have the same length,
 where one vetor has 2 more points.
 -shall I renormalize the LIF signal by divind it by the number of averages and number of points? Boh. 
 TODO: normalize properly the lif signal in both classes. What changeds if u take a dataset with 40 or 10 avgs? Please check
 TODO: extract the amplitude of the laser peak and try to use it to monitor laser power fluctuations. 
"""
import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal._peak_finding import find_peaks

from get_dataset_path import get_path_of_one_dataset  # need only for __main__

# very ugly but it used in both classes as default
default_sign_int_time_win = [0e-6, 10e-6]
default_bckgrnd_int_time_win = [0e-6, 10e-6]


class SingleDataset:
    '''
    This class loads and store ONE SINGLE dataset.
    Default is INTEGRATED TIME SIGNAL of background and 

    Variables:

    Methods:


    reminder for h5py: they have two kind of objects, datasets and groups.
    Groups work like dictionaries, and datasets work like NumPy arrays
    '''

    # class static variables, bounded to the class and not to the instance of the class

    def __init__(self, path_of_dataset, file_number,
                 signal_integration_time_win=default_sign_int_time_win,  # default signal integration time
                 background_integration_time_win=default_bckgrnd_int_time_win,  # default background integration time
                 bin_value=0,  # 0 = no bin_value, 2, 4 or 8 means I generated a copy of the dataset, binned.
                 smoothing_value=0):  # 0 = no smoothing, > 0 smooths the dataset (not the raw data)
        # extract the dataset (kind of degenerate case)
        [self.files_in_folder_full_path, _, _, self.file_numbers] = get_path_of_one_dataset(path_of_dataset,
                                                                                            file_number, file_number,
                                                                                            debug_mode=False,
                                                                                            list_files_nums_to_skip=[])
        if len(self.files_in_folder_full_path) > 1:  # to be removed
            print('***Ahhhh something wrong in single dataset init!***')

        # works like Python dictionary
        # File object is itself a group (root group).        
        self.file = h5py.File(self.files_in_folder_full_path[0], 'r')
        # attrbutes of this h5 file are saved in self.attributes
        self.attributes = self.file.attrs
        # with lst(self.attributes): ['averages', 'bgchannel', 'osci_tbins',
        # 'osci_tstart', 'osci_tstop', 'relchannel', 'scan0', 'scan0_time',
        # 'scan1', 'scan1_time', 'scanchannel', 'scanrel0', 'scanrel1']
        # to extract each of them, use:
        #            self.bgchannel = self.file['inputs/bgchannel'][()]
        self.scan_length = len(list(self.file['signal'].keys()))
        self.signal_integration_time_win = signal_integration_time_win
        self.background_integration_time_win = background_integration_time_win
        self.inputs = self.file['inputs']

        # each dataset has this structure:
        #    'background'      (no attrs)
        #         '0.0023'     (shape 25002, attributes x1, x2, xb1, xb2, y, y_oc)
        # they are 'background', 'background_pd', 'signal', 'signal_pd'

        # the timebase is the same for all the 4 data super-sets: signal, signal_pd, background, background_pd
        self.timebase_np = np.array(list(self.file['signal'].keys()),
                                    dtype=np.float64)  # stored in a proper numpy array

        # The "main" timebase is the timebase of the laser delay scan. The other
        # one is the timebase of the signal acquired by the scope.
        # which is necessary to, for example, gate the dataset and plot LIF traces.
        self.lif_timebase = np.linspace(self.attributes['osci_tstart'],
                                        self.attributes['osci_tstop'], self.attributes['osci_tbins'])

        # check for nan values in the dataset
        # normally commented has it takes quite some time (half of total running time)
        # self.nan_values = self.check_for_nan_values()

        # load the signal nd the background for the given integration values
        [self.signal_np, self.background_np] = self.integrate_in_given_time_window(
            self.signal_integration_time_win, self.background_integration_time_win)

        # load the scope traces of the lif signal (to be plot against self.lif_timebase)
        [self.lif_signal_np, self.lif_background_np, self.lif_difference_np] = self.get_lif_signal()

        if bin_value > 0 and smoothing_value > 0:  # bad but not too bad: simply binning operation gets overwritten
            print('Moderate warning: you are both binning and smoothing.\nCurrently only smoothing will succeed.')
        # data bin_value (if bin_value = 0, analyzed_data
        self.bin_value = bin_value
        self.analyzed_signal, self.analyzed_background = self.get_binned_dataset(self.bin_value)

        # smooth the dataset
        if smoothing_value > 0:  # otherwise I overwrite the binned value, also then smoothing_value=0
            self.smoothing_value = smoothing_value
            self.analyzed_signal, self.analyzed_background = self.get_smoothed_dataset(self.smoothing_value)

        # compute the SNR for this dataset. For his definition, check under code und sim -> int. win. analysis
        self.SNR = np.mean(self.analyzed_signal)/np.mean(self.analyzed_background)-1

        # some plots setup
        self.xlabel = 'laser delay (us)'
        self.ylabel = 'PMT signal (a.u.)'
        print('Loaded dataset\t ...', self.files_in_folder_full_path[0][-36:], '\r')
        print('SNR for this dataset is:\t\t', round(self.SNR, 4), '\r')

    def check_for_nan_values(self):
        bool_nan_found = False
        nan_i_ve_found = 0
        # check nan directly in the raw data, not in the integrated data.
        # this may be a bit long
        for laser_delay_str in self.file['background']:
            # this two lines get the oscilloscope signal
            if np.isnan(np.sum(np.array(self.file['signal/' + laser_delay_str][()], np.float64))) or np.isnan(
                    np.sum(np.array(self.file['signal/' + laser_delay_str][()], np.float64))):  # ugly, srry
                bool_nan_found = True
        if bool_nan_found:  # then I have to count how many of them, spit an error and name of dataset
            nan_i_ve_found = 1  # arbitrary, still to be implemented
            print('*** I have found ', nan_i_ve_found, ' nan values in dataset ', self.files_in_folder_full_path,
                  ' ***\n***  AHHHHHHH Panic  ***\n')
        return nan_i_ve_found

    def integrate_in_given_time_window(self, signal_integration_time_win, background_integration_time_win):
        '''
        Parameters
        ----------
        signal_integration_time_win : list, [0, 10e-6] in seconds
        background_integration_time_win : ist, [0, 10e-6] in seconds

        Returns
        -------
        signal and background vectors, as np_array with a length of self.scan_length
        May fail if nan are encountered.
        The mean is computed forcing float64 operation, float32 is default and may be inaccurate. Worth a check.
        '''
        # generate two boolean gate arrays, one for the signal one for the background
        gate_for_signal = np.logical_and(self.lif_timebase >= signal_integration_time_win[0],
                                         self.lif_timebase <= signal_integration_time_win[1])
        gate_for_background = np.logical_and(self.lif_timebase >= background_integration_time_win[0],
                                             self.lif_timebase <= background_integration_time_win[1])

        signal_np = np.empty(self.scan_length, np.float64)  # proper init of empty signal and background
        background_np = np.empty(self.scan_length, np.float64)
        # cycle in signal database, but the counter is good also for the other 3 datasets
        i = 0
        for laser_delay_str in self.file['background']:
            # this two lines get the oscilloscope signal and do the gating too
            signal_trace = np.array(self.file['signal/' + laser_delay_str][()], np.float64)[gate_for_signal]
            background_trace = np.array(self.file['background/' + laser_delay_str][()], np.float64)[gate_for_background]
            #            signal_np[i] = np.mean( signal_trace - background_trace, dtype=np.float64) #does not work, different length even with same integration window
            signal_np[i] = np.mean(signal_trace, dtype=np.float64)  # forcing float64 for more accuracy, see docs
            background_np[i] = np.mean(background_trace, dtype=np.float64)
            i += 1
        return signal_np, background_np
        # note: I don not return self.signal_np but signal_np, s.t. the function
        # can be used without overwriting the variables of the dataset.

    def get_lif_signal(self):
        """
        Sum up the whole dataset, regardless of the laser delay time stamp, to
        get a LIF signal.
        NOTE: SHOULD I RENORMALIZE THE THREE VECTORS??

        Returns
        -------
        lif_signal_np : numpy array, of length len(self.lif_timebase), float 64
            total integrated signal.
        lif_background_np : idem
            total integrated background.
        lif_difference_np : idem
            direct term-by-term difference of the two vectors above.

        """
        # proper init, to be sure everything is float64.
        lif_signal_np = np.zeros(len(self.lif_timebase), np.float64)
        lif_background_np = np.zeros(len(self.lif_timebase), np.float64)
        for laser_delay_str in self.file['signal']:  # I cycle on signal but works for all
            lif_signal_np += np.array(self.file['signal/' + laser_delay_str][()], np.float64)
            lif_background_np += np.array(self.file['background/' + laser_delay_str][()], np.float64)
        lif_difference_np = lif_signal_np - lif_background_np  # direct background subtraction
        return lif_signal_np, lif_background_np, lif_difference_np

    def plot_signal_and_background(self):
        """
        By default plots the analyzed dataset.
        If not binning/smoothing/further analysis is being made, it will plot automatically the raw data
        """
        plt.figure()
        legend = ['signal', 'background', 'difference']
        plt.plot(self.timebase_np * 1e6, self.analyzed_signal,
                 self.timebase_np * 1e6, self.analyzed_background,
                 self.timebase_np * 1e6, self.analyzed_signal - self.analyzed_background)
        plt.fill_between(self.timebase_np * 1e6, 0, self.analyzed_signal - self.analyzed_background, alpha=0.5)
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.legend(legend)

    def plot_lif_signal(self, print_mean_values=False):
        plt.figure()
        legend = ['signal', 'background', 'difference']
        plt.plot(self.lif_timebase * 1e6, self.lif_signal_np, self.lif_timebase * 1e6,
                 self.lif_background_np, self.lif_timebase * 1e6,
                 self.lif_difference_np, lw=.8, ms=0.2)
        plt.xlabel(self.xlabel)
        plt.ylabel('LIF signal')
        plt.legend(legend)
        plt.ylim([-0.2, np.max(self.lif_difference_np) * 3])
        if print_mean_values:
            print('Mean of total integrated  signal: \t', np.round(np.mean(self.lif_signal_np), 4), ' +- ',
                  np.round(np.std(self.lif_signal_np), decimals=2), '\r')
            print('Mean of total integrated  background: \t', np.round(np.mean(self.lif_background_np), decimals=4),
                  ' +- ', np.round(np.std(self.lif_background_np), decimals=2), '\r')
            print('Mean of total integrated  difference: \t', np.round(np.mean(self.lif_difference_np), decimals=4),
                  ' +- ', np.round(np.std(self.lif_timebase), decimals=4), '\r')
            print('Signal to noise ratio: \t\t\t',
                  np.round(np.mean(self.lif_signal_np) / np.mean(self.lif_background_np), decimals=4), '\r\n')
            print("***This datas for sure are not normalized!!!!!1***")

    def get_binned_dataset(self, bin_value):
        if bin_value not in 2 ** np.arange(8) and bin_value != 0:
            print("Your bin_value value is not in [0, 2, 4, 8]. It is", bin_value)
        analyzed_signal, analyzed_background = self.signal_np, self.background_np  # default is non-processed dataset
        if bin_value in 2 ** np.arange(8):  # I do the bin_value; zero is excluded
            excess = len(self.signal_np) - len(self.signal_np) // bin_value * bin_value
            # get a shorter vector but with mean values
            analyzed_signal = np.mean(
                np.reshape(self.signal_np[:-excess], (len(self.signal_np) // bin_value, bin_value)), 1)
            analyzed_background = np.mean(
                np.reshape(self.background_np[:-excess], (len(self.background_np) // bin_value, bin_value)), 1)
            # re-shape a vector with the same size at the self.signal_np
            analyzed_signal = np.append(np.repeat(analyzed_signal, bin_value), self.signal_np[-excess:])
            analyzed_background = np.append(np.repeat(analyzed_background, bin_value), self.background_np[-excess:])
            print("Dataset binned by", bin_value, " bins")
        return analyzed_signal, analyzed_background

    def get_smoothed_dataset(self, smoothing_value):
        analyzed_signal, analyzed_background = self.signal_np, self.background_np  # default is non-processed dataset
        if smoothing_value > 0:
            print('Smoothing the dataset with smoothing_val of ', smoothing_value)
            analyzed_signal = gaussian_filter1d(analyzed_signal, smoothing_value)
            analyzed_background = gaussian_filter1d(analyzed_background, smoothing_value)
        return analyzed_signal, analyzed_background

    def find_peaks(self, peak_window_indexes=[150, 250], my_prominence=1, my_width=2, plotting=True, verbose=True):
        [signal_peaks, _] = find_peaks(self.analyzed_signal[peak_window_indexes[0]:peak_window_indexes[1]], prominence=my_prominence, width=my_width)
        signal_peaks = signal_peaks + peak_window_indexes[0] # sum what was subtracted
        signal_peaks_heights = self.analyzed_signal[signal_peaks]
        [difference_peaks, _] = find_peaks((self.analyzed_signal-self.analyzed_background)[peak_window_indexes[0]:peak_window_indexes[1]], prominence=my_prominence, width=my_width)
        difference_peaks = difference_peaks + peak_window_indexes[0]
        difference_peaks_heights = (self.analyzed_signal-self.analyzed_background)[difference_peaks]
        if plotting == True:
            plt.figure()
            plt.plot(self.analyzed_signal, alpha=0.5)
            plt.scatter(signal_peaks, self.analyzed_signal[signal_peaks], lw=5)
            plt.plot(self.analyzed_signal-self.analyzed_background, alpha=0.5)
            plt.scatter(difference_peaks, (self.analyzed_signal-self.analyzed_background)[difference_peaks], lw=3)
            plt.xlabel('array index (can be converted to laser delay)')
            plt.ylabel('signal intensity')
            plt.legend(['signal', 'signal-background'])
        if verbose:
            print('Peaks in signal:\t', signal_peaks, '\nwith values:\t\t', signal_peaks_heights)
            print('Peaks in difference:\t', difference_peaks, '\nwith values:\t\t', difference_peaks_heights)

        return signal_peaks, signal_peaks_heights, difference_peaks, difference_peaks_heights

class MultiplesCoherentDataset(SingleDataset):

    def __init__(self, path_of_dataset, first_file_number,
                 last_file_number, list_files_nums_to_skip=[],
                 signal_integration_time_win=default_sign_int_time_win,
                 background_integration_time_win=default_bckgrnd_int_time_win,
                 bin_value=0,
                 smoothing_value=0):  # bin value will be passed to SingleDataset constructor
        print('__init__ of multiples coherent dataset')

        # retrieve the proper file structure
        [self.files_in_folder_full_path, _, _, self.file_numbers] = get_path_of_one_dataset(path_of_dataset,
                                                                                            first_file_number,
                                                                                            last_file_number,
                                                                                            debug_mode=False,
                                                                                            list_files_nums_to_skip=list_files_nums_to_skip)
        self.number_of_datasets = len(self.files_in_folder_full_path)  # how many files

        # this acquires the whole dataset and stores it in only one variable
        self.complete_dataset = [SingleDataset(path_of_dataset, self.file_numbers[i], signal_integration_time_win,
                                               background_integration_time_win, bin_value, smoothing_value) for i in
                                 range(self.number_of_datasets)]
        # kaboom

        self.is_timebase_coherent = self.check_timebase_coherence()  # just store the bool

        # data that pertains to the whole dataset
        [self.signal_matrix_y_val, self.background_matrix_y_val] = self.get_matrices()
        # lif data
        self.lif_timebase = self.complete_dataset[0].lif_timebase
        [self.lif_signal_np, self.lif_background_np, self.lif_difference_np] = self.get_lif_signal()

        # plotting
        self.xlabel = 'laser delay (ms)'
        self.ylabel_intensity = 'signal on PMT (a.u.)'
        self.ylabel_2dscan = 'filenames'

    def check_timebase_coherence(self):
        '''
        compare every self.timebase_np of every dataset among each other.
        If coherent ==> return True, otherwise, return False
        TODO: Just trow an expection, non fatal such that one can still look at what happened. 
        '''
        for i in range(self.number_of_datasets - 1):
            if not np.array_equal(self.complete_dataset[i].timebase_np, self.complete_dataset[i + 1].timebase_np):
                print('I have found a mismatch between timebase of this dataset.')
                print('Problem between ', self.files_in_folder_full_path[i], 'and ',
                      self.files_in_folder_full_path[i + 1], sep='\n')
                print('Ahihaihai quanti guai in Paraguay\n')
                return False
        return True

    def get_matrices(self):
        self.signal_matrix_y_val = np.zeros([self.number_of_datasets, len(self.complete_dataset[0].signal_np)],
                                            np.float64)
        self.background_matrix_y_val = np.zeros([self.number_of_datasets, len(self.complete_dataset[0].background_np)],
                                                np.float64)
        for i in range(self.number_of_datasets):
            self.signal_matrix_y_val[i, :] = self.complete_dataset[i].analyzed_signal
            self.background_matrix_y_val[i, :] = self.complete_dataset[i].analyzed_background

        # integrate over x axis (laser delay) and over y axis (filename)
        self.signal_integrated_over_filenames = np.mean(self.signal_matrix_y_val, axis=0)
        self.background_integrated_over_filenames = np.mean(self.background_matrix_y_val, axis=0)
        self.signal_integrated_over_laser_delay = np.mean(self.signal_matrix_y_val, axis=1)
        self.signal_integrated_over_laser_delay_std = np.std(self.signal_matrix_y_val, axis=1)
        self.background_integrated_over_laser_delay = np.mean(self.background_matrix_y_val, axis=1)
        self.background_integrated_over_laser_delay_std = np.std(self.background_matrix_y_val, axis=1)
        return self.signal_matrix_y_val, self.background_matrix_y_val

    def get_lif_signal(self):
        '''
        I inherit the method from single classes and just sum up everthing toghether.
        Returns
        -------
        lif_signal_np : numpy array, of length len(self.lif_timebase), float 64
            total integrated signal.
        lif_background_np : idem
            total integrated background.
        lif_difference_np : idem
            direct term-by-term difference of the two vectors above.

        '''
        # proper init, to be sure everything is float64.
        lif_signal_np = np.zeros(len(self.lif_timebase), np.float64)
        lif_background_np = np.zeros(len(self.lif_timebase), np.float64)
        lif_difference_np = np.zeros(len(self.lif_timebase), np.float64)
        for i in range(self.number_of_datasets):
            lif_signal_np += self.complete_dataset[i].lif_signal_np
            lif_background_np += self.complete_dataset[i].lif_background_np
        lif_difference_np = lif_signal_np - lif_background_np  # direct background substraction
        return lif_signal_np, lif_background_np, lif_difference_np

    def plot_signal_and_background(self):
        plt.figure()

        # plot signal integrated over the filenames
        plt.subplot(2, 2, 3)
        plt.plot(self.complete_dataset[0].timebase_np * 1e3, self.signal_integrated_over_filenames, lw=3, ms=8)
        plt.xlabel(self.xlabel), plt.ylabel(self.ylabel_intensity), plt.title('Signal integrated over all the files')

        # plot the signal integrated over the laser delays
        plt.subplot(2, 2, 4)
        legend = ['signal', 'background', 'differnce']
        plt.plot(self.file_numbers, self.signal_integrated_over_laser_delay, lw=3)
        plt.plot(self.file_numbers, self.background_integrated_over_laser_delay)
        plt.plot(self.file_numbers,
                 self.signal_integrated_over_laser_delay - self.background_integrated_over_laser_delay)
        plt.xlabel(self.ylabel_2dscan), plt.ylabel(self.ylabel_intensity), plt.title(
            'Signal integrated over laser delay')
        plt.legend(legend)
        # 2D plot
        plt.subplot(2, 1, 1)
        plt.imshow(self.signal_matrix_y_val, origin='lower', aspect='auto',
                   interpolation='none',
                   extent=[self.complete_dataset[0].timebase_np[0] * 1e3,
                           self.complete_dataset[0].timebase_np[-1] * 1e3,
                           self.file_numbers[0], self.file_numbers[-1]])
        plt.grid(None), plt.colorbar()
        plt.xlabel(self.xlabel), plt.ylabel(self.ylabel_2dscan)

        plt.figure()
        for i in range(self.number_of_datasets):
            #            plt.fill_between(self.complete_dataset[0].timebase_np * 1e3, 0, self.signal_matrix_y_val[i, :],
            #                             alpha=0.25)
            plt.plot(self.complete_dataset[0].timebase_np * 1e3,
                     self.signal_matrix_y_val[i, :] - self.background_matrix_y_val[i, :] + i * 0.002, lw=4,
                     alpha=-0.05 * i + 0.7, ms=20)
        #        plt.ylim([np.my_min(self.signal_matrix_y_val) - 1E-4, np.my_max(self.signal_matrix_y_val) + 1e-4])
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel_2dscan)
        legend = self.file_numbers
        plt.legend(legend)

        # plot only the background. Does it correlates with laser fluctuations??
        plt.figure()
        plt.title('Mean and std of laser background of the whole filename')
        plt.errorbar(self.file_numbers, self.signal_integrated_over_laser_delay,
                     yerr=self.signal_integrated_over_laser_delay_std,
                     alpha=0.5, ms=5, lw=4, capsize=10)
        plt.errorbar(self.file_numbers, self.background_integrated_over_laser_delay,
                     yerr=self.background_integrated_over_laser_delay_std, alpha=0.8, ms=5, lw=4, capsize=40)
        plt.errorbar(self.file_numbers,
                     self.signal_integrated_over_laser_delay - self.background_integrated_over_laser_delay,
                     alpha=0.5, ms=5, lw=4, capsize=10)
        plt.legend(['signal', 'background', 'difference'])
        plt.xlabel(self.ylabel_2dscan), plt.ylabel(self.ylabel_intensity)

        # plot only the 2d image
        plt.figure()
        plt.imshow(self.signal_matrix_y_val, origin='lower', aspect='auto',
                   interpolation='none',
                   extent=[self.complete_dataset[0].timebase_np[0] * 1e3,
                           self.complete_dataset[0].timebase_np[-1] * 1e3,
                           self.file_numbers[0], self.file_numbers[-1]])
        plt.grid(None), plt.colorbar()
        plt.xlabel(self.xlabel), plt.ylabel(self.ylabel_2dscan)

    def plot_lif_signal(self, print_mean_values=False):
        plt.figure()
        legend = ['signal', 'background', 'difference']
        plt.plot(self.lif_timebase * 1e6, self.lif_signal_np, self.lif_timebase * 1e6, self.lif_background_np,
                 self.lif_timebase * 1e6, self.lif_difference_np, lw=.8, ms=0.2)
        plt.xlabel(self.complete_dataset[0].xlabel)
        plt.ylabel('LIF signal'), plt.title('LIF signal from the whole dataset')
        plt.legend(legend)
        plt.ylim([-0.2, np.max(self.lif_signal_np) * .1])

        # print the mean values, std and ratio of total integrated signal and background
        if print_mean_values:
            print('Mean of total integrated  signal: \t', np.round(np.mean(self.lif_signal_np), 4), ' +- ',
                  np.round(np.std(self.lif_signal_np), decimals=2), '\r')
            print('Mean of total integrated  background: \t', np.round(np.mean(self.lif_background_np), decimals=4),
                  ' +- ', np.round(np.std(self.lif_background_np), decimals=6), '\r')
            print('Mean of total integrated  difference: \t', np.round(np.mean(self.lif_difference_np), decimals=4),
                  ' +- ', np.round(np.std(self.lif_timebase), decimals=6), '\r')
            print('Signal to noise ratio: \t\t\t',
                  np.round(np.mean(self.lif_signal_np) / np.mean(self.lif_background_np), decimals=4), '\r\n')
            print("***This datas for sure are not normalized!!!!!1***")


if __name__ == "__main__":
    first_file_number = 40
    last_file_number = 52
    files_nums_to_skip = [46]
    path_of_dataset = r'C:\Users\Gruppe Willitsch\Desktop\Lab315\Lab315\DA_scripts\da_hyt_two\test_dataset/2021-01-07-040.h5'
    #    path_of_dataset = r'/Users/pietro/Documents/PhD/coding/Lab315/Lab315/DA_scripts/da_hyt_two/test_dataset/2021-01-07-041.h5'

    my_first_dataset = SingleDataset(path_of_dataset, 41, bin_value=2, smoothing_value=0)
    my_first_dataset.plot_signal_and_background()
    my_first_dataset.plot_lif_signal()

    iterator = False
    if iterator:  # cycle onto the integration windows
        signal_time_wnd = [[0, 5e-6], [0.2e-6, 5e-6], [0.4e-6, 5e-6], [1e-6, 5]]
        plt.figure()
        for i in range(len(signal_time_wnd[:][:])):
            print(i)
            print(signal_time_wnd[1][:])
            x, y = my_first_dataset.integrate_in_given_time_window(signal_time_wnd[:][i], [0, 10e-6])
            plt.plot(my_first_dataset.timebase_np * 1000, x, linewidth=4)
            legend = ['0', '1', '2', '3']
            plt.legend(legend)

    # MULTIPLE COHERENT DATASET:
    mine_multiples_coherent_dataset = MultiplesCoherentDataset(
        path_of_dataset, 40, 50, [46],
        bin_value=0, smoothing_value=0.5)
    mine_multiples_coherent_dataset.check_timebase_coherence()
    mine_multiples_coherent_dataset.plot_signal_and_background()
    #    mine_multiples_coherent_dataset.plot_lif_signal()
    plt.show()
