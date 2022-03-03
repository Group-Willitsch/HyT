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
 TODO: check how normalization works when integrating the signal without single photons
 TODO: re-check integratin time window with a newer dataset at low velocity.
 TODO: put smaller (1 Mb) dataset in the test dirname_pressure_dataset and a small description in txt
 TODO: redefine SNR, print the mean number of photons for s/b
 TODO: add second photodiode channel for SingleCoherentDataset
 TODO: better integrate the single photon counting into the MultipleCoherentDataset
 TODO: add single photons peak height distribution for MultipleCoherentDataset
"""
import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal._peak_finding import find_peaks

from get_dataset_path import get_path_of_one_dataset  # need only for __main__
import platform

if platform.system() == 'Darwin':
    import matplotlib

    # matplotlib.use('TkAgg')  # fix a PyCharm bug on Mac

# very ugly but it is used in both classes as default
default_sign_int_time_win = [0.195e-6, 5e-6]
default_bckgrnd_int_time_win = [0.195e-6, 5e-6]
default_sign_int_time_win = [0.3e-6, 5e-6] # use this since Jan 2022
default_bckgrnd_int_time_win = [0.3e-6, 5e-6]

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
                 smoothing_value=0,  # 0 = no smoothing, > 0 smooths the dataset (not the raw data)
                 photon_counting=False,
                 second_da_channel=False):  # photodiode channel
        # extract the dataset (kind of degenerate case)
        [self.files_in_folder_full_path, _, _, self.file_numbers] = get_path_of_one_dataset(path_of_dataset,
                                                                                            file_number, file_number,
                                                                                            debug_mode=False,
                                                                                            list_files_nums_to_skip=[])
        if len(self.files_in_folder_full_path) > 1:  # to be removed
            print('***Ahhhh something wrong in single dataset init!***')

        # save input DA values into class variables
        self.bin_value = bin_value
        self.smoothing_value = smoothing_value
        self.photon_counting = photon_counting

        # works like Python dictionary
        # File object is itself a group (root group).        
        self.file = h5py.File(self.files_in_folder_full_path[0], 'r')

        # attributes of this h5 file are saved in self.attributes
        self.attributes = self.file.attrs
        # with list(self.attributes): ['averages', 'bgchannel', 'osci_tbins',
        # 'osci_tstart', 'osci_tstop', 'relchannel', 'scan0', 'scan0_time',
        # 'scan1', 'scan1_time', 'scanchannel', 'scanrel0', 'scanrel1']
        # to extract each of them, use:
        #            self.bgchannel = self.attributes['averages'][()]
        self.trigger_freqs = self.file['inputs/freqs'][()]  # trigger frequencies
        self.averages = self.attributes['averages'][()]  # extract the number of averages

        self.scan_length = len(list(self.file['signal'].keys()))
        self.signal_integration_time_win = signal_integration_time_win
        self.background_integration_time_win = background_integration_time_win
        self.inputs = self.file['inputs']

        # each dataset has this structure:
        #    'background'      (no attrs)
        #         '0.0023'     (shape 25002, attributes x1, x2, xb1, xb2, y, y_oc)
        # they are 'background', 'background_pd', 'signal', 'signal_pd'

        # the timebase is the same for all the 4 data super-sets: signal, signal_pd, background, background_pd
        #        print(max(self.file['signal'].keys()))
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
        # added later: refined photon counting
        if self.photon_counting:
            [self.signal_np, self.background_np] = self.count_single_photons()
        else:
            [self.signal_np, self.background_np] = self.integrate_in_given_time_window(
                self.signal_integration_time_win, self.background_integration_time_win)

        # load the scope traces of the lif signal (to be plot against self.lif_timebase)
        [self.lif_signal_np, self.lif_background_np, self.lif_difference_np] = self.get_lif_signal()

        if self.bin_value > 0 and self.smoothing_value > 0:  # bad but not too bad: simply binning operation gets overwritten
            print('Moderate warning: you are both binning and smoothing.\nCurrently only smoothing will succeed.')
        # data bin_value (if bin_value = 0, analyzed_data
        self.analyzed_signal, self.analyzed_background = self.get_binned_dataset(self.bin_value)

        # smooth the dataset
        if self.smoothing_value > 0:  # otherwise I overwrite the binned value, also then smoothing_value=0
            self.analyzed_signal, self.analyzed_background = self.get_smoothed_dataset(self.smoothing_value)



        # define errors
        self.yerr = np.zeros(len(self.analyzed_signal))  # currently not implemented; implemented in Multiples...

        # acquire the photodiode trace and make one avg for every laser delay point
        if second_da_channel:
            [self.photodiode_signal_np, self.photodiode_background_np] = self.get_photodiode_signal()

        # compute the SNR for this dataset. For his definition, check under code und sim -> int. win. analysis
        self.SNR = np.mean(self.analyzed_signal) / np.mean(self.analyzed_background) - 1

        # some plots setup
        self.xlabel = 'laser delay (us)'
        self.ylabel = 'integrated PMT signal (a.u.)'
        if self.photon_counting:
            self.ylabel = '# photons (/1 laser shot)'
        print('Loaded dataset\t ...', self.files_in_folder_full_path[0][-36:], '\r')
        print('This dataset is taken with ', self.averages, ' averages')

    #        print('SNR for this dataset is:\t\t', round(self.SNR, 4), '\r') # obsolete, wrong SNR definition

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
        # initialize two more timebases, which are gated properly. Used mainly to plot.
        self.gated_signal_timebase = self.lif_timebase[gate_for_signal]
        self.gated_background_timebase = self.lif_timebase[gate_for_background]

        signal_np = np.zeros(self.scan_length, np.float64)  # proper init of empty signal and background
        background_np = np.zeros(self.scan_length, np.float64)
        # cycle in signal database, but the counter is good also for the other 3 datasets
        i = 0
        self.jitter_arrival_times = np.empty(self.scan_length, np.float64)
        for laser_delay_str in self.file['background']:
            # this two lines get the oscilloscope signal and do the gating too
            self.jitter_arrival_times[i] = np.abs(
                np.array(self.file['signal/' + laser_delay_str][()], np.float64)[0:720] - 0.1).argmin()
            #           self.ts = np.array(self.file['signal/' + laser_delay_str][()], np.float64)
            #           self.tb = np.array(self.file['background/' + laser_delay_str][()], np.float64)
            #          print('Full signal trace', ts, 'with length', ts.shape, '\n')
            signal_trace = np.array(self.file['signal/' + laser_delay_str][()], np.float64)[gate_for_signal]
            try:
                background_trace = np.array(self.file['background/' + laser_delay_str][()], np.float64)[
                    gate_for_background]
            except IndexError:
                self.broken_background_trace = np.array(self.file['background/' + laser_delay_str][()], np.float64)
                print(
                    'AHIAHIAHI I have found an IndexError while reshaping the backgroudn signal!\nQuanti guai in Paraguay!!')
                print('Self.broken_background_trace is ', self.broken_background_trace, 'with length',
                      self.broken_background_trace.shape)
                print('Signal trace is ', signal_trace.shape, 'length gate_for_background', gate_for_background.shape)

            #            signal_np[i] = np.mean( signal_trace - background_trace, dtype=np.float64) #does not work, different length even with same integration window
            signal_np[i] = np.mean(signal_trace, dtype=np.float64)  # forcing float64 for more accuracy, see docs
            background_np[i] = np.mean(background_trace, dtype=np.float64)
            i += 1

        self.jitter_evaluation = False  # false to skip the jitter eval

        if self.jitter_evaluation:
            for i in range(len(self.jitter_arrival_times)):
                self.jitter_arrival_times[i] = self.lif_timebase[int(self.jitter_arrival_times[i])]

            plt.figure()
            plt.title('Jitter evaluation')
            plt.hist(self.jitter_arrival_times, bins=10)
            plt.xlabel('Arrival time wrt trigger (us)')
        return signal_np, background_np
        # note: I don not return self.signal_np but signal_np, s.t. the function
        # can be used without overwriting the variables of the dataset.

    def count_single_photons(self, verbose=False, my_prominence=140e-3):
        '''
        This function is meant to be used in the init instead of "integrated in given time window"
        One can also call it later tough, as it does not overwrite the self.signal_np and self.background_np by default.
        Those are overwritten in the constructor.
        Since 18th May 2021, the raw data is multiplied by the number of averages taken! This is quite a difference if
        one wants to compare different datasets taken with different averages. Is much better as, without introducing
        further noise (we checked this), results in having always the same real voltages on the vertical axis, which
        can be compared with the scope and among different averages.

        Parameters
        ----------
        verbose - print and plot like hell
        my_prominence - very delicate parameter for peak finding. Sets the vertical threshold. Is in VOLTS

        Returns
        -------
        signal np, background np

        TODO:
        - fix normalization of the photon counted - done, jsut decide which route to take
        - plot histogram with number of events/bin
        - plot histograms with loglogscale.
        '''
        # these are the same as usual function
        signal_np = np.zeros(self.scan_length, np.float64)  # proper init of empty signal and background
        background_np = np.zeros(self.scan_length, np.float64)
        # gating as above, maybe not a good idea
        gate_for_signal = np.logical_and(self.lif_timebase >= self.signal_integration_time_win[0],
                                         self.lif_timebase <= self.signal_integration_time_win[1])
        gate_for_background = np.logical_and(self.lif_timebase >= self.background_integration_time_win[0],
                                             self.lif_timebase <= self.background_integration_time_win[1])
        # two more timebase as above, maybe also not a good idea.
        self.gated_signal_timebase = self.lif_timebase[gate_for_signal]
        self.gated_background_timebase = self.lif_timebase[gate_for_background]

        # *** very important parameter for photon counting ***
        # my_prominence = 40e-3   # in Volts! Typ. a peak is 50 mV to 200 mV high
        my_width = 10 # in vector index units, typ. FWHM around 10 points
        take_avg_of_number_of_photons = False  # change between two differnt counting possibilities

        hist_of_all_peaks_s = np.zeros([0], np.float64)  # this keeps histogram of single photons HEIGHT on signal channel
        hist_of_all_peaks_b = np.zeros([0], np.float64)  # same but for background channel
        hist_arrival_times_s = np.zeros([0], np.float64)  # this keeps histrogram of the arrival times, for signal
        hist_arrival_times_b = np.zeros([0], np.float64)  # for background

        i = 0
        for laser_delay_str in self.file['background']:
            # I keep the raw and smoothed traces separate. When on verbose, I plot a comparison between some of them
            # *** Here I also mutliply it by the number of averages, to get real units (*self.averages) ***
            raw_signal_trace = np.array(self.file['signal/' + laser_delay_str][()], np.float64)[
                                   gate_for_signal] * self.averages
            raw_background_trace = np.array(self.file['background/' + laser_delay_str][()], np.float64)[
                                       gate_for_background] * self.averages
            signal_trace = gaussian_filter1d(raw_signal_trace, 2.5)  # smoothing
            background_trace = gaussian_filter1d(raw_background_trace, 2.5)

            [signal_peaks, _] = find_peaks(signal_trace, prominence=my_prominence,
                                           width=my_width)  # peak finding function call
            [background_peaks, _] = find_peaks(background_trace, prominence=my_prominence, width=my_width)
            if take_avg_of_number_of_photons:
                # route number 1: sum up all the voltages and, at the end, divide by mean photon height
                signal_np[i] = np.sum(signal_trace[
                                          signal_peaks])  # here I sum up all the peak intensities. This is not photon counting, but not too far
                background_np[i] = np.sum(background_trace[background_peaks])
            else:
                # route number 2: just count the peaks you found (probably correct one)
                signal_np[i] = len(signal_peaks)
                background_np[i] = len(background_peaks)

            # append the peaks I found to the total peak, to plot a hist
            hist_of_all_peaks_s = np.append(hist_of_all_peaks_s, signal_trace[signal_peaks])
            hist_of_all_peaks_b = np.append(hist_of_all_peaks_b, background_trace[background_peaks])
            # similar for arrival times. Watch out: you have to use the GATED timebase
            hist_arrival_times_s = np.append(hist_arrival_times_s, self.gated_signal_timebase[signal_peaks])
            hist_arrival_times_b = np.append(hist_arrival_times_b, self.gated_background_timebase[background_peaks])

            if verbose and i < 10:  # plot first i raw data
                plt.figure(6)
                plt.hist(self.gated_signal_timebase[signal_peaks], bins=18, alpha=.5)

                plt.figure(102)
                plt.plot(signal_trace * 1e3, ms=0, lw=2)  # signal trace
                plt.plot(signal_peaks, signal_trace[signal_peaks] * 1e3, 'k', ms=8, lw=0)
                plt.xlabel('Array index')  # cause the photon counting algorithm works with the index
                plt.ylabel('Voltage on the scope (mV)')

                plt.figure(101)
                print('Figure 1010 BROKEN, OVERWRITING ITSELF')
                plt.subplot(2, 2, 1)  # raw vs smooth
                plt.title('Comparison RAW vs SMOOTH, signal traces')
                plt.plot(self.gated_signal_timebase * 1e6, raw_signal_trace * 1e3)
                plt.plot(self.gated_signal_timebase * 1e6, signal_trace * 1e3)
                plt.xlabel('Time after trigger (us)')
                plt.ylabel('Voltage on scope (mV)')

                plt.subplot(2, 2, 2)  # Did you found the peaks in signal channel?
                plt.title('(Smoothed!) Signal traces with single photon counting')
                plt.plot(signal_trace * 1e3, ms = 0, lw = 2)  # signal trace
                plt.plot(signal_peaks, signal_trace[signal_peaks] * 1e3, 'k', ms=8, lw=0)
                plt.xlabel('Array index')  # cause the photon counting algorithm works with the index
                plt.ylabel('Voltage on the scope (mV)')
                #                plt.ylim([0, 5.e-2])

                plt.subplot(2, 2, 3)
                plt.title('Background traces with single photon counting')
                #               plt.ylim([0, 5.e-2])
                plt.plot(background_trace * 1e3)
                plt.plot(background_peaks, background_trace[background_peaks] * 1e3, 'k', ms=8, lw=0)
                plt.xlabel('Array index')  # cause the photon counting algorithm works with the index
                plt.ylabel('Voltage on the scope (mV)')
            i += 1
        # fix the normalization of the photons counted (IF route 1)
        print('Mean of a single photon signal \t\t+- std (mV): ', np.mean(hist_of_all_peaks_s) * 1e3, '\t+- ',
              np.std(hist_of_all_peaks_s) * 1e3)
        print('Mean of a single photon background \t+- std (mV): ', np.mean(hist_of_all_peaks_b) * 1e3, '\t+- ',
              np.std(hist_of_all_peaks_b) * 1e3)
        if take_avg_of_number_of_photons:
            # this gives, roughly, the signal expressed in single photons, which can be also non-integer
            signal_np = signal_np / np.mean(hist_of_all_peaks_s)
            background_np = background_np / np.mean(hist_of_all_peaks_b)

        # this gives the number of photons per shot, regardless of route 1 or 2
        signal_np = signal_np / self.averages
        background_np = background_np / self.averages
        # plot here what you just computed, instead of plotting it also later
        if verbose:
            print('Total number of photons in signal channel: \t', len(signal_peaks))
            print('Total number of photons in backgroundtt channel: \t', len(background_peaks))
            print(hist_of_all_peaks_s, len(hist_of_all_peaks_s))
            # figure with the scan, can be removed
            plt.figure()
            plt.plot(self.timebase_np * 1e6, signal_np)
            plt.plot(self.timebase_np * 1e6, background_np)
            plt.plot(self.timebase_np * 1e6, signal_np - background_np, lw=1)
            if self.smoothing_value != 0:
                plt.plot(self.timebase_np * 1e6, gaussian_filter1d(signal_np - background_np, self.smoothing_value),
                         'k', lw=4)
            else:
                plt.plot(self.timebase_np * 1e6, signal_np - background_np, 'k', lw=4)
            plt.xlabel('Laser delay (us)')
            plt.ylabel('Single photons (number, per laser shot)')
            legend = ['signal', 'background', 'difference', 'difference + smoothing']
            plt.legend(legend)
            plt.title('Signal and background with photon counting')

            # figure with histograms
            plt.figure()
            plt.subplot(1, 3, 1)  # plot the histogram of all photons
            # watch out: units are all in Volts, but for plots better mV
            plt.plot(my_prominence * 1e3, np.mean(signal_np), ms=10)
            plt.plot(1600, np.mean(signal_np), ms=10)  # scope theshold, no point above (on paper)
            plt.hist(hist_of_all_peaks_s * 1e3, bins=100, alpha=.5)
            plt.hist(hist_of_all_peaks_b * 1e3, bins=100, alpha=.5)
            legend = ['s.p. threshold', 'scope saturates at...', 'photons on signal channel',
                      'photons on background channel', ]
            plt.legend(legend)
            plt.xlabel('Height of single photon peak (mV)')
            plt.ylabel('number of photons')
            plt.title('Histogram of single photons detected')

            # plot the arrival times of all the photons - lin scale
            plt.subplot(1, 3, 2)
            plt.hist(hist_arrival_times_s * 1e6, bins=100, alpha=.6)
            plt.hist(hist_arrival_times_b * 1e6, bins=100, alpha=.6)
            legend = ['signal', 'background']
            plt.title('Arrival times of single photons - linlin')
            plt.xlabel('Arrival time (us)')
            plt.ylabel('Number of photons in the bin')
            plt.legend(legend)

            # plot the arrival times of all the photons - lin scale
            plt.subplot(1, 3, 3)
            plt.hist(hist_arrival_times_s * 1e6, bins=100, alpha=.6, log=True)
            plt.hist(hist_arrival_times_b * 1e6, bins=100, alpha=.6, log=True)
            legend = ['signal', 'background']
            plt.title('Arrival times of single photons - lin log')
            plt.xlabel('Arrival time (us)')
            plt.ylabel('Number of photons in the bin')
            plt.legend(legend)
        return signal_np, background_np

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

    def get_photodiode_signal(self, verbose=False):
        '''

        Args:
            verbose: plotting

        Returns:
            The two traces, to be saved than in the __init__
            Ugly, instead of another function this could just be the initi itself.
            The function can be called afterward to repeat the plotting e.g.
        '''
        photodiode_signal_np = np.zeros(self.scan_length, np.float64)  # proper init of empty signal and background
        photodiode_background_np = np.zeros(self.scan_length, np.float64)
        i = 0
        if verbose:
            plt.figure()
        for laser_delay_str in self.file['background']:
            photodiode_signal_trace = np.array(self.file['signal_pd/' + laser_delay_str][()], np.float64)
            photodiode_background_trace = np.array(self.file['background_pd/' + laser_delay_str][()], np.float64)

            photodiode_signal_np[i] = np.mean(photodiode_signal_trace,
                                              dtype=np.float64)  # forcing float64 for more accuracy, see docs
            photodiode_background_np[i] = np.mean(photodiode_background_trace, dtype=np.float64)
            i += 1
            if i < 4 and verbose:  # plot some traces
                plt.plot(photodiode_signal_trace)
                plt.plot(photodiode_background_trace)
        if verbose:
            plt.title('Some photodiode traces (unknown timebase currently)')
            plt.xlabel('Index of vector')
            plt.ylabel('Oscilloscope voltage (mV? probably, to check)')
            print('The photodiode trace is ', len(photodiode_signal_trace), ' points long')
            plt.figure()
            plt.plot(self.timebase_np * 1e3, photodiode_signal_np)
            plt.plot(self.timebase_np * 1e3, photodiode_background_np)
            plt.title('Averaged photdiode signal VS laser delay')
            plt.xlabel('Laser delay (ms)')
            plt.ylabel('Mean voltage (mV probably)')

        return photodiode_signal_np, photodiode_background_np

    def plot_signal_and_background(self, wavelength_scan=False):
        """
        By default plots the analyzed dataset.
        If not binning/smoothing/further analysis is being made, it will plot automatically the raw data.
        Added later option to change timebase and xlabel to correctly plot wavelength scans
        """
        plt.figure()
        legend = ['signal', 'background', 'difference']
        # smoothed/binned data
        if not wavelength_scan:
            # smoothed/binned data
            plt.plot(self.timebase_np * 1e6, self.analyzed_signal,
                     self.timebase_np * 1e6, self.analyzed_background,
                     self.timebase_np * 1e6, self.analyzed_signal - self.analyzed_background)
            # raw data
            plt.plot(self.timebase_np * 1e6, self.signal_np,
                     self.timebase_np * 1e6, self.background_np,
                     self.timebase_np * 1e6, self.signal_np - self.background_np, alpha=0.6, ms=1)
            # fill between the two
            plt.fill_between(self.timebase_np * 1e6, 0, self.analyzed_signal - self.analyzed_background, alpha=.3,
                             color='green', lw=2)
            plt.xlabel(self.xlabel)
        else:  # plot wavelength scan timebase
            ax = plt.gca()
            ax.format_coord = lambda x, y: '%8f, %8f' % (x, y)  # fix problem displaying digits...
            # smoothed/binned data
            plt.plot(self.timebase_np, self.analyzed_signal,
                     self.timebase_np, self.analyzed_background,
                     self.timebase_np, self.analyzed_signal - self.analyzed_background)
            # raw data
            plt.plot(self.timebase_np, self.signal_np,
                     self.timebase_np, self.background_np,
                     self.timebase_np, self.signal_np - self.background_np, alpha=0.6, ms=1)
            # fill between the two
            plt.fill_between(self.timebase_np, 0, self.analyzed_signal - self.analyzed_background, alpha=.3,
                             color='green', lw=2)
            plt.xlabel('Laser wavelength (cm**-1)')

        plt.title('file n. ' + str(self.file_numbers[0]))
        plt.legend(legend)
        plt.ylabel(self.ylabel)  # this change automatically depending on single photon mode
        jitter = False
        if jitter == True:
            plt.figure()
            plt.hist(self.analyzed_background, bins=15)

    def plot_lif_signal(self, print_mean_values=False):
        plt.figure()
        legend = ['signal', 'background', 'difference']
        plt.plot(self.lif_timebase * 1e6, self.lif_signal_np,
                 self.lif_timebase * 1e6, self.lif_background_np,
                 self.lif_timebase * 1e6, self.lif_difference_np, lw=.8, ms=0.2)
        plt.xlabel(self.xlabel)
        plt.ylabel('LIF signal')
        plt.legend(legend)
        plt.ylim([-0.01, 1])
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
            if excess != 0:
                print('Binning the dataset, excess is ', excess)
                # get a shorter vector but with mean values
                analyzed_signal = np.mean(
                    np.reshape(self.signal_np[:-excess], (len(self.signal_np) // bin_value, bin_value)), 1)
                analyzed_background = np.mean(
                    np.reshape(self.background_np[:-excess], (len(self.background_np) // bin_value, bin_value)), 1)
                # re-shape a vector with the same size at the self.signal_np
                analyzed_signal = np.append(np.repeat(analyzed_signal, bin_value), self.signal_np[-excess:])
                analyzed_background = np.append(np.repeat(analyzed_background, bin_value), self.background_np[-excess:])
            else:  # counting vector positions form 0 causes troubles...
                print('Binning the dataset, excess is zero')
                # get a shorter vector but with mean values
                analyzed_signal = np.mean(np.reshape(self.signal_np, (len(self.signal_np) // bin_value, bin_value)), 1)
                analyzed_background = np.mean(
                    np.reshape(self.background_np, (len(self.background_np) // bin_value, bin_value)), 1)
                # re-shape a vector with the same size at the self.signal_np
                analyzed_signal = np.repeat(analyzed_signal, bin_value)
                analyzed_background = np.repeat(analyzed_background, bin_value)

            print("Dataset binned by", bin_value, " bins")
        return analyzed_signal, analyzed_background

    def get_smoothed_dataset(self, smoothing_value):
        analyzed_signal, analyzed_background = self.signal_np, self.background_np  # default is non-processed dataset
        if smoothing_value > 0:
            print('Smoothing the dataset with smoothing_val of ', smoothing_value)
            analyzed_signal = gaussian_filter1d(analyzed_signal, smoothing_value)
            analyzed_background = gaussian_filter1d(analyzed_background, smoothing_value)
        return analyzed_signal, analyzed_background

    def find_peaks(self, peak_window_indexes=[44, 60], my_prominence=0, my_width=2, plotting=True, verbose=True):
        [signal_peaks, _] = find_peaks(self.analyzed_signal[peak_window_indexes[0]:peak_window_indexes[1]],
                                       prominence=my_prominence, width=my_width)
        signal_peaks = signal_peaks + peak_window_indexes[0]  # sum what was subtracted
        signal_peaks_heights = self.analyzed_signal[signal_peaks]
        [difference_peaks, _] = find_peaks(
            (self.analyzed_signal - self.analyzed_background)[peak_window_indexes[0]:peak_window_indexes[1]],
            prominence=my_prominence, width=my_width)
        difference_peaks = difference_peaks + peak_window_indexes[0]
        difference_peaks_heights = (self.analyzed_signal - self.analyzed_background)[difference_peaks]
        if plotting == True:
            plt.figure()
            plt.plot(self.analyzed_signal, alpha=0.5)
            plt.scatter(signal_peaks, self.analyzed_signal[signal_peaks], lw=5)
            plt.plot(self.analyzed_signal - self.analyzed_background, alpha=0.5)
            plt.scatter(difference_peaks, (self.analyzed_signal - self.analyzed_background)[difference_peaks], lw=3)
            plt.xlabel('array index (can be converted to laser delay)')
            plt.ylabel('signal intensity')
            plt.legend(['signal', 'signal-background'])
        if verbose:
            print('Peaks in signal:\t', signal_peaks, '\nwith values:\t\t', signal_peaks_heights)
            print('Peaks in difference:\t', difference_peaks, '\nwith values:\t\t', difference_peaks_heights)

        return signal_peaks, signal_peaks_heights, difference_peaks, difference_peaks_heights

    def fit_OH_lifetime(self, offset=False, weights=False):
        '''
        This functions make an exponential lifetime fit of a SingleDataset.
        The fit is made with x = self.timebase, y = self.analysed_signal - self.analysed_background
        Data is fitted to the function:
            amplitude * exp( -t /tau ) + offset
        The fit also make a plot where the 3 fit parameters (ampl, tau, offset) with uncertainties and the reduced
        chi squared are reported.
        I use the lmfit package. I have tried other packages like scilearnt and scipy.optimize, but they proved either |
        to give different results or to have too poor statistical analysis.
        I have tested Matlab against lmfit on one single exponential decay (as of July 2021) and they give the same
        numerical results, both for parameters and their uncertainties.

        One can add an automitized way to save the fit results in a txt file. This is very good and should be implemented.

        The weights are weight = 1/self.yerr, i.e. 1/sigma. Funny, I tough it should have been 1/sigma**2.
        '''
        print('Doing exponential decay fit.')
        from lmfit import Model, fit_report
        def exponential_decay_model(x, amplitude, tau, offset):
            return amplitude * np.exp(-x / tau) + offset

        my_model = Model(exponential_decay_model)
        print('Model: amplitude * exp(-t /tau) + offset')
        # these contains the initial parameters
        fit_params = my_model.make_params(amplitude=.5, tau=.5, offset=0, verbose=False,
                                          scale_covar=True)
        if not offset:
            print('Fitting forcing zero offset')
            fit_params['offset'].vary = False
            fit_params['offset'].value = 0.

        # fit_params['amplitude'].vary = False
        # fit_params['amplitude'].value = 0.522

        if weights:
            print('Fit done using weights, stored in self.yerr')
            result = my_model.fit(self.analyzed_signal - self.analyzed_background, fit_params, x=self.timebase_np,
                                  weights=1 / self.yerr, scale_covar=True)
            print('scale_covar set to False')
            # extra important option: scale_covar False or True. Most likely False.
        else:
            print('Fit done without weights, stored in self.yerr')
            result = my_model.fit(self.analyzed_signal - self.analyzed_background, fit_params, x=self.timebase_np)

        result.params.pretty_print()
        self.fitted_signal_minus_background = result.eval()

        # fancy plotting
        plt.figure()
        plt.errorbar(self.timebase_np, self.analyzed_signal - self.analyzed_background,
                     yerr=self.yerr, capsize=3, elinewidth=0.5, capthick=1.5, label='analysed signal-background',
                     lw=1.5)
        plt.plot(self.timebase_np, self.fitted_signal_minus_background, label='fit', lw=1.5)
        plt.plot(self.timebase_np, self.analyzed_signal, label='raw signal', lw=1.5)
        plt.plot(self.timebase_np, self.analyzed_background, label='raw background', lw=1.5)

        if offset:
            plt.plot(self.timebase_np, self.timebase_np * 0 + result.params['offset'].value, label='offset', ms=0)

        plt.legend()
        plt.xlabel('laser delay (s)')
        plt.ylabel('signal (photons/laser shot)')
        if offset:
            if weights:
                plt.title('Fitting with A * exp( -t /tau) + offset (w. weights)')

            else:
                plt.title('Fitting with A * exp( -t /tau) + offset (w.o. weights)')
        else:
            if weights:
                plt.title('Fitting with A * exp( -t /tau) and weights')
            else:
                plt.title('Fitting with A * exp( -t /tau) without weights')

        # add results to plot for better labbok
        plt.figtext(.15, .15, 'red. chis. = ' + str("{:.4f}".format(result.redchi)))
        plt.figtext(.15, .2, 'chi-sq. = ' + str("{:.3f}".format(result.chisqr)))
        plt.figtext(.15, .25, 'tau = ' + str("{:.3f}".format(result.params['tau'].value)) + ' s +- ' + str(
            "{:.3f}".format(result.params['tau'].stderr)) + 's')
        plt.figtext(.15, .3, 'A = ' + str("{:.3f}".format(result.params['amplitude'].value)) + ' s +- ' + str(
            "{:.3f}".format(result.params['amplitude'].stderr)) + 's')
        if offset:
            plt.figtext(.15, .35, 'offset = ' + str("{:.3f}".format(result.params['offset'].value)) + ' s +- ' + str(
                "{:.3f}".format(result.params['offset'].stderr)) + 's')
        print(fit_report(result))
        # saving the fit statistics to a txt file to be drag and drop in OneNote.
        # with open('temp_fit_result.txt', 'w') as filename:
        #     filename.write(fit_report(result))

    def fit_a_single_peak(self, inf_index=0, sup_index=-1, plotting=True):
        '''

        :param inf_index: position of the lower start of the peak
        :param sup_index: position of the upper start of the peak
        :param plotting:
        :return: nothing, just prints
        '''

        from lmfit.models import GaussianModel
        x = self.timebase_np[inf_index: sup_index]
        y = self.analyzed_signal[inf_index: sup_index] - self.analyzed_background[inf_index: sup_index]

        mod = GaussianModel()
        pars = mod.guess(y, x=x)
        out = mod.fit(y, pars, x=x)
        if plotting:
            plt.figure()
            plt.subplot(1, 3, 1)
            plt.plot(self.analyzed_signal-self.analyzed_background)
            plt.scatter(inf_index, (self.analyzed_signal-self.analyzed_background)[inf_index], color='red')
            if sup_index==-1:
                plt.scatter(len(self.analyzed_signal), (self.analyzed_signal-self.analyzed_background)[sup_index], color='red')
            else:
                plt.scatter(sup_index, (self.analyzed_signal-self.analyzed_background)[sup_index], color='red')

            plt.subplot(1, 3, 2)
            plt.plot(x, y, label='experim. peak')
            plt.plot(x, out.eval(), label='fit')
            plt.legend()
            plt.xlabel(self.xlabel)
            plt.ylabel(self.ylabel)

            plt.subplot(1, 3, 3)
            plt.plot(self.timebase_np, self.analyzed_signal-self.analyzed_background, label='experim. peak')
            plt.plot(self.timebase_np, out.eval(x=self.timebase_np), label='fit')
            plt.legend()
            plt.xlabel(self.xlabel)
            plt.ylabel(self.ylabel)


            print(out.fit_report(min_correl=0.25))

        print('AREA UNDER THE PEAK [ph/shot * micro seconds]:\t' , 1e6*out.best_values['amplitude'], '\t+- ', out.params['amplitude'].stderr)
        peak_height = out.best_values['amplitude']/(out.best_values['sigma'] * np.sqrt(2*np.pi))
        peak_stderr = np.sqrt((peak_height*out.params['amplitude'].stderr/out.params['amplitude'].value)**2+(peak_height*out.params['sigma'].stderr/out.params['amplitude'].value)**2)
        print('PEAK HEIGHT [ph/shot]\t\t\t\t\t\t\t', peak_height, '\t+- ', peak_stderr)
        out.conf_interval()

    def fancy_plotting(self):
        print("Still to be implemented")
        Milano_orange_lighter = [254 / 256, 156 / 256, 0 / 256]
        Milano_orange = [214 / 256, 106 / 256, 0 / 256]
        darker_gray = [80 / 256, 80 / 256, 80 / 256]
        light_gray = [128 / 256, 128 / 256, 128 / 256]
        fig = plt.figure()
        ax = plt.axes()
        # # smoothed/binned data
        # plt.plot(self.timebase_np * 1e6, self.analyzed_signal,
        #          self.timebase_np * 1e6, self.analyzed_background,
        #          self.timebase_np * 1e6, self.analyzed_signal - self.analyzed_background)
        plt.plot(self.timebase_np * 1e6, self.analyzed_signal, label='bla', color=Milano_orange)

        # fill between the two
        # plt.fill_between(self.timebase_np * 1e6, 0, self.analyzed_signal - self.analyzed_background, alpha=.3,
        #                  color='green', lw=2)
        plt.xlabel(self.xlabel, color=light_gray)
        plt.ylabel(self.ylabel, color=light_gray)

        plt.grid(None)
        plt.rc('axes', edgecolor=light_gray)
        plt.gca().get_xticklabels()[1].set_color(light_gray)
        plt.gca().yaxis.label.set_color(light_gray)
        plt.gca().xaxis.label.set_color(light_gray)
        plt.gca().tick_params(axis='x', labelcolor=light_gray)  # Works
        plt.gca().tick_params(axis='y', labelcolor=light_gray)  # Works
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.label.set_color(light_gray)
        ax.tick_params(axis='x', colors=light_gray)
        ax.yaxis.label.set_color(light_gray)
        ax.tick_params(axis='y', colors=light_gray)

        # legend = ['10 kV', '12.5 kV']
        # plt.legend(legend, frameon=False, labelcolor=light_gray)

        # plt.savefig('temp.pdf')


class MultiplesCoherentDataset(SingleDataset):
    '''
    This class collects a set of Single dataset.
    Works only if the timebase of the laser scan is exactly the same and is not supposed to work otherwise.
    It currently inherits from the SingleDataset class, but is probably useless/stupid. Few methods are overwritten.
    You can access the single i-th datasaet with self.complete_dataset[i]
    '''

    def __init__(self, path_of_dataset, first_file_number,
                 last_file_number, list_files_nums_to_skip=[],
                 signal_integration_time_win=default_sign_int_time_win,
                 background_integration_time_win=default_bckgrnd_int_time_win,
                 bin_value=0,  # bin value will be passed to SingleDataset constructor
                 smoothing_value=0,
                 photon_counting=False,
                 second_da_channel=False):
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
                                               background_integration_time_win, bin_value, smoothing_value,
                                               photon_counting, second_da_channel) for i in
                                 range(self.number_of_datasets)]
        # kaboom

        # this check if the timebase of the whole dataset is consistent.
        # if not, it will print out "AHIHAHIAHI" before Python crashes.
        self.is_timebase_coherent = self.check_timebase_coherence()  # just store the bool

        # data that pertains to the whole dataset
        '''
        s/b matrix raw: keeps all the signal_np and background_np. A bit useless imho.
        s/b matrix analysed: keeps all the analyzed_signal and analyzed_background
        ..._integrated_over_laser_delay (s and b): matrices integrated over laser delay
        ..._integrated over filename: matrices integrated over filename
        REMEMBER: ALL OF THESE ARE NEVER THE DIFFERENCE BETWEEN SIGNAL AND BACKGROUND!
        '''
        self.signal_matrix_raw = np.zeros([self.number_of_datasets, len(self.complete_dataset[0].signal_np)],
                                          np.float64)
        self.background_matrix_raw = np.zeros([self.number_of_datasets, len(self.complete_dataset[0].background_np)],
                                              np.float64)
        self.signal_matrix_analyzed = np.zeros([self.number_of_datasets, len(self.complete_dataset[0].signal_np)],
                                               np.float64)
        self.background_matrix_analyzed = np.zeros(
            [self.number_of_datasets, len(self.complete_dataset[0].background_np)],
            np.float64)
        for i in range(self.number_of_datasets):
            self.signal_matrix_raw[i, :] = self.complete_dataset[i].signal_np  # raw
            self.background_matrix_raw[i, :] = self.complete_dataset[i].background_np
            self.signal_matrix_analyzed[i, :] = self.complete_dataset[i].analyzed_signal  # analysed
            self.background_matrix_analyzed[i, :] = self.complete_dataset[i].analyzed_background

        # mean over x axis (laser delay) and over y axis (filename), of analyzed data
        # these being means, they are reliable to compare absolute values
        self.signal_avg_over_filenames = np.mean(self.signal_matrix_analyzed, axis=0)
        self.background_avg_over_filenames = np.mean(self.background_matrix_analyzed, axis=0)
        self.signal_avg_over_laser_delay = np.mean(self.signal_matrix_analyzed, axis=1)
        self.background_avg_over_laser_delay = np.mean(self.background_matrix_analyzed, axis=1)

        # in single photon counting mode, it makes sense to sum up all the photons, i.e. to integrate the signal
        # over laser delays to get the whole peak height
        # units will be # photons. (at variance with #photons/laser shot)
        self.signal_integrated_over_laser_delay = np.sum(self.signal_matrix_analyzed, axis=1)
        self.background_integrated_over_laser_delay = np.sum(self.background_matrix_analyzed, axis=1)

        # lif data
        self.lif_timebase = self.complete_dataset[0].lif_timebase  # re-cycle the timebase
        [self.lif_signal_np, self.lif_background_np, self.lif_difference_np] = self.get_lif_signal()

        # recycle the timebases
        self.timebase_np = self.complete_dataset[0].timebase_np
        self.gated_signal_timebase = self.complete_dataset[0].gated_signal_timebase
        self.gated_background_timebase = self.complete_dataset[0].gated_background_timebase


        if photon_counting:
            # This errors assume Poisson statistics.
            # sigma_rate = sqrt(rate) / sqrt(N_experiments)
            # such that sigma_rate / rate = 1 / sqrt( N_photons)
            self.yerr_signal = np.sqrt(self.signal_avg_over_filenames) / np.sqrt(
                self.complete_dataset[0].averages * self.number_of_datasets)
            self.yerr_background = np.sqrt(self.background_avg_over_filenames) / np.sqrt(
                self.complete_dataset[0].averages * self.number_of_datasets)
        else: # noise is given by intensity fluctuations across different repetitions
            self.yerr_signal = np.std(self.signal_matrix_analyzed, 0) / np.sqrt(self.number_of_datasets)
            self.yerr_background = np.std(self.background_matrix_analyzed, 0) / np.sqrt(self.number_of_datasets)
        # summed up independently
        self.yerr = np.sqrt(self.yerr_signal ** 2 + self.yerr_background ** 2) # error on difference is the quadrature

        # plotting
        self.xlabel = self.complete_dataset[0].xlabel  # very very ugly, my apologize
        self.ylabel = self.complete_dataset[0].ylabel
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

    def get_lif_signal(self):
        '''
        I inherit the method from single classes and just sum up everthing toghether.
        Could be moved to init
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
        for i in range(self.number_of_datasets):
            lif_signal_np += self.complete_dataset[i].lif_signal_np
            lif_background_np += self.complete_dataset[i].lif_background_np
        lif_difference_np = lif_signal_np - lif_background_np  # direct background subtraction
        return lif_signal_np, lif_background_np, lif_difference_np

    def plot_signal_and_background(self):
        # main plot
        plt.figure()
        plt.rcParams['font.size'] = '12'

        # 2D subplot
        plt.subplot(2, 1, 1)
        plt.imshow(self.signal_matrix_analyzed - self.background_matrix_analyzed, origin='lower', aspect='auto',
                   interpolation='none',
                   extent=[self.timebase_np[0] * 1e6,
                           self.timebase_np[-1] * 1e6,
                           self.file_numbers[0], self.file_numbers[-1]])
        ax = plt.gca()
        ax.set_yticks(np.arange(self.file_numbers[0], self.file_numbers[-1] + 1, 1))  # fix filename axis
        ax.set_yticklabels(np.arange(self.file_numbers[0], self.file_numbers[-1] + 1, 1))
        plt.colorbar(label=self.ylabel)
        plt.grid(None)
        plt.xlabel(self.xlabel), plt.ylabel(self.ylabel_2dscan)

        # plot signal, background and difference integrated over the filenames
        plt.subplot(2, 2, 3)
        # analyzed data
        plt.plot(self.timebase_np * 1e6, self.signal_avg_over_filenames,
                 self.timebase_np * 1e6, self.background_avg_over_filenames,
                 self.timebase_np * 1e6,
                 self.signal_avg_over_filenames - self.background_avg_over_filenames, lw=1, ms=0.3)
        plt.xlabel(self.xlabel), plt.ylabel(self.ylabel), plt.title('Signal integrated over all the files')
        legend = ['signal', 'background', 'difference']
        plt.legend(legend)
        # raw data
        plt.plot(self.timebase_np * 1e6, np.mean(self.signal_matrix_raw, axis=0),
                 self.timebase_np * 1e6, np.mean(self.background_matrix_raw, axis=0),
                 self.timebase_np * 1e6,
                 np.mean(self.signal_matrix_raw, axis=0) - np.mean(self.background_matrix_raw, axis=0), alpha=0.6, ms=1)
        # fill between the two
        plt.fill_between(self.timebase_np * 1e6, 0,
                         self.signal_avg_over_filenames - self.background_avg_over_filenames, alpha=.3,
                         color='green', lw=2)

        # plot the signal integrated over the laser delays vs filenames
        plt.subplot(2, 2, 4)
        legend = ['signal', 'background', 'difference']
        plt.plot(self.file_numbers, self.signal_integrated_over_laser_delay, lw=1)
        plt.plot(self.file_numbers, self.background_integrated_over_laser_delay)
        plt.plot(self.file_numbers,
                 self.signal_integrated_over_laser_delay - self.background_integrated_over_laser_delay)
        plt.xlabel(self.ylabel_2dscan), plt.ylabel('Total # photons'), plt.title(
            'Signal integrated over laser delay')
        plt.xticks(self.file_numbers)
        plt.legend(legend)

        plt.rcParams['font.size'] = '12'

        # plot all the traces together, separately. If offset==True, it adds an offset
        offset = False
        # colors = plt.cm.jet(np.linspace(0, 1, self.number_of_datasets))
        plt.figure()
        for i in range(self.number_of_datasets):
            if offset:
                plt.plot(self.timebase_np * 1e6,
                         self.signal_matrix_analyzed[i, :] - self.background_matrix_analyzed[i, :] + .33 * i * np.mean(
                             self.signal_avg_over_filenames),
                         alpha=.5, ms=3, lw=2)
            else:
                plt.plot(self.timebase_np * 1e6,
                         self.signal_matrix_analyzed[i, :] - self.background_matrix_analyzed[i, :],
                         alpha=.8, ms=3, lw=5)
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.title('(Signal - background) for all the filenames')
        legend = self.file_numbers
        plt.legend(legend)

        # plot all the traces AVERAGED together (both raw and analysed, both signal, background, difference)
        # useful when the signal is low
        plt.figure()
        # avg of analysed data
        plt.plot(self.timebase_np * 1e6, self.signal_avg_over_filenames,
                 self.timebase_np * 1e6, self.background_avg_over_filenames,
                 self.timebase_np * 1e6,
                 self.signal_avg_over_filenames - self.background_avg_over_filenames, lw=1, ms=0.3)
        plt.xlabel(self.xlabel), plt.ylabel(self.ylabel), plt.title(
            'Averaged of all the files, raw/analysed,\n signal/background/difference')
        legend = ['signal', 'background', 'difference']
        plt.legend(legend)
        # fill between the two
        plt.fill_between(self.timebase_np * 1e6, 0,
                         self.signal_avg_over_filenames - self.background_avg_over_filenames, alpha=.3,
                         color='green', lw=2)
        # avg of raw data
        plt.plot(self.timebase_np * 1e6, np.mean(self.signal_matrix_raw, axis=0),
                 self.timebase_np * 1e6, np.mean(self.background_matrix_raw, axis=0),
                 self.timebase_np * 1e6,
                 np.mean(self.signal_matrix_raw, axis=0) - np.mean(self.background_matrix_raw, axis=0), alpha=0.6, ms=1)

        # plot only the 2d image
        plt.figure()
        plt.imshow(self.signal_matrix_analyzed - self.background_matrix_analyzed, origin='lower', aspect='auto',
                   interpolation='none',
                   extent=[self.timebase_np[0] * 1e3,
                           self.timebase_np[-1] * 1e3,
                           self.file_numbers[0], self.file_numbers[-1]])
        ax = plt.gca()
        ax.set_yticks(np.arange(self.file_numbers[0], self.file_numbers[-1] + .3, 1))  # fix filename axis
        ax.set_yticklabels(np.arange(self.file_numbers[0], self.file_numbers[-1] + 1, 1))
        plt.colorbar(label=self.ylabel)
        plt.grid(None)
        plt.xlabel(self.xlabel), plt.ylabel(self.ylabel_2dscan)
        plt.title('Analysed (signal-background) plot of the 2D scan')

        # plot many additional informations in only one single plot
        additional_plots = True
        if additional_plots:
            add_info = plt.figure()
            # quick and dirty way to plot histogram of the background. Is it gaussian distributed?
            plt.rcParams['font.size'] = '12'
            plt.subplot(2, 3, 1)
            for i in range(self.number_of_datasets):
                plt.hist(self.background_matrix_raw[i, :], bins=10, alpha=0.4)
            plt.hist(np.reshape(self.background_matrix_raw, [self.background_matrix_raw.size, 1]),
                     bins=10 * self.number_of_datasets, alpha=0.6, fill='None')
            plt.title("Distribution of RAW background values,\n integrated in given int window")
            plt.xlabel(self.ylabel)
            plt.ylabel('Occurence')
            plt.grid(None)

            # quick and dirty plot the SNR against the filenumber
            plt.subplot(2, 3, 2)
            for i in range(self.number_of_datasets):
                plt.scatter(self.file_numbers[i], self.complete_dataset[i].SNR, color='magenta')
            plt.title('SNR vs filename number')
            plt.xlabel('Filename number')
            plt.ylabel('SNR')
            plt.xticks(self.file_numbers)

            plt.subplot(2, 3, 3)  # plot of mean background against file number
            plt.grid(True, which='major', lw=.2)
            for i in range(self.number_of_datasets):
                plt.scatter(self.file_numbers[i], np.mean(self.background_matrix_raw[i]), color='red')
            plt.title('Mean raw background vs filename')
            plt.xlabel('Filename number')
            plt.xticks(self.file_numbers)
            plt.ylabel(self.ylabel)

            # 3 image plot of signal, background-mean(background), signal-background.
            # Does the signal correlates with the background? If yes, may be due to OH molecules around
            plt.subplot(2, 3, 4)
            plt.imshow(self.signal_matrix_analyzed, origin='lower', aspect='auto',
                       interpolation='none',
                       extent=[self.timebase_np[0] * 1e3,
                               self.timebase_np[-1] * 1e3,
                               self.file_numbers[0], self.file_numbers[-1]])
            plt.title('Signal (all analysed)'), plt.colorbar()
            plt.xlabel('laser delay (us)'), plt.ylabel('filenumber')
            ax = plt.gca()
            ax.set_yticks(np.arange(self.file_numbers[0], self.file_numbers[-1] + .3, 1))  # fix filename axis
            ax.set_yticklabels(np.arange(self.file_numbers[0], self.file_numbers[-1] + 1, 1))
            plt.grid(None)

            plt.subplot(2, 3, 5)
            plt.imshow(self.background_matrix_analyzed, origin='lower',
                       aspect='auto',
                       interpolation='none', cmap='hot',
                       extent=[self.timebase_np[0] * 1e3,
                               self.timebase_np[-1] * 1e3,
                               self.file_numbers[0], self.file_numbers[-1]])
            plt.title('Background minus mean \nof ' + str(np.mean(self.background_matrix_raw))[:7]), plt.colorbar()
            plt.xlabel('laser delay (us)')
            plt.grid(None)
            ax = plt.gca()
            ax.set_yticks(np.arange(self.file_numbers[0], self.file_numbers[-1] + .3, 1))  # fix filename axis
            ax.set_yticklabels(np.arange(self.file_numbers[0], self.file_numbers[-1] + 1, 1))

            plt.subplot(2, 3, 6)
            plt.imshow(self.signal_matrix_analyzed - self.background_matrix_analyzed, origin='lower', aspect='auto',
                       interpolation='none',
                       extent=[self.timebase_np[0] * 1e3,
                               self.timebase_np[-1] * 1e3,
                               self.file_numbers[0], self.file_numbers[-1]])
            plt.title('Difference'), plt.xlabel('laser delay (us)')
            plt.colorbar(label=self.ylabel)
            plt.grid(None)

    def plot_lif_signal(self, print_mean_values=False):
        plt.figure()
        legend = ['signal', 'background', 'difference']
        plt.subplot(1, 2, 1)
        plt.plot(self.lif_timebase * 1e6, self.lif_signal_np, self.lif_timebase * 1e6, self.lif_background_np,
                 self.lif_timebase * 1e6, self.lif_difference_np, lw=.8, ms=0.2)
        plt.xlabel(self.complete_dataset[0].xlabel)
        plt.ylabel('LIF signal'), plt.title('LIF signal from the whole dataset - linlin')
        plt.legend(legend)
        plt.ylim([-0.2, np.max(self.lif_signal_np) * .1])
        plt.subplot(1, 2, 2)  # same but in loglog
        plt.loglog(self.lif_timebase * 1e6, self.lif_signal_np, self.lif_timebase * 1e6, self.lif_background_np,
                   self.lif_timebase * 1e6, self.lif_difference_np, lw=.8, ms=0.2)
        plt.xlabel(self.complete_dataset[0].xlabel + ' (log)')
        plt.ylabel('LIF signal - log'), plt.title('LIF signal from the whole dataset -loglog')
        plt.legend(legend)

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

    def plot_repetitions(self):
        # two plots: the standard one, with yerr  + a big one with many more details about underlaying statistics
        plt.rcParams['font.size'] = '12'
        colors = plt.cm.jet(np.linspace(0, 1, self.number_of_datasets))
        # for brevity
        sig_minus_back_int_ov_a_fl = self.signal_avg_over_filenames - self.background_avg_over_filenames

        plt.figure()  # first standard plot
        # analyzed data
        plt.plot(self.timebase_np * 1e6, self.signal_avg_over_filenames,
                 self.timebase_np * 1e6, self.background_avg_over_filenames,
                 self.timebase_np * 1e6, sig_minus_back_int_ov_a_fl,
                 lw=3, ms=0.3)
        # fill an area of +- yerr
        plt.fill_between(self.timebase_np * 1e6, sig_minus_back_int_ov_a_fl - self.yerr,
                         sig_minus_back_int_ov_a_fl + self.yerr, alpha=.4,
                         color='green', lw=1)
        plt.fill_between(self.timebase_np * 1e6, 0., sig_minus_back_int_ov_a_fl, alpha=.3, lw=0)
        plt.plot(self.timebase_np * 1e6, self.timebase_np * 0., 'k', lw=1, alpha=.5, ms=0,
                 label='_nolegend_')  # zero baseline
        plt.xlabel(self.xlabel), plt.ylabel(self.ylabel), plt.title(
            'Signal averaged over all the files\n+- possion error ( = sqrt(rate) / (number of experiment) )')
        legend = ['signal', 'background', 'difference', '+- self.yerr']
        plt.legend(legend)

        # second plot
        plt.figure()
        plt.subplot(1, 2, 1)
        # analyzed data
        plt.plot(self.timebase_np * 1e6, self.signal_avg_over_filenames,
                 self.timebase_np * 1e6, self.background_avg_over_filenames,
                 self.timebase_np * 1e6, sig_minus_back_int_ov_a_fl,
                 lw=3, ms=0.3)
        # fill an area of +- 1 std * prefactor
        plt.fill_between(self.timebase_np * 1e6, sig_minus_back_int_ov_a_fl - self.yerr,
                         sig_minus_back_int_ov_a_fl + self.yerr, alpha=.4,
                         color='green', lw=1)
        # without prefactor, just +- 1std
        plt.fill_between(self.timebase_np * 1e6, sig_minus_back_int_ov_a_fl - self.yerr,
                         sig_minus_back_int_ov_a_fl + self.yerr, alpha=.2,
                         lw=0)
        plt.plot(self.timebase_np * 1e6, self.timebase_np * 0., 'k', lw=1, alpha=.5, ms=0,
                 label='_nolegend_')  # zero baseline
        plt.xlabel(self.xlabel), plt.ylabel(self.ylabel), plt.title('Signal averaged over all the files\n+- yerr')
        legend = ['signal', 'background', 'difference', '+- yerr ', '+- std']
        plt.legend(legend)

        # plot all the traces + (thick, red) the average
        plt.subplot(2, 2, 2)
        for i in range(self.number_of_datasets):
            plt.plot(self.timebase_np * 1e6, self.signal_matrix_analyzed[i, :] - self.background_matrix_analyzed[i, :],
                     'k', lw=1.2, ms=.01, alpha=.5, color=colors[i])
        plt.plot(self.timebase_np * 1e6, self.signal_avg_over_filenames - self.background_avg_over_filenames, 'red',
                 lw=2.8, ms=0)
        plt.xlabel(self.xlabel), plt.ylabel(self.ylabel), plt.title('Analysed signal-background, all the files + avg')

        cumulative_diff = np.zeros([0, ], np.float64)  # save all the difference to plot histograms
        plt.subplot(2, 4, 7)  # plot a trace
        for i in range(self.number_of_datasets):
            difference = self.signal_matrix_analyzed[i, :].T - self.background_matrix_analyzed[i, :].T \
                         - self.signal_avg_over_filenames + self.background_avg_over_filenames
            plt.plot(self.timebase_np * 1e6, difference, 'k', lw=1.2, ms=.01, alpha=.5, color=colors[i])
            cumulative_diff = np.append(cumulative_diff, difference)

        plt.plot(self.timebase_np * 1e6, 0. * self.timebase_np, 'red', lw=2.8, ms=0)
        plt.xlabel(self.xlabel), plt.ylabel(self.ylabel), plt.title('analysed (signal-background)\n - mean of dataset')

        plt.subplot(2, 4, 8)  # plot histograms
        plt.hist(np.reshape(cumulative_diff, [len(self.timebase_np), self.number_of_datasets]), bins=20, stacked=True,
                 alpha=.7, orientation='horizontal')
        plt.hist(cumulative_diff, bins=20, histtype='step', lw=3, edgecolor='black', orientation='horizontal')
        plt.xlabel('occurrences'), plt.title('Distributions of deviation from the mean')

    def fit_OH_lifetime(self, offset=False):
        '''
        See fit_OH_lifetime function in SingleDataset for header
        '''
        print('Doing exponential decay fit.')
        from lmfit import Model, fit_report
        def exponential_decay_model(x, amplitude, tau, offset):
            return amplitude * np.exp(-x / tau) + offset

        my_model = Model(exponential_decay_model)
        print('Model: amplitude * exp(-t /tau) + offset')
        # these contains the initial parameters
        fit_params = my_model.make_params(amplitude=.5, tau=.5, offset=0., verbose=False,
                                          scale_covar=True)
        if not offset:
            print('Fitting forcing zero offset')
            fit_params['offset'].vary = False
            fit_params['offset'].value = 0.
        #        fit_params['amplitude'].vary = False
        #        fit_params['amplitude'].value = .45

        result = my_model.fit(self.signal_avg_over_filenames - self.background_avg_over_filenames, fit_params,
                              x=self.timebase_np,
                              weights=1 / self.yerr)

        result.params.pretty_print()  # print everything with little effort.
        self.fitted_signal_minus_background = result.eval()

        # fancy plotting
        plt.figure()
        plt.errorbar(self.timebase_np, self.signal_avg_over_filenames - self.background_avg_over_filenames,
                     yerr=self.yerr, capsize=3, elinewidth=0.5, capthick=1.5, label='analysed signal-background')
        plt.plot(self.timebase_np, self.fitted_signal_minus_background, label='fit')
        if offset:  # plot horizontal line with the fitted offset
            plt.plot(self.timebase_np, self.timebase_np * 0 + result.params['offset'].value, label='offset', ms=0)

        plt.legend()
        plt.xlabel('laser delay (s)')
        plt.ylabel('signal (photons/laser shot)')
        if not offset:
            plt.title('Fitting with A * exp( -t /tau)')
        else:
            plt.title('Fitting with A * exp( -t /tau) + offset')

        # add results to plot for better labbook
        plt.figtext(.15, .15, 'red. ch.s.' + str("{:.4f}".format(result.redchi)))
        plt.figtext(.15, .2, 'chi-sq. = ' + str("{:.3f}".format(result.chisqr)))
        plt.figtext(.15, .25, 'tau = ' + str("{:.3f}".format(result.params['tau'].value)) + ' s +- ' + str(
            "{:.3f}".format(result.params['tau'].stderr)) + 's')
        plt.figtext(.15, .3, 'A = ' + str("{:.3f}".format(result.params['amplitude'].value)) + ' s +- ' + str(
            "{:.3f}".format(result.params['amplitude'].stderr)) + 's')
        plt.figtext(.15, .35, 'offset = ' + str("{:.3f}".format(result.params['offset'].value)) + ' s +- ' + str(
            "{:.3f}".format(result.params['offset'].stderr)) + 's')
        print(fit_report(result))
        # saving the fit statistics to a txt file to be drag and drop in OneNote.
        # with open('temp_fit_result.txt', 'w') as filename:
        #     filename.write(fit_report(result))

    def fit_a_single_peak(self, inf_index=0, sup_index=-1, plotting=True, weights=False):
        '''

        :param inf_index: position of the lower start of the peak
        :param sup_index: position of the upper start of the peak
        :param plotting:
        :return: nothing, just prints
        '''

        from lmfit.models import GaussianModel
        x = self.timebase_np[inf_index: sup_index]
        y = self.signal_avg_over_filenames[inf_index: sup_index] - self.background_avg_over_filenames[inf_index: sup_index]
        weight = 1/self.yerr[inf_index: sup_index]

        mod = GaussianModel()
        pars = mod.guess(y, x=x)

        if weights:
            print('Fit done using weights, stored in self.yerr')
            out = mod.fit(y, pars, x=x, weights=weight, scale_covar=False)
            print('scale_covar set to False')
            # extra important option: scale_covar False or True. Most likely False.
        else:
            out = mod.fit(y, pars, x=x)

        if plotting:
            plt.figure()
            plt.subplot(1, 3, 1)
            plt.plot(self.signal_avg_over_filenames-self.background_avg_over_filenames)
            plt.scatter(inf_index, (self.signal_avg_over_filenames-self.background_avg_over_filenames)[inf_index], color='red')
            if sup_index==-1:
                plt.scatter(len(self.signal_avg_over_filenames), (self.signal_avg_over_filenames-self.background_avg_over_filenames)[sup_index], color='red')
            else:
                plt.scatter(sup_index, (self.signal_avg_over_filenames-self.background_avg_over_filenames)[sup_index], color='red')

            plt.subplot(1, 3, 2)
            plt.plot(x, y, label='experim. peak')
            plt.plot(x, out.eval(), label='fit')
            plt.legend()
            plt.xlabel(self.xlabel)
            plt.ylabel(self.ylabel)

            plt.subplot(1, 3, 3)
            plt.plot(self.timebase_np, self.signal_avg_over_filenames-self.background_avg_over_filenames, label='experim. peak')
            plt.plot(self.timebase_np, out.eval(x=self.timebase_np), label='fit')
            plt.legend()
            plt.xlabel(self.xlabel)
            plt.ylabel(self.ylabel)

            plt.figure()
            plt.plot(self.timebase_np*1e6, self.signal_avg_over_filenames-self.background_avg_over_filenames, color='indianred', label='exp. data')
            plt.plot(self.timebase_np*1e6, out.eval(x=self.timebase_np), color='silver', label='gaussian fit', ms=1,)
            plt.legend()
            plt.xlabel('time of flight (\u03bcs)')
            plt.ylabel(self.ylabel)

            print(out.fit_report(min_correl=0.25))

        peak_height = out.best_values['amplitude'] / (out.best_values['sigma'] * np.sqrt(2 * np.pi))
        peak_stderr = np.sqrt((peak_height * out.params['amplitude'].stderr / out.params['amplitude'].value) ** 2 + (peak_height * out.params['sigma'].stderr / out.params['sigma'].value) ** 2)

        print('AREA\t\t\t', 'AREA ERROR\t\t', 'HEIGHT\t\t\t', 'HEIGHT ERROR')
        print(1e6 * out.best_values['amplitude'], 1e6 * out.params['amplitude'].stderr, peak_height, peak_stderr)

        # print('AREA UNDER THE PEAK [ph/shot * micro seconds]:\t' , 1e6*out.best_values['amplitude'], '\t+- ', 1e6*out.params['amplitude'].stderr)
        # peak_height = out.best_values['amplitude']/(out.best_values['sigma'] * np.sqrt(2*np.pi))
        # peak_stderr = np.sqrt((peak_height*out.params['amplitude'].stderr/out.params['amplitude'].value)**2+(peak_height*out.params['sigma'].stderr/out.params['amplitude'].value)**2)
        # print('PEAK HEIGHT [ph/shot]\t\t\t\t\t\t\t', peak_height, '\t+- ', peak_stderr)
        out.conf_interval()


class FakeDataset():
    # TODO: implement it maybe on day
    def __init__(self, path_of_dataset, file_number,
                 signal_integration_time_win=default_sign_int_time_win,
                 background_integration_time_win=default_bckgrnd_int_time_win,
                 bin_value=0,
                 smoothing_value=0):  # bin value will be passed to SingleDataset constructor
        print('__init__ of fake dataset')

        # retrieve the proper file structure
        [self.files_in_folder_full_path, _, _, self.file_numbers] = get_path_of_one_dataset(path_of_dataset,
                                                                                            file_number, file_number,
                                                                                            debug_mode=False,
                                                                                            list_files_nums_to_skip=[])

        self.file = h5py.File(self.files_in_folder_full_path[0], 'r')

        # attrbutes of this h5 file are saved in self.attributes
        self.attributes = self.file.attrs
        self.number_of_averages = self.attributes['averages']  # a bit useless
        self.relchannel = self.attributes['relchannel']
        self.relchannel = self.attributes['scanchannel']
        # note to myself: the bgchannel is always 4 even if the real bg channel is set to, let's say, 1 ...
        # never never trust
        #        self.bg_channel = self.attributes['bgchannel'] #bg channel, can be usefull
        print(list(self.attributes))
        # FakeScan has only these attributes:
        # ['averages', 'bgchannel', 'relchannel', 'scanchannel']
        # Missing: osci_tbins, osci_tstart, osci_tstop, scan0, scan0_time, scan1, scan1_time, scanchannel, scanrel0, scanrel1

        # to extract each of them, use:
        #            self.bgchannel = self.file['inputs/bgchannel'][()]


#        self.scan_length = len(list(self.file['signal'].keys()))
#        self.signal_integration_time_win = signal_integration_time_win
#        self.background_integration_time_win = background_integration_time_win
#        self.inputs = self.file['inputs']

# each dataset has this structure:
#    'background'      (no attrs)
#         '0.0023'     (shape 25002, attributes x1, x2, xb1, xb2, y, y_oc)
# they are 'background', 'background_pd', 'signal', 'signal_pd'

# the timebase is the same for all the 4 data super-sets: signal, signal_pd, background, background_pd
#        print(max(self.file['signal'].keys()))
#        self.timebase_np = np.array(list(self.file['signal'].keys()),
#                                   dtype=np.float64)  # stored in a proper numpy array


if __name__ == "__main__":
    plt.close('all')
    # path_of_dataset = r'C:\Users\Gruppe Willitsch\Desktop\Lab315\Lab315\DA_scripts\da_hyt_two\test_dataset\2021-03-25-043.h5'
    # path_of_dataset = r'/home/pietro/Documents/PhD/coding/Lab315/Lab315/DA_scripts/da_hyt_two/test_dataset/2021-03-25-043.h5'
    path_of_dataset = r'/Users/pietro/Documents/PhD/coding/Lab315/Lab315/DA_scripts/da_hyt_two/test_dataset/2021-01-07-045.h5'
    tt = MultiplesCoherentDataset(path_of_dataset, 40, 41, smoothing_value=1., photon_counting=True,
                                  second_da_channel=False)
    tt.plot_signal_and_background()
    tt.plot_repetitions()
    plt.show()
