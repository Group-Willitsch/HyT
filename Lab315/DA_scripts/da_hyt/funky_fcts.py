import numpy as np
import sys
import matplotlib.pyplot as plt
import h5py
import traceback
from scipy.signal import argrelextrema
import os


def bug():
    plt.show()
    sys.exit()


def smooth_1d_fct(data, smooth_val):
    from scipy.ndimage import gaussian_filter1d
    data = gaussian_filter1d(data, smooth_val)
    return (data)


def smooth_2d_row(array_2d, gating, smooth_val, norm):
    from scipy.ndimage import gaussian_filter1d

    if gating:
        print('NOTE: ALL NAN valvues in the data set will be set to 0!!!')
        print('NOTE: ALL NAN valvues in the data set will be set to 0!!!')
        print('NOTE: ALL NAN valvues in the data set will be set to 0!!!')
        array_2d[np.isnan(array_2d)] = 0
    if norm:
        for i in range(len(array_2d[:, 0])):
            array_2d[i, :] = gaussian_filter1d(array_2d[i, :], smooth_val)
            array_2d[i, :] = array_2d[i, :] / np.max(array_2d[i, :])
    else:
        try:
            for i in range(len(array_2d[:, 0])):
                array_2d[i, :] = gaussian_filter1d(array_2d[i, :], smooth_val)
        except:
            from colorama import Fore
            from colorama import Style
            print(f'{Fore.RED}You have only one scan but try a 2D smooth... fool...{Style.RESET_ALL}')
            bug()
    return (array_2d)


class get_data:
    def __init__(self, file_list):
        self.file_list = file_list
        self.laser = False
        self.lif_gate_sig = [2E-6, 100E-6]  # [0.4E-6, 4.8E-6]#
        self.lif_gate_bac = [0E-6, 2E-7]
        self.test_pc = False

        ###--- Local variables for collecting the data
        self.y_arr = 0
        self.signal = 0
        self.back = 0
        self.laser_arr_sig = 0
        self.laser_arr_bac = 0
        self.scan1_t_att = 0
        self.scan2_t_att = 0
        self.check_t1 = 0
        self.check_t2 = 0
        self.h5_file = 0
        self.check_att = 0

    def collect_only_data(self, i, y_s, y_b, laser_sig, laser_bac):
        if i == 0:
            self.y_arr = y_s - (y_b)
            self.signal = y_s
            self.back = y_b

        else:
            self.y_arr = np.vstack([self.y_arr, (y_s * 1 - (y_b) * 1)])
            self.signal = np.vstack([self.signal, y_s])
            self.back = np.vstack([self.back, y_b])

    def classic_only_data(self):
        for i in range(len(self.file_list)):
            self.h5_file = h5py.File(self.file_list[i], 'r')
            scan_var = list(self.h5_file['signal'].keys())

            ###---Check scan parameters---
            self.check_att = self.h5_file.attrs

            y_s = np.zeros((len(scan_var)))
            y_b = np.zeros((len(scan_var)))
            laser_sig = np.zeros((len(scan_var)))
            laser_bac = np.zeros((len(scan_var)))

            ####-----Get the y-values for the data-----####
            for j in range(len(scan_var)):
                file_attrs = self.h5_file['signal/' + str(scan_var[j])].attrs
                y_s[j] = file_attrs['y_oc']
                file_attrs = self.h5_file['background/' + str(scan_var[j])].attrs
                y_b[j] = file_attrs['y_oc']

                if self.laser:
                    file_attrs = self.h5_file['signal_pd/' + str(scan_var[j])].attrs
                    laser_sig[j] = file_attrs['ypd']
                    file_attrs = self.h5_file['background_pd/' + str(scan_var[j])].attrs
                    laser_bac[j] = file_attrs['ypd']

            get_data.collect_only_data(self, i, y_s, y_b, laser_sig, laser_bac)
        return (self.y_arr, self.signal, self.back)

    def collect(self, i, y_s, y_b, laser_sig, laser_bac):
        if i == 0:
            self.y_arr = y_s - (y_b)
            self.signal = y_s
            self.back = y_b
            if self.laser:
                self.laser_arr_sig = laser_sig
                self.laser_arr_bac = laser_bac

            ###---Check scan parameters---
            self.scan1_t_att = self.check_att['scan0_time']
            self.scan2_t_att = self.check_att['scan1_time']
            try:
                self.check_t1 = self.h5_file['inputs/time1'].value[3]
                self.check_t2 = self.h5_file['inputs/time2'].value[3]
            except:
                try:
                    self.check_t1 = self.h5_file['inputs/self.guithread.time1'].value[3]
                    self.check_t2 = self.h5_file['inputs/self.guithread.time2'].value[3]
                except:
                    check_t1 = self.h5_file['inputs/read.time1'].value[3]
                    check_t2 = self.h5_file['inputs/read.time2'].value[3]

        else:
            self.y_arr = np.vstack([self.y_arr, (y_s * 1 - (y_b) * 1)])
            self.signal = np.vstack([self.signal, y_s])
            self.back = np.vstack([self.back, y_b])

            if self.laser:
                self.laser_arr_sig = np.vstack([self.laser_arr_sig, laser_sig])
                self.laser_arr_bac = np.vstack([self.laser_arr_bac, laser_bac])

            ###---Check scan parameters---
            self.scan1_t_att = np.vstack([self.scan1_t_att, self.check_att['scan0_time']])
            self.scan2_t_att = np.vstack([self.scan2_t_att, self.check_att['scan1_time']])
            try:
                self.check_t1 = np.vstack([self.check_t1, self.h5_file['inputs/time1'].value[3]])
                self.check_t2 = np.vstack([self.check_t2, self.h5_file['inputs/time2'].value[3]])
            except:
                try:
                    self.check_t1 = np.vstack([self.check_t1, self.h5_file['inputs/self.guithread.time1'].value[3]])
                    self.check_t2 = np.vstack([self.check_t2, self.h5_file['inputs/self.guithread.time2'].value[3]])
                except:
                    self.check_t1 = np.vstack([self.check_t1, self.h5_file['inputs/read.time1'].value[3]])
                    self.check_t2 = np.vstack([self.check_t2, self.h5_file['inputs/read.time1'].value[3]])

    def classic(self):
        for i in range(len(self.file_list)):
            self.h5_file = h5py.File(self.file_list[i], 'r')
            scan_var = list(self.h5_file['signal'].keys())

            ###---Check scan parameters---
            self.check_att = self.h5_file.attrs

            y_s = np.zeros((len(scan_var)))
            y_b = np.zeros((len(scan_var)))
            laser_sig = np.zeros((len(scan_var)))
            laser_bac = np.zeros((len(scan_var)))

            ####-----Get the y-values for the data-----####
            for j in range(len(scan_var)):
                file_attrs = self.h5_file['signal/' + str(scan_var[j])].attrs
                y_s[j] = file_attrs['y_oc']
                file_attrs = self.h5_file['background/' + str(scan_var[j])].attrs
                y_b[j] = file_attrs['y_oc']

                if self.laser:
                    file_attrs = self.h5_file['signal_pd/' + str(scan_var[j])].attrs
                    laser_sig[j] = file_attrs['ypd']
                    file_attrs = self.h5_file['background_pd/' + str(scan_var[j])].attrs
                    laser_bac[j] = file_attrs['ypd']

            get_data.collect(self, i, y_s, y_b, laser_sig, laser_bac)

        scan_var = [float(i) for i in scan_var]  ####---Convert to float, note there is still a problem with rounding..
        return (self.y_arr, self.signal, self.back, self.laser_arr_sig, self.laser_arr_bac, scan_var, self.scan1_t_att,
                self.scan2_t_att, self.check_t1, self.check_t2)

    def int_lif(self):
        for i in range(len(self.file_list)):
            self.h5_file = h5py.File(self.file_list[i], 'r')
            scan_var = list(self.h5_file['signal'].keys())

            ###---Check scan parameters---
            self.check_att = self.h5_file.attrs

            y_s = np.zeros((len(scan_var)))
            y_b = np.zeros((len(scan_var)))
            laser_sig = np.zeros((len(scan_var)))
            laser_bac = np.zeros((len(scan_var)))

            ####-----Get the t-values for the data-----####
            lif_time = np.linspace(self.check_att['osci_tstart'], self.check_att['osci_tstop'],
                                   self.check_att['osci_tbins'])
            lif_time = np.round(lif_time, 10)
            gate = np.logical_and(lif_time > self.lif_gate_sig[0], lif_time < self.lif_gate_sig[1])
            gate_offset = np.logical_and(lif_time > self.lif_gate_bac[0], lif_time < self.lif_gate_bac[1])

            ####-----Get the y-values for the data-----####
            for j in range(len(scan_var)):
                signal_trace = self.h5_file['signal/' + str(scan_var[j])][:]
                back_trace = self.h5_file['background/' + str(scan_var[j])][:]
                if False:  # debugging purposes
                    plt.figure(100)
                    print('gate:\t\t', gate)
                    print('gate_offset:\t\t', gate_offset)
                    plt.plot(lif_time, signal_trace)  # Check again for jitter!!!
                try:
                    # original from Thomas
                    y_s[j] = np.mean(signal_trace[gate]) - np.mean(signal_trace[gate_offset])
                    y_b[j] = np.mean(back_trace[gate]) - np.mean(back_trace[gate_offset])
                #                    y_s[j]       = np.mean(signal_trace[gate])
                #                    y_b[j]       = np.mean(back_trace[gate])

                except:
                    y_s[j] = np.nan
                    y_b[j] = np.nan

                if self.laser:
                    file_attrs = self.h5_file['signal_pd/' + str(scan_var[j])].attrs
                    laser_sig[j] = file_attrs['ypd']
                    file_attrs = self.h5_file['background_pd/' + str(scan_var[j])].attrs
                    laser_bac[j] = file_attrs['ypd']

            get_data.collect(self, i, y_s, y_b, laser_sig, laser_bac)

        scan_var = [float(i) for i in scan_var]  ####---Convert to float, note there is still a problem with rounding..
        return (self.y_arr, self.signal, self.back, self.laser_arr_sig, self.laser_arr_bac, scan_var, self.scan1_t_att,
                self.scan2_t_att, self.check_t1, self.check_t2)

    def pho_cnt(self):
        from scipy.ndimage import gaussian_filter1d
        from scipy.signal import argrelextrema

        for i in range(len(self.file_list)):
            self.h5_file = h5py.File(self.file_list[i], 'r')
            scan_var = list(self.h5_file['signal'].keys())

            ###---Check scan parameters---
            self.check_att = self.h5_file.attrs

            y_s = np.zeros((len(scan_var)))
            y_b = np.zeros((len(scan_var)))
            laser_sig = np.zeros((len(scan_var)))
            laser_bac = np.zeros((len(scan_var)))

            ###--- Adjust photon amplitude to the nbr of avg
            ###See OneNote/Schmit/Thomas/LIF signal analysis/1.3.2019 Single photon analysis
            photon_low = 0.19 / self.check_att['averages']
            photon_sig = 0.21 / self.check_att['averages']

            ####-----Get the t-values for the data-----####
            lif_time = np.linspace(self.check_att['osci_tstart'], self.check_att['osci_tstop'],
                                   self.check_att['osci_tbins'])
            lif_time = np.round(lif_time, 10)
            gate = np.logical_and(lif_time > self.lif_gate_sig[0], lif_time < self.lif_gate_sig[1])
            gate_offset = np.logical_and(lif_time > self.lif_gate_bac[0], lif_time < self.lif_gate_bac[1])

            ####-----Get the y-values for the data-----####
            for j in range(len(scan_var)):
                try:
                    signal_trace = self.h5_file['signal/' + str(scan_var[j])][:]
                    back_trace = self.h5_file['background/' + str(scan_var[j])][:]
                    if self.laser:
                        file_attrs = self.h5_file['signal_pd/' + str(scan_var[j])].attrs
                        laser_sig[j] = file_attrs['ypd']
                        file_attrs = self.h5_file['background_pd/' + str(scan_var[j])].attrs
                        laser_bac[j] = file_attrs['ypd']

                    ##lif_gate_bac

                    signal_gf = gaussian_filter1d(signal_trace[gate] - np.mean(signal_trace[gate_offset]), 5)
                    backgr_gf = gaussian_filter1d(back_trace[gate] - np.mean(back_trace[gate_offset]), 5)

                    ##--Remove values below threshold
                    signal_gf = signal_gf * (signal_gf > photon_low).astype(int)
                    backgr_gf = backgr_gf * (backgr_gf > photon_low).astype(int)

                    ##--Mark peaks
                    # argrelextrema
                    sig_peaks = np.array(argrelextrema(signal_gf, np.greater))
                    back_peaks = np.array(argrelextrema(backgr_gf, np.greater))

                    ##--Get amplitudes and peaks
                    sig_amps = signal_gf[sig_peaks]
                    back_amps = backgr_gf[back_peaks]

                    test_pc = False
                    if test_pc:
                        print('Test photon counting')
                        fgc += 1
                        gs = plt.GridSpec(1, 2)
                        fig = plt.figure(fgc)
                        ax1 = fig.add_subplot(gs[0, 0])
                        plt.title('Signal')
                        plt.plot(signal_gf)
                        plt.scatter(sig_peaks, sig_amps)
                        ax1 = fig.add_subplot(gs[0, 1])
                        plt.title('Background')
                        plt.plot(backgr_gf)
                        plt.scatter(back_peaks, back_amps)
                        ff.bug()

                    y_s[j] = np.sum(np.round(sig_amps / photon_sig))
                    y_b[j] = np.sum(np.round(back_amps / photon_sig))




                except:
                    print('Except')
                    traceback.print_exc()
                    y_s[j] = np.nan
                    y_b[j] = np.nan

            get_data.collect(self, i, y_s, y_b, laser_sig, laser_bac)

        scan_var = [float(i) for i in scan_var]  ####---Convert to float, note there is still a problem with rounding..
        return (self.y_arr, self.signal, self.back, self.laser_arr_sig, self.laser_arr_bac, scan_var, self.scan1_t_att,
                self.scan2_t_att, self.check_t1, self.check_t2)


class sort_laser():
    def __init__(self, laser_sig, laser_bac):
        self.laser_sig = laser_sig
        self.laser_bac = laser_bac
        self.laser_t = np.zeros(len(laser_bac.ravel()) * 2)

    def sort(self):

        laser_sig_1d = self.laser_sig.ravel()
        laser_bac_1d = self.laser_bac.ravel()

        for i in range(len(self.laser_t)):
            if i % 2 == 0:
                self.laser_t[i] = laser_sig_1d[int(i / 2)]
            if i % 2 == 1:
                self.laser_t[i] = laser_bac_1d[int(i / 2)]  # laser_bac_1d[int(i/2)]

        return (self.laser_t)


class gating_two_d():
    def __init__(self, y_arr, sig, bac, l_sig, l_bac, laser_t, sigma_g):
        self.y_arr = y_arr
        self.sig = sig
        self.bac = bac
        self.l_sig = l_sig
        self.l_bac = l_bac
        self.sigma_g = sigma_g
        self.laser_t = laser_t

    def gate_laser(self):
        laser_mean = np.nanmean(self.laser_t)
        laser_std = np.nanstd(self.laser_t)

        laser_min = laser_mean - self.sigma_g * laser_std
        laser_max = laser_mean + self.sigma_g * laser_std

        gate_laser_sig = (self.l_sig > laser_min) * (self.l_sig < laser_max)
        gate_laser_bac = (self.l_bac > laser_min) * (self.l_bac < laser_max)
        gate_laser_t = (self.laser_t > laser_min) * (self.laser_t < laser_max)

        return (gate_laser_sig, gate_laser_bac, gate_laser_t)

    def gate_data(self):
        ###--- Get boarders for the signal measurment
        data_mean_s = np.nanmean(self.sig, axis=0)
        data_std_s = np.nanstd(self.sig, axis=0)
        data_min_s = data_mean_s - self.sigma_g * data_std_s
        data_max_s = data_mean_s + self.sigma_g * data_std_s

        ###--- Get boarders for the background measumernt
        data_mean_b = np.nanmean(self.bac, axis=0)
        data_std_b = np.nanstd(self.bac, axis=0)
        data_min_b = data_mean_b - self.sigma_g * data_std_b
        data_max_b = data_mean_b + self.sigma_g * data_std_b

        ###--- Create 2D array for gating
        gate_sig = np.zeros(np.shape(self.sig))
        gate_bac = np.zeros(np.shape(self.sig))

        for i in range(len(gate_sig[0, :])):
            gate_sig[:, i] = (self.sig[:, i] > data_min_s[i]) * (self.sig[:, i] < data_max_s[i])
            gate_bac[:, i] = (self.bac[:, i] > data_min_b[i]) * (self.bac[:, i] < data_max_b[i])

        gate_sig[gate_sig == 0] = np.nan
        gate_bac[gate_bac == 0] = np.nan

        return self.y_arr * gate_sig * gate_bac, self.sig * gate_sig, self.bac * gate_bac


class analyse():
    def __init__(self, y_arr):
        self.y_arr = y_arr
        self.ente = 6

    def get_ente(self):
        print('Extract each ' + str(self.ente) + ' shot')

        nte_arr = np.zeros((self.ente, len(self.y_arr[0, :])))
        nte_arr_err = np.zeros((self.ente, len(self.y_arr[0, :])))
        for i in range(self.ente):
            nte_arr[i] = np.nanmean(self.y_arr[i::self.ente, :], axis=0)
            # print(str(i)+"x Eukalyptus")
            # print(nte_arr)
            nte_arr_err[i] = np.nanstd(self.y_arr[i::self.ente, :], axis=0) / np.sqrt(len(self.y_arr[i::3]))

        return (nte_arr, nte_arr_err)


class functions():
    def __init__(self):
        self.x = 0
        self.a = 0
        self.b = 0
        self.c = 0

    def exp(self):
        return self.a + self.b * np.exp(-self.x / self.c)


class fit():

    def __init__(self, x, y):
        from lmfit import Model
        self.x = x
        self.y = y
        self.func = functions()

    ###Fits
    def expfit(self):
        from lmfit import Model
        gmodel = Model()
        result = gmodel.fit(self.y, x=self.x, a=1e-5, b=1e-5, c=3e5)
        # print(result.best_values)
        # print('1/2 pressure is equal to '+str(np.round(result.best_values['c']/86400./1.35914091423,2))+' days')
        # print('1/e pressure is equal to '+str(np.round(result.best_values['c']/86400.,2))+' days')
        # return(result.best_values['a'],result.best_values['b'],result.best_values['c'])


def get_dataset_full_path(path_of_dataset, first_file_number, last_file_number, debug_mode=False,
                          list_files_nums_to_skip=[]):
    '''
    Function that gives you the full path of each data you wanna load.
    Takes into account:
        -existance of dirname_pressure_dataset
        -proper number position
        -workd for one file or one big dataset (in only one dirname_pressure_dataset)
        -How to do if you have more than one different dirname_pressure_dataset (for example for
                                                a overnight measurment?)

    Parameters
    ----------
    path_of_dataset : list with path to the wanted dataset or to the wanted dirname_pressure_dataset with the dataset in it
    first_file_number: smallerst filename, number between 1 and 999, extrema included, uint
    last_file_number: biggest filename, number between 1 and 999, extrema included, uint. Can be the same as firdt file num
    

    Returns
    -------
    files_in_folder_full_path :what originally called "files"
    dirname_dataset: name of directory where this dataset is placed
    files_in_folder: name of files loaded 
    '''
    dirname_dataset = os.path.dirname(path_of_dataset)  # this gives same output
    # for two different inputs, either the path of the .h5 file or the dirname_pressure_dataset just above
    # will give the same output.
    files_in_folder, files_in_folder_full_path, data_number = [], [], []
    for file_cycler in sorted(os.listdir(dirname_dataset)):
        if file_cycler.endswith(".h5"):
            temp, _ = os.path.splitext(file_cycler)  # temp stores the int number of file
            temp = int(temp[-3:])
            if not (temp in list_files_nums_to_skip):  # removes files to be skipped
                data_number.append(temp)  # this get the number of the dataset
                files_in_folder.append(file_cycler)

    if debug_mode:
        print("Debug mode in file loading")
        print('first_file_number: \t\t', first_file_number)
        print('last_file_number: \t\t', last_file_number)
        print('path_of_dataset:\t\t', path_of_dataset)
        print('dirname_dataset:\t\t', dirname_dataset)
        print('files_in_folder:\t\t', files_in_folder)
        print('data_number:\t\t', data_number)

    if os.path.exists(dirname_dataset) == False:
        print('Non existing path:\t\t:', path_of_dataset)
        print('Given path does not exists (or not given access to). \nClosing.\r\n')
        bug()

    # follows ugly piece of code to resize files_in_folder to the given limits
    data_number = np.array(data_number)  # ugly
    if first_file_number in data_number and last_file_number in data_number and first_file_number <= last_file_number and first_file_number >= 0 and last_file_number >= 1:
        temp = np.array(files_in_folder)
        files_in_folder = list(temp[np.logical_and(data_number >= first_file_number, data_number <= last_file_number)])
    else:
        print(
            'Given file numbers do not match the file numbers I found in this dirname_pressure_dataset, or some other trivial error. Scemo.\nQuit.\n')
        bug()
    if debug_mode:
        print('dirname_dataset:\t\t', dirname_dataset)
        print('Files I found:\t\t', data_number)
        print('Files founded: \t\t\t\t', files_in_folder)

    for file_cycler in files_in_folder:
        files_in_folder_full_path.append(os.path.join(dirname_dataset, file_cycler))

        # on Win can help:
    #    print('Any change? (check on Win)\t\t', os.path.normpath(path_of_dataset)) #check if helps on win
    print('Given dataset path:\t\t\t', path_of_dataset)
    print('Dataset directory exists:\t\t', dirname_dataset)

    return files_in_folder_full_path, dirname_dataset, files_in_folder


if __name__ == "__main__":
    x = np.linspace(0, 5, 6)

#    y = np.empty(np.shape(x))
#    for i in range(len(y)):
#        y[i] = len(y)-np.exp(i/len(y))
#    da = fit(x,y)
#    
#    a,b,c = da.expfit()


# print(da.)
