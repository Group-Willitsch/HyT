import numpy as np
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import dates
from datetime import datetime
from scipy.optimize import curve_fit
from lmfit import Model
import glob
import pandas as pd
import os

'''
For other users: Ctrl+F "CHANGE" to see all the stuff you have to change to make it work for you pressures.
'''

folder_of_pressure_dataset =r'C:\Users\Gruppe Willitsch\Desktop\data\pressure_log'
last_N_days_to_plot = 10

plt.close('all')


def get_list_of_pressure_files(dirname_pressure_dataset):
    """
    CHANGE DIRNAME_PRESSURE_DATASET TO MATCH YOU FOLDER, the rest should work like a charm if you pressure log files
    are named like pressure_2021-11-07.txt. It will fail if the filenames have different patterns.
    List all the files ending in .txt and containing pressure_ in the file name.
    Get the full path for each of the file, the filename only and a tensor containing the day, month and year of each
    file.
    Files are sorted (fast and not messy).

    input:
    directory name of the pressure dataset, i.e. the single folder in which all the log files are kept.
    The dafault one is the one currently in use in Lab 3.15

    @return
        files_in_folder: list containting all the file names
        files_in_folder_full_path: same but with dirname_pressure_dataset in front
        dd_mm_yyyy: N*3 tensor containing dat, month, year for each dataset present in files_in_folder. Sorted, 1:1
        correspondence.
    """
    # default dirname_pressure_dataset is the correct one
    print('Folder is:', dirname_pressure_dataset)

    files_in_folder, files_in_folder_full_path = [], []
    dd_mm_yyyy = np.zeros((0, 3), dtype=int)  # store a N * 3 tensor to keep day, month, year
    for file_cycler in sorted(os.listdir(dirname_pressure_dataset)):  # sorted very imp. for speed
        if file_cycler.endswith(".txt") and 'pressure_' in file_cycler:  # hopefully get only pressure logs
            files_in_folder.append(file_cycler)
            temp, _ = os.path.splitext(file_cycler)  # temp stores the int number of file
            day, month, year = int(temp[-2:]), int(temp[-5:-3]), int(temp[-10:-6])
            dd_mm_yyyy = np.append(dd_mm_yyyy, [np.array([day, month, year])], axis=0)  # very ugly but works

    for file_cycler in files_in_folder:
        files_in_folder_full_path.append(os.path.join(dirname_pressure_dataset, file_cycler))

    #    print('files_in_folder', files_in_folder, '\n files_in_folder_full_path', files_in_folder_full_path, '\n data_number', dd_mm_yyyy)
    return files_in_folder_full_path, files_in_folder, dd_mm_yyyy


class pressure_dataset():

    def __init__(self, folder=r'C:\Users\Gruppe Willitsch\Desktop\data\pressure_log',
                 last_N_days = 5):
        #        print(files_in_folder_full_path)
        self.last_N_days = last_N_days
        [self.files_in_folder_full_path, self.files_in_folder, self.dd_mm_yyyy] = get_list_of_pressure_files(folder)
        self.files_in_folder_full_path = self.files_in_folder_full_path[-last_N_days:]
        print(self.files_in_folder_full_path)
        i = 0
        for file_cycler in self.files_in_folder_full_path:
            if i == 0:
                self.data = pd.read_csv(file_cycler, infer_datetime_format=True, delimiter='\t', header=None,
                                        skiprows=1)
            else:
                data_tmp = pd.read_csv(file_cycler, infer_datetime_format=True, delimiter='\t', header=None, skiprows=1)
                self.data = pd.DataFrame.append(self.data, data_tmp, ignore_index=True)
            i += 1
        print(self.data[0])
        self.data[0] = self.data[0].str.strip('W. Europe Standard')
        self.data[0] = pd.to_datetime(self.data[0], format='%Y-%m-%d %H:%M:%S', errors='ignore')

        print('here', self.data[0])

        # CHANGE HERE IF YOU HAVE DIFFERENT NUMBERS OF CHANNELS/ DIFFERENT NAMES
        self.timebase_str = np.array((self.data[0]).astype(str))
        print('timebase_str', self.timebase_str)
        self.source = np.array((self.data[1]).astype(float))
        self.dec = np.array((self.data[2]).astype(float))
        self.trap = np.array((self.data[3]).astype(float))
        self.for_dec = np.array((self.data[4]).astype(float))
        self.for_tra = np.array((self.data[5]).astype(float))

        # # Get relative time x in seconds
        # self.relative_times = np.zeros(len(self.data[0]))
        # for i in range(1, len(self.data[0])):
        #     self.relative_times[i] = self.relative_times[i - 1] + (
        #                 pd.Timestamp(self.data[0][i]) - pd.Timestamp(self.data[0][i - 1])).total_seconds()
        # print('relative_times', self.relative_times)

        # start_time = pd.Timedelta
        # df['totalseconds'] =

        # plotting
        self.xlabel = 'time (a.u.)'
        self.ylabel = 'pressure (mbar)'

    def plot_all(self):
        plt.figure()
        plt.title('Chambers pressures over last ' + str(self.last_N_days) + ' days')
        plt.subplot(2, 2, 1)
        plt.plot(self.source)
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.title('source chamber')

        plt.subplot(2, 2, 2)
        plt.plot(self.dec)
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.title('decelerator chamber')

        plt.subplot(2, 2, 3)
        plt.plot(self.trap)
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.title('trap chamber')

        plt.subplot(2, 2, 4)
        plt.plot(self.for_dec)
        plt.plot(self.for_tra)
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.title('forlines')

    # locs, labels = plt.xticks()
        # near_t = []
        #
        # for i in range(len(locs)):
        #     variable = (np.abs(dt - locs[i])).argmin()
        #     if i == 0:
        #         near_t = str(t[variable])
        #     else:
        #         near_t = np.append(near_t,str(t[variable]))
        #
        # ###--- Finalize plot
        # plt.xticks(locs[::2],near_t[::2],rotation = 45)
        # plt.grid()
        # plt.legend()
        # plt.tight_layout(pad=1.08, h_pad=None, w_pad=None, rect=None)
        # plt.yscale('Log')
        #


if __name__ == "__main__":
    ####------Parameters-----
    ####------ All parameters im config file
    bound = last_N_days_to_plot * (24 * 60 ** 2)
    my_pressure_dataset = pressure_dataset(folder_of_pressure_dataset, last_N_days=last_N_days_to_plot)
    my_pressure_dataset.plot_all()

#
#
#
#
#
#
#
#
# def expfit(x,a,b,c):
# 	return abs(a)+b*np.exp(-x/c)
#
# def expfit2para(x,b,c):
# 	return b*np.exp(-x/c)
#
# def get_data(files,folder):
#     ###4 Matter 4 time, #https://www.tutorialspoint.com/python/time_strptime.htm
#     formatter = mdates.DateFormatter('%Y-%m-%d %H')
#
#     print()
#
#
#
#     for i in range(len(files)):
#         file = str(folder[0])+str(files[i])
#         if i == 0:
#             data = pd.read_csv(file,infer_datetime_format=True, delimiter = '\t', header = None, skiprows = 1)
#         else:
#             data_tmp = pd.read_csv(file,infer_datetime_format=True, delimiter = '\t', header = None, skiprows = 1)
#             data = pd.DataFrame.append(data,data_tmp,ignore_index=True)
#
#
#
#     data[0] = data[0].str.strip('W. Europe Standard')
#     data[0] = pd.to_datetime(data[0], format='%Y-%m-%d %H:%M:%S', errors='ignore')
#
#     source   = np.array((data[1]).astype(float))
#     dec      = np.array((data[2]).astype(float))
#     trap     = np.array((data[3]).astype(float))
#     for_dec  = np.array((data[4]).astype(float))
#     for_tra  = np.array((data[5]).astype(float))
#
#     ####--- Get relative time x in seconds
#     r_time = np.zeros(len(data[0]))
#     for i in range(1,len(data[0])):
#         r_time[i] = r_time[i-1]+(pd.Timestamp(data[0][i])-pd.Timestamp(data[0][i-1])).total_seconds()
#     t = np.array(data[0].astype(str))
#
#     return(t,r_time,source,dec,trap,for_dec,for_tra)
#
# def gate_data(t,dt,source,dec,trap,for_dec,for_tra,gate):
#     t       = t[gate]
#     dt      = dt[gate]
#     dt      = dt-np.min(dt)
#     source  = source[gate]
#     dec     = dec[gate]
#     trap    = trap[gate]
#     for_dec = for_dec[gate]
#     for_tra = for_tra[gate]
#     return(t,dt,source,dec,trap,for_dec,for_tra)
#
#
#


#
# ##-- Get all data
# # t,dt,source,dec,trap,for_dec,for_tra = get_data(files,dirname_pressure_dataset)
#
# fgc+=1
# plt.figure(fgc)
#
# ###-- Gate all data
# if plot_t == 0:
#     print('Plot all times')
# if plot_t == 1:
#     gate = (dt>(np.max(dt)-bound))
#     t,dt,source,dec,trap,for_dec,for_tra = gate_data(t,dt,source,dec,trap,for_dec,for_tra,gate)
# if plot_t == 2:
#     gate = (dt<(np.max(dt)-bound))
#     t,dt,source,dec,trap,for_dec,for_tra = gate_data(t,dt,source,dec,trap,for_dec,for_tra,gate)
#
# ###-- Plot the selected channels
# for i in range(len(channels)):
#     if channels[i]=='source'  or channels[i] == '0':
#         data = source
#     if channels[i]=='dec'     or channels[i] == '1':
#         data = dec
#     if channels[i]=='trap'    or channels[i] == '2':
#         data = trap
#     if channels[i]=='for_src' or channels[i] == '3':
#         data = for_dec
#     if channels[i]=='for_dec' or channels[i] == '4':
#         data = for_tra
#
#     ##-- Fit data
#     if fit:
#         gmodel = Model(expfit2para)
#         #result = gmodel.fit(data, x=dt, a = 1e-11,b = 1e-5,c = 1e2)
#         result = gmodel.fit(data, x=dt,b = 1e-5,c = 1e2)
#
#         print(result.best_values)
#         print('1/2 pressure is equal to '+str(np.round(result.best_values['c']/86400./1.35914091423,2))+' days')
#         print('1/e pressure is equal to '+str(np.round(result.best_values['c']/86400.,2))+' days')
#         #params = gmodel.make_params()
#         ##result.best_values
#         ##result.best_fit
#
#     ##-- Plot data
#     plt.plot(dt,data,label=str(channels[i]))#,marker='o')
#     #plt.xlim(0,)
#     plt.ylabel('Pressure (mBar)')
#
#     ##-- If fit plot fit
#     if fit:
#         x = np.linspace(np.min(dt),np.max(dt)*3,1000)
#         #a = result.best_values['a']
#         b = result.best_values['b']
#         c = result.best_values['c']
#         #plt.plot(x,expfit(x,a,b,c),label='Final pressure is '+str(np.round(a*1e9,5))+str(' e-9 mBar'),alpha=0.5)
#         plt.plot(x,expfit2para(x,b,c),label='Final pressure is 0',alpha=0.5)
#
#
# ###--- Get the right ticks
# locs, labels = plt.xticks()
# near_t = []
#
# for i in range(len(locs)):
#     variable = (np.abs(dt - locs[i])).argmin()
#     if i == 0:
#         near_t = str(t[variable])
#     else:
#         near_t = np.append(near_t,str(t[variable]))
#
# ###--- Finalize plot
# plt.xticks(locs[::2],near_t[::2],rotation = 45)
# plt.grid()
# plt.legend()
# plt.tight_layout(pad=1.08, h_pad=None, w_pad=None, rect=None)
# plt.yscale('Log')
#
#
# plt.show()
#
#
#
#
#
