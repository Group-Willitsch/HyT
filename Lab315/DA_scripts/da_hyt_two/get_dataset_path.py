#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 10:06:00 2021

@author: pietro
"""
import numpy as np
import os
import sys


def get_path_of_one_dataset(path_of_dataset, first_file_number, last_file_number,
                          debug_mode=False, list_files_nums_to_skip=[]):
    '''
    class that contains the full path of each data you wanna load.
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
    

    Variables
    -------
    files_in_folder_full_path: what originally called "files"
    dirname_dataset: name of directory where this dataset is placed
    files_in_folder: name of files loaded 

    Args:
        debug_mode:
    '''
    dirname_dataset = os.path.dirname(path_of_dataset) #this gives same output
    #for two different inputs, either the path of the .h5 file or the dirname_pressure_dataset just above
    #will give the same output. NOPE does not work currently.
    files_in_folder, files_in_folder_full_path, data_number = [], [], []
    for file_cycler in sorted(os.listdir(dirname_dataset)):
        if file_cycler.endswith(".h5"):
            temp, _ = os.path.splitext(file_cycler) #temp stores the int number of file
            temp = int(temp[-3:])
            if not(temp in list_files_nums_to_skip): #removes files to be skipped
                data_number.append(temp) #this get the number of the dataset
                files_in_folder.append(file_cycler)


    if debug_mode:
        print("Debug mode in file loading")
        print('first_file_number: \t\t', first_file_number)
        print('last_file_number: \t\t', last_file_number)
        print('path_of_dataset:\t\t', path_of_dataset)
        print('dirname_dataset:\t\t', dirname_dataset)
        print('files_in_folder:\t\t', files_in_folder)
        print('data_number:\t\t', data_number)


    if os.path.exists(dirname_dataset)==False:
        print('Non existing path:\t\t:', path_of_dataset)
        print('Given path does not exists (or not given access to). \nClosing.\r\n')
        sys.exit()




    #follows ugly piece of code to resize files_in_folder to the given limits
    data_number = np.array(data_number) #ugly
    if first_file_number in data_number and last_file_number in data_number and first_file_number <= last_file_number and first_file_number >=0 and last_file_number>=1:
        temp = np.array(files_in_folder)
        files_in_folder = list(temp[ np.logical_and(data_number >= first_file_number, data_number <= last_file_number) ] )
        data_number = list(data_number[np.logical_and(data_number >= first_file_number, data_number <= last_file_number)] )
        #I also update data_number, cause keeping alredy the number is useful for plotting
    else:
        print('Given file numbers do not match the file numbers I found in this dirname_pressure_dataset, or some other trivial error. Scemo.\nQuit.\n')
        sys.exit()
    if debug_mode:
        print('dirname_dataset:\t\t', dirname_dataset)
        print('Files I found:\t\t', data_number)
        print('Files founded: \t\t\t\t', files_in_folder)

    for file_cycler in files_in_folder:
        files_in_folder_full_path.append( os.path.join(dirname_dataset, file_cycler) )

    #on Win can help:
#    print('Any change? (check on Win)\t\t', os.path.normpath(path_of_dataset)) #check if helps on win
    if debug_mode:
        print('Given dataset path:\t\t\t', path_of_dataset)
        print('Dataset directory exists:\t', dirname_dataset)
    return files_in_folder_full_path, dirname_dataset, files_in_folder, data_number


if __name__ == "__main__":
    first_file_number = 41
    last_file_number = 44
    files_nums_to_skip = []
    path_of_dataset = r'/Users/pietro/Documents/PhD/coding/Lab315/Lab315/DA_scripts/da_hyt_two/test_dataset/2021-01-07-041.h5'
    path_of_dataset = r'C:\Users\Gruppe Willitsch\Desktop\Lab315\Lab315\DA_scripts\da_hyt_two\test_dataset\2021-01-07-040.h5'


    [files_in_folder_full_path, dirname_dataset, files_in_folder, data_number] = get_path_of_one_dataset(path_of_dataset,
                                    first_file_number, last_file_number, debug_mode=True,
                                    list_files_nums_to_skip=files_nums_to_skip)
    verbose = False
    if verbose:
        print("Running config file as main\r\n")
        print('getcwd:      ', os.getcwd())
        print('__file__:    ', __file__)
        print('basename:    ', os.path.basename(__file__))
        print('dirname:     ', os.path.dirname(__file__))
        print('[set target path 1]')
        target_path_1 = os.path.join(os.path.dirname(__file__), 'target_1.txt')
        tt = os.path.join( os.path.dirname(__file__), 'test_dataset')
        print('target_path_1: ', target_path_1)
        print('normalize    : ', os.path.normpath(target_path_1))
        print(os.sep)
        print(os.sep is os.path.sep)

        print('basename of dataset:\t\t\t', os.path.basename(path_of_dataset) )
        basename_without_ext = os.path.splitext(os.path.basename(path_of_dataset))[0]
        print('basename without extension:\t', basename_without_ext)
        print(os.path.dirname(path_of_dataset))

        subdirname = os.path.basename(os.path.dirname(path_of_dataset))
        print(subdirname)





