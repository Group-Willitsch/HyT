####-----------This is ur config file to analyse the data
####-----------U can do single files or multi-file analysis
import os 
import funky_fcts as ff



    
first_file_number = 21
last_file_number = 21
path_of_dataset = r'/Users/pietro/Documents/PhD/code_scripts/Lab315/Lab315/DA_scripts/da_hyt/test_dataset/'
path_of_dataset = r'C:\Users\Gruppe Willitsch\Desktop\data\2019-06-21/2019-06-21-001.h5'
[files, _, _] = ff.get_dataset_full_path( path_of_dataset, first_file_number, last_file_number, False, list_files_nums_to_skip=[])


start = first_file_number #looks like word start is used somewhere else...
end = last_file_number

files_mode  = 2       #### 0 = certain file(s),
                      #### 1 = files between two values, i.e., 1-10,
                      #### 2 = Generate ur own list......:(

#####------- Laser and parameters
laser  = False     ###Did you record the laser power on a photodiode or similar?

####----- Scan parameters
fake_scan   = False      ###Did you do a fake scan? Scan, for instance,x times the same delay with 10 averages before moving to the next delay


####----- Analysis parameters parameters
####----- Analysis parameters parameters
da_method   = 1     #### 0 = classic, 1 = variable gates, 2 = photon counting, 666 = test mode, plot directly the data
plot_format = 0       #### 0 = plot all for, for example, search of the incoupling time,
                      #### 1 = add everything like in lifetime measurments
plot_pietro = True

plot_ente   = 0
entee       = 0
gif         = False

two_d_plots = False   #### 2D analyses of the data


binning     = False
bin_val     = 2
smooth_1d   = False
smooth_2d   = False   ### Smooth the rows in 2d
smooth_val  = 2.5
norm        = False

#####------- gating
gating   = False
sigma_g  = 5        ####For the laser and data
gate_laser = False
gate_data =  False





#####------- If you want single shots from a 2D plot
shots_from_2d = [3,4, 5, 6, 7]#,'104','105']
sumgate = [1700, 2100] # for performance plot
lif_gate_sig = [0E-7, 5.0E-6]#[0.4E-6, 4.8E-6]#[0.4E-6, 2.53E-6]#
lif_gate_bac = [0.0E-6, .0E-7] # fuer offset correction


colors = ['b','g','r','c','m','y','k','w']









if __name__ == "__main__":
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
    files_nums_to_skip = [47, 48]
    [files_in_folder_full_path, dirname_dataset, files_in_folder] = ff.get_dataset_full_path(path_of_dataset, 
                                    first_file_number, last_file_number, debug_mode=False, 
                                    list_files_nums_to_skip=files_nums_to_skip)
    