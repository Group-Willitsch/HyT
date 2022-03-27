%% Script to read all the out_a*.dat files, load them correctly and re-save
% them as .mat files, which occupies much less space (16 Mb shrinks to 1
% Mb) and should be much faster to load too. 

%% NOTE: THE SYMMETRIZATION OF THE FIELDS IS NOT DONE IN THIS SCRIPT!!
% YOU MUST DO IT IN THE INPUTpARAMETERS CLASS (JUST TO BE SURE, so you can change it later on)

% based on the loadAccelerationFields function of the InputParameter class

%% normal mode fields
voltages = [10, 12.5];
SIMION_ni = 151; % ni number of grid points x (along beam axis
SIMION_nj = 41; % number of grid points y,z (perpendicular to beam axis)
SIMION_nk = 41; % number of grid points y,z (perpendicular to beam axis)

[ax_norm, ay_norm, az_norm] = deal([],[],[]);

fprintf('Loading normal mode fields ...\t')
for i =1:2
filename = '../dec_Norm_' + strrep( string(voltages(i)), '.', 'p') + 'kV/output/'
paren = @(x, varargin) x(varargin{:}); % in-line function to reshape a matrix without a temporary variable

ax_norm = permute( reshape( table2array( paren( ...
    readtable(filename + 'outax.dat'), ':', 4) ), ...
    [SIMION_nj, SIMION_nk, SIMION_ni]), [3 2 1]);
ax_norm=ax_norm(21:131,:,:); %reshape martrix from n_begin(21)-n_end(131) in x direction

ay_norm = permute( reshape( table2array( paren( ...
    readtable(filename + 'outay.dat'), ':', 4) ), ...
    [SIMION_nj, SIMION_nk, SIMION_ni]), [3 2 1]);
ay_norm=ay_norm(21:131,:,:);

az_norm = permute( reshape( table2array( paren( ...
    readtable(filename + 'outaz.dat'), ':', 4) ), ...
    [SIMION_nj, SIMION_nk, SIMION_ni]), [3 2 1]);
az_norm=az_norm(21:131,:,:);


% symmetrize the y, z fields NOT DONE HERE!!
%  ay_norm = (ay_norm - flip(ay_norm, 2))/2;
%  az_norm = (az_norm - flip(az_norm, 3))/2;

% save the variables ax_norm, ay_norm, az_norm in new mat files
new_filename = './acc/dec_norm_' + strrep( string(voltages(i)), '.', 'p') + 'kV/'
save(new_filename + 'a_norm', 'ax_norm', 'ay_norm', 'az_norm')
clearvars ax_norm ay_norm az_norm
end



%% Focusing mode field
fprintf('Loading focusing mode accelerations...')
voltages = [10, 12.5, 13.5];
for i=1:length(voltages)
    fieldFolder = '../dec_FM_' + strrep( string(voltages(i)), '.', 'p') + 'kV/output/';

    
    
    % Here again ax, ay, az were reshapded to n-begin-n-end
    % (21-131)
    ax_norm = permute( reshape( table2array( paren( ...
    readtable(fieldFolder + 'outax.dat'), ':', 4) ), ...
    [SIMION_nj, SIMION_nk, SIMION_ni]), [3 2 1]);
    ax_norm=ax_norm(21:131,:,:);
    
    
    ax_pos = permute( reshape( table2array( paren( ...
    readtable(fieldFolder + 'outax_fm_pos.dat'), ':', 4) ), ...
    [SIMION_nj, SIMION_nk, SIMION_ni]), [3 2 1]);
    ax_pos=ax_pos(21:131,:,:);
    
    
    ax_neg = permute( reshape( table2array( paren( ...
    readtable(fieldFolder + 'outax_fm_neg.dat'), ':', 4) ), ...
    [SIMION_nj, SIMION_nk, SIMION_ni]), [3 2 1]);
    ax_neg=ax_neg(21:131,:,:);
    
    
    ay_norm = permute( reshape( table2array( paren( ...
    readtable(fieldFolder + 'outay.dat'), ':', 4) ), ...
    [SIMION_nj, SIMION_nk, SIMION_ni]), [3 2 1]);
    ay_norm=ay_norm(21:131,:,:);
    
    ay_pos = permute( reshape( table2array( paren( ...
    readtable(fieldFolder + 'outay_fm_pos.dat'), ':', 4) ), ...
    [SIMION_nj, SIMION_nk, SIMION_ni]), [3 2 1]);
    ay_pos=ay_pos(21:131,:,:);
    
    ay_neg = permute( reshape( table2array( paren( ...
    readtable(fieldFolder + 'outay_fm_neg.dat'), ':', 4) ), ...
    [SIMION_nj, SIMION_nk, SIMION_ni]), [3 2 1]);
    ay_neg=ay_neg(21:131,:,:);
    
    
    az_norm = permute( reshape( table2array( paren( ...
    readtable(fieldFolder + 'outaz.dat'), ':', 4) ), ...
    [SIMION_nj, SIMION_nk, SIMION_ni]), [3 2 1]);
    az_norm=az_norm(21:131,:,:);
    
    
    az_pos = permute( reshape( table2array( paren( ...
    readtable(fieldFolder + 'outaz_fm_pos.dat'), ':', 4) ), ...
    [SIMION_nj, SIMION_nk, SIMION_ni]), [3 2 1]);
    az_pos=az_pos(21:131,:,:);
    
    az_neg = permute( reshape( table2array( paren( ...
    readtable(fieldFolder + 'outaz_fm_neg.dat'), ':', 4) ), ...
    [SIMION_nj, SIMION_nk, SIMION_ni]), [3 2 1]);
    az_neg=az_neg(21:131,:,:);
    
    % % symmetrize the y, z fields NOT DONE IN HERE
    % ay_norm = (ay_norm - flip(ay_norm, 2))/2;
    % az_norm = (az_norm - flip(az_norm, 3))/2;
    % 
    % % %            
    % % %                  also for focusing mode
    % ay_pos = (ay_pos + flip(ay_neg, 3))/2;
    % ay_neg = flip( ay_pos, 3);
    % 
    % az_pos = (az_pos - flip(az_neg, 3))/2;
    % az_neg = - flip(az_pos, 3);
    new_filename = './acc/dec_foc_' + strrep( string(voltages(i)), '.', 'p') + 'kV/'
    save(new_filename + 'a_norm', 'ax_norm', 'ay_norm', 'az_norm')
    save(new_filename + 'a_pos', 'ax_pos', 'ay_pos', 'az_pos')
    save(new_filename + 'a_neg', 'ax_neg', 'ay_neg', 'az_neg')

    clearvars ax_norm ax_pos ax_neg ay_norm ay_neg ay_pos az_norm az_pos az_neg
end
