classdef InputParameters < handle 
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        params % class containing all relavant parameters

        ax_norm, ax_pos, ax_neg % acceleration matrices; normal/focusing/focusing
        ay_norm, ay_pos, ay_neg
        az_norm, az_pos, az_neg
        ax_norm_1d_interpl, ax_neg_1d_interpl % acceleration along beam axis, only for generating time sequence
        
        
        ax_norm_extended, ay_norm_extended, az_norm_extended
        ax_norm_interpl,ay_norm_interpl,az_norm_interpl
        ax_norm_H_extended, ay_norm_H_extended, az_norm_H_extended
        ax_norm_H_interpl, ay_norm_H_interpl, az_norm_H_interpl
        ax_neg_extended, ay_neg_extended, az_neg_extended
        ax_neg_interpl,ay_neg_interpl,az_neg_interpl
        ax_neg_H_extended, ay_neg_H_extended, az_neg_H_extended
        ax_neg_H_interpl, ay_neg_H_interpl, az_neg_H_interpl
        ax_pos_extended, ay_pos_extended, az_pos_extended
        ax_pos_interpl,ay_pos_interpl,az_pos_interpl
        ax_pos_H_extended, ay_pos_H_extended, az_pos_H_extended
        ax_pos_H_interpl, ay_pos_H_interpl, az_pos_H_interpl
        
        
        % FORTRAN trigger sequence imported from the T2Jumps files
        T2Jump_time_vec, T2Jump_trigger_pattern, T2Jump_stage_number 

        % MATLAB trigger sequence (M=matlab), to be compared with Fortran
        M_time_vec, M_trigger_pattern, M_stage_number 
        M_synch_position, M_synch_velocity, M_synch_time % position, velocity and coordinate of the synchrounous molecule
        M_sequence_path % holds the path to the Matlab sequence (either generated or loaded)
    
        verbose % plots and print everything, false by default
        fortran_seq_bool % if false, ignore Fortran sequence (do not load it etc etc)
        always_generate_M_seq % force re-generation of M_seq regardless of whether it already exists
    end

    methods
        %% constructor of the class
        function obj = InputParameters( voltage, focusing_mode_bool, vel_synch_mol, target_velocity, options )
            
            arguments
                voltage                     double
                focusing_mode_bool          logical
                vel_synch_mol               double
                target_velocity             double
                options.Phase               double = 10    % default phase     
                options.Verbose             logical = false % default behaviour of verbose
                options.FortranSeqBool      logical = true % default behaviour for Fortran sequences
                options.AlwaysGenerateMSeq  logical = false % if true, always force new generation of M sequence
            end
            if nargin < 4
                 fprintf('Usage: InputParamters( voltage, focusing_mode_bool, V_synch_mol, V_target, options)\n(Options are: \t Phase\tVerbose\tFortranSeqBool\tAlwaysGenerateMSeq)\n')
            end

            % assign input args to class variables
            % mandatory args
            obj.params.FLY_voltage_on_electrodes = voltage; % choose between 10, 12.5, 13. Will fail if fields do not exists!
            obj.params.FLY_focusing_mode_bool = focusing_mode_bool; % flag for focusing mode operation
            obj.params.CALC_vel_synch_mol = vel_synch_mol; % vx velocity synchronous molecule (m/s)
            obj.params.FLY_target_velocity = target_velocity; % target velocity (mhhh, overlap with the phase angle, but is needed to load appropriate velocity sequence)
            % optional args
            obj.params.CALC_phase_degrees = options.Phase;
            obj.verbose = options.Verbose;
            obj.fortran_seq_bool = options.FortranSeqBool;
            obj.always_generate_M_seq = options.AlwaysGenerateMSeq;


            %% MANUAL DEFINITION OF ALL THE INPUT PARAMETRS
            % EVERYTHING MUST BE IN SI UNITS! 
            
            % Simion array -- DO NOT CHANGE
            obj.params.SIMION_ni = 151; % ni number of grid points x (along beam axis
            obj.params.SIMION_nj = 41; % number of grid points y,z (perpendicular to beam axis)
            obj.params.SIMION_nk = 41; % number of grid points y,z (perpendicular to beam axis)
            obj.params.SIMION_nbegin = 21; % Used part of array along beam axis, rest is for fitting. Start at nbegin...
            obj.params.SIMION_nend = 131; % ... ends at nend
            obj.params.SIMION_grid_units_p_meter = 20000.0; %gu #gridunits/meter (2*6666.6667)

            % Parameters of the PHYSICAL DECELERATOR
            obj.params.PHYS_valve_to_skimmer = 253.0e-3; % LA Nozzle to skimmer (m)
            obj.params.PHYS_skimmer_radius = 1.5e-3; %	r_s radius skimmer (m)
            obj.params.PHYS_skimmer_to_dec = 66.6e-3; % LB Skimmer to Hexapole or Decelerator (m)
            obj.params.PHYS_exit_to_detection = 11.52e-3; %L4 last stage decelerator to detection (m) % set to 10.2 mm 
            obj.params.PHYS_number_of_electrodes = 123; % nt (Number of electrodes in dec.-1) size!! 
            obj.params.PHYS_number_of_switching_times = 124; %# switching times in decelerator
            obj.params.PHYS_distance_stages = (obj.params.SIMION_nend - obj.params.SIMION_nbegin)/obj.params.SIMION_grid_units_p_meter; % distance between two stages, 5.5 mm
            obj.params.PHYS_valve_to_dec = obj.params.PHYS_valve_to_skimmer + obj.params.PHYS_skimmer_to_dec;
            obj.params.PHYS_length_dec = obj.params.PHYS_distance_stages * obj.params.PHYS_number_of_electrodes; % total length of decelerator
            obj.params.PHYS_valve_to_detection = obj.params.PHYS_valve_to_skimmer + obj.params.PHYS_skimmer_to_dec + obj.params.PHYS_length_dec + obj.params.PHYS_exit_to_detection;
            obj.params.PHYS_seperation_pins = (obj.params.SIMION_nj - 1)/obj.params.SIMION_grid_units_p_meter;
            
            % Parameters used in theoretical calculations
            obj.params.CALC_phase_distance = obj.params.CALC_phase_degrees * obj.params.PHYS_distance_stages / 180;

            obj.params.CALC_OH_mass_amu = 17; %	mass Mass of molecule in AMU
            obj.params.CALC_dipole_moment = 1.668; % mu dipole moment (Debye)
            obj.params.CALC_lambda = 1.649e9; % Lambda-doublet splitting (Hz) (it was in GHz before)
            obj.params.CALC_Beff = 0.58; % B effective value of MK/J(J+1)
            
            % Parameters used for the molecular beam
            obj.params.BEAM_avg_velocity_beam = 450; % Average velocity incoming package (m/s) probably wrong
            obj.params.BEAM_radius_of_nozzle = 0.25e-3; % radius of the valve nozzle (m)
            obj.params.BEAM_long_pos_spread = 11.5e-3; % longitudinal position spread (m) - along x aixs or beam propagation
            obj.params.BEAM_long_vel_spread = 0.20; % relative velocity spread along beam axis 0.112 0.12
            obj.params.BEAM_trans_velocity_spread = 0.1; % 	velocity spread perp. to beam axis/average velocity 0.10

            % Parameters used in fly and in the simulation
            obj.params.FLY_incoupling_time = 710.2e-6; % valve - decelerator incoupling time (s)
            obj.params.FLY_detection_laser_diameter = 1e-3;
            obj.params.FLY_simulated_target_vel = []; % to be assigned while loading or generating the Matlab sequence.
            % Stores the rounded value of the last velocity vector. Ideally
            % equals obj.params.FLY_target_velocity
            
            % add here any other parameter of interest. Can be an array
            % too. 

            % load the acceleration fields
            obj.loadAccelerationFields();
            
            % load Fortran Time sequence. 
            if obj.fortran_seq_bool
                obj.loadFortranTimeSequence();
            end 



            if obj.checkIfMatlabSequenceAlreadyExists && ~obj.always_generate_M_seq % if the Matlab sequence already exists, just load it
                fprintf('Matlab sequence already exists, I will just load it.\n')
                obj.loadMatlabSequence(); % it loads all the M_something variables
            else 
                fprintf('Matlab sequence not found. Generating a new one.\n')
                obj.generateMatlabTimeSequence();
            end
            
            obj.InterpolateAccField()

        end

        function setPhase(obj, new_phase) % mandatory to update the phase distance too!!
            obj.params.CALC_phase_degrees = new_phase;
            obj.params.CALC_phase_distance = obj.params.CALC_phase_degrees * obj.params.PHYS_distance_stages / 180;
        end
        
        %% Load the acceleration fields
        % (TODO: maybe the negative one can be obtained via a flip of the postive
        % one)
        function loadAccelerationFields(obj)
            if obj.params.FLY_focusing_mode_bool == false % load normal fields
                fprintf('Loading normal mode fields ...\t')
                fieldFolder = '../dec_Norm_' + strrep( string(obj.params.FLY_voltage_on_electrodes), '.', 'p') + 'kV/output/';
                paren = @(x, varargin) x(varargin{:}); % in-line function to reshape a matrix without a temporary variable
                obj.ax_norm = permute( reshape( table2array( paren( ...
                    readtable(fieldFolder + 'outax.dat'), ':', 4) ), ...
                    [obj.params.SIMION_nj, obj.params.SIMION_nk, obj.params.SIMION_ni]), [3 2 1]);
                obj.ax_norm=obj.ax_norm(21:131,:,:); %reshape martrix from n_begin(21)-n_end(131) in x direction

                obj.ay_norm = permute( reshape( table2array( paren( ...
                    readtable(fieldFolder + 'outay.dat'), ':', 4) ), ...
                    [obj.params.SIMION_nj, obj.params.SIMION_nk, obj.params.SIMION_ni]), [3 2 1]);
                obj.ay_norm=obj.ay_norm(21:131,:,:);

                obj.az_norm = permute( reshape( table2array( paren( ...
                    readtable(fieldFolder + 'outaz.dat'), ':', 4) ), ...
                    [obj.params.SIMION_nj, obj.params.SIMION_nk, obj.params.SIMION_ni]), [3 2 1]);
                obj.az_norm=obj.az_norm(21:131,:,:);
                % set the focusing mode fields to empty. Necessary if you
                % first load focusing mode fields and after you re-load
                % normal mode fields.
                if not(isempty(obj.ax_pos))
                    [obj.ax_pos, obj.ax_neg, obj.ay_pos, obj.ay_neg, obj.az_pos, obj.az_neg] = deal([], [], [], [], [], []);
                end
                
                % symmetrize the y, z fields
                 obj.ay_norm = (obj.ay_norm - flip(obj.ay_norm, 2))/2;
                 obj.az_norm = (obj.az_norm - flip(obj.az_norm, 3))/2;
                
            else % load focusing mode fields
                fprintf('Loading focusing mode accelerations...')
                fieldFolder = '../dec_FM_' + strrep( string(obj.params.FLY_voltage_on_electrodes), '.', 'p') + 'kV/output/';
%                fieldFolder = strrep(fieldFolder, '.', 'p'); % replace dot with p
                
                % follows ugly but fast code to read all the acc files. In
                % steps:
                % out = readtable('./inputs/outax_fm_neg.dat');
                % fields = table2array(out(:, 4));
                % acczyx = reshape(fields, [41, 41, 151]);
                % accxyz = permute(acczyx, [3 2 1]);
                paren = @(x, varargin) x(varargin{:}); % in-line function to reshape a matrix without a temporary variable
                % Here again ax, ay, az were reshapded to n-begin-n-end
                % (21-131)
                obj.ax_norm = permute( reshape( table2array( paren( ...
                    readtable(fieldFolder + 'outax.dat'), ':', 4) ), ...
                    [obj.params.SIMION_nj, obj.params.SIMION_nk, obj.params.SIMION_ni]), [3 2 1]);
                obj.ax_norm=obj.ax_norm(21:131,:,:);


                obj.ax_pos = permute( reshape( table2array( paren( ...
                    readtable(fieldFolder + 'outax_fm_pos.dat'), ':', 4) ), ...
                    [obj.params.SIMION_nj, obj.params.SIMION_nk, obj.params.SIMION_ni]), [3 2 1]);
                obj.ax_pos=obj.ax_pos(21:131,:,:);


                obj.ax_neg = permute( reshape( table2array( paren( ...
                    readtable(fieldFolder + 'outax_fm_neg.dat'), ':', 4) ), ...
                    [obj.params.SIMION_nj, obj.params.SIMION_nk, obj.params.SIMION_ni]), [3 2 1]);
                obj.ax_neg=obj.ax_neg(21:131,:,:);


                obj.ay_norm = permute( reshape( table2array( paren( ...
                    readtable(fieldFolder + 'outay.dat'), ':', 4) ), ...
                    [obj.params.SIMION_nj, obj.params.SIMION_nk, obj.params.SIMION_ni]), [3 2 1]);
                obj.ay_norm=obj.ay_norm(21:131,:,:);

                obj.ay_pos = permute( reshape( table2array( paren( ...
                    readtable(fieldFolder + 'outay_fm_pos.dat'), ':', 4) ), ...
                    [obj.params.SIMION_nj, obj.params.SIMION_nk, obj.params.SIMION_ni]), [3 2 1]);
                obj.ay_pos=obj.ay_pos(21:131,:,:);

                obj.ay_neg = permute( reshape( table2array( paren( ...
                    readtable(fieldFolder + 'outay_fm_neg.dat'), ':', 4) ), ...
                    [obj.params.SIMION_nj, obj.params.SIMION_nk, obj.params.SIMION_ni]), [3 2 1]);
                obj.ay_neg=obj.ay_neg(21:131,:,:);


                obj.az_norm = permute( reshape( table2array( paren( ...
                    readtable(fieldFolder + 'outaz.dat'), ':', 4) ), ...
                    [obj.params.SIMION_nj, obj.params.SIMION_nk, obj.params.SIMION_ni]), [3 2 1]);
                obj.az_norm=obj.az_norm(21:131,:,:);

               
                obj.az_pos = permute( reshape( table2array( paren( ...
                    readtable(fieldFolder + 'outaz_fm_pos.dat'), ':', 4) ), ...
                    [obj.params.SIMION_nj, obj.params.SIMION_nk, obj.params.SIMION_ni]), [3 2 1]);
                obj.az_pos=obj.az_pos(21:131,:,:);

                obj.az_neg = permute( reshape( table2array( paren( ...
                    readtable(fieldFolder + 'outaz_fm_neg.dat'), ':', 4) ), ...
                    [obj.params.SIMION_nj, obj.params.SIMION_nk, obj.params.SIMION_ni]), [3 2 1]);
                obj.az_neg=obj.az_neg(21:131,:,:);
                clearvars fieldFolder
                
                % symmetrize the y, z fields
                 obj.ay_norm = (obj.ay_norm - flip(obj.ay_norm, 2))/2;
                 obj.az_norm = (obj.az_norm - flip(obj.az_norm, 3))/2;
                
                
            end
            fprintf('\tloaded\n')
        end
        
        %% Interpolate acceleration field
        function InterpolateAccField(obj)
            num_grids_x = 111 + 110 * 122 + 20;
            num_grids_y = 41 + 20;
            num_grids_z = 41 + 20;
            gridded_x = linspace(-10/obj.params.SIMION_grid_units_p_meter, obj.params.PHYS_length_dec + 10/obj.params.SIMION_grid_units_p_meter, num_grids_x);
            gridded_y = linspace(-10/obj.params.SIMION_grid_units_p_meter-obj.params.PHYS_seperation_pins/2.0, obj.params.PHYS_seperation_pins/2.0 + 10/obj.params.SIMION_grid_units_p_meter, num_grids_y);
            gridded_z =	linspace(-10/obj.params.SIMION_grid_units_p_meter-obj.params.PHYS_seperation_pins/2.0, obj.params.PHYS_seperation_pins/2.0 + 10/obj.params.SIMION_grid_units_p_meter, num_grids_z);

            obj.ax_norm_extended = zeros(num_grids_x, num_grids_y, num_grids_z); %Vertical ones
            obj.ay_norm_extended = zeros(num_grids_x, num_grids_y, num_grids_z);
            obj.az_norm_extended = zeros(num_grids_x, num_grids_y, num_grids_z);
            % below changed ax_norm(21:131,:,:) to ax_norm since already cut 
            % not sure about flip argument 22:130 --> 1:109 in cut matrix
            % correct?
            
            obj.ax_norm_extended(11:num_grids_x-10,11:51,11:51) = cat(1, repmat(cat(1,obj.ax_norm, flip(-obj.ax_norm(1:109,:,:),1)), 61, 1, 1), obj.ax_norm);  % before  obj.ax_norm_extended(11:num_grids_x-10,11:51,11:51) = cat(1, repmat(cat(1,obj.ax_norm, flip(-obj.ax_norm(22:130,:,:),1)), 61, 1, 1), obj.ax_norm);
            obj.ay_norm_extended(11:num_grids_x-10,11:51,11:51) = cat(1, repmat(cat(1,obj.ay_norm, flip(obj.ay_norm(1:109,:,:),1)), 61, 1, 1), obj.ay_norm);
            obj.az_norm_extended(11:num_grids_x-10,11:51,11:51) = cat(1, repmat(cat(1,obj.az_norm, flip(obj.az_norm(1:109,:,:),1)), 61, 1, 1), obj.az_norm);
            obj.ax_norm_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.ax_norm_extended,'linear','linear');
            obj.ay_norm_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.ay_norm_extended,'linear','linear');
            obj.az_norm_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.az_norm_extended,'linear','linear');

            obj.ax_norm_H_extended = zeros(num_grids_x, num_grids_y, num_grids_z); %Horizontal
            obj.ay_norm_H_extended = zeros(num_grids_x, num_grids_y, num_grids_z);
            obj.az_norm_H_extended = zeros(num_grids_x, num_grids_y, num_grids_z);
            obj.ax_norm_H_extended(11:num_grids_x-10,11:51,11:51) = cat(1, repmat(cat(1,flip(-obj.ax_norm,1), obj.ax_norm(1:109,:,:)), 61, 1, 1), flip(-obj.ax_norm,1));
            obj.ay_norm_H_extended(11:num_grids_x-10,11:51,11:51) = cat(1, repmat(cat(1,flip(obj.ay_norm,1), obj.ay_norm(1:109,:,:)), 61, 1, 1), flip(obj.ay_norm,1));
            obj.az_norm_H_extended(11:num_grids_x-10,11:51,11:51) = cat(1, repmat(cat(1,flip(obj.az_norm,1), obj.az_norm(1:109,:,:)), 61, 1, 1), flip(obj.az_norm,1));
            obj.ax_norm_H_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.ax_norm_H_extended,'linear','linear');
            obj.ay_norm_H_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.ay_norm_H_extended,'linear','linear');% exchange z and y
            obj.az_norm_H_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.az_norm_H_extended,'linear','linear');
            
            if obj.params.FLY_focusing_mode_bool
                obj.ax_neg_extended = zeros(num_grids_x, num_grids_y, num_grids_z);
                obj.ay_neg_extended = zeros(num_grids_x, num_grids_y, num_grids_z);
                obj.az_neg_extended = zeros(num_grids_x, num_grids_y, num_grids_z);
                obj.ax_neg_extended(11:num_grids_x-10,11:51,11:51) = cat(1, repmat(cat(1,obj.ax_neg, flip(-obj.ax_neg(1:109,:,:),1)), 61, 1, 1), obj.ax_neg);
                obj.ay_neg_extended(11:num_grids_x-10,11:51,11:51) = cat(1, repmat(cat(1,obj.ay_neg, flip(obj.ay_neg(1:109,:,:),1)), 61, 1, 1), obj.ay_neg);
                obj.az_neg_extended(11:num_grids_x-10,11:51,11:51) = cat(1, repmat(cat(1,obj.az_neg, flip(obj.az_neg(1:109,:,:),1)), 61, 1, 1), obj.az_neg);
                obj.ax_neg_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.ax_neg_extended,'linear','linear');
                obj.ay_neg_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.ay_neg_extended,'linear','linear');
                obj.az_neg_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.az_neg_extended,'linear','linear');

                obj.ax_neg_H_extended = zeros(num_grids_x, num_grids_y, num_grids_z);
                obj.ay_neg_H_extended = zeros(num_grids_x, num_grids_y, num_grids_z);
                obj.az_neg_H_extended = zeros(num_grids_x, num_grids_y, num_grids_z);
                obj.ax_neg_H_extended(11:num_grids_x-10,11:51,11:51) = cat(1, repmat(cat(1,flip(-obj.ax_neg,1), obj.ax_neg(1:109,:,:)), 61, 1, 1), flip(-obj.ax_neg,1));
                obj.ay_neg_H_extended(11:num_grids_x-10,11:51,11:51) = cat(1, repmat(cat(1,flip(obj.ay_neg,1), obj.ay_neg(1:109,:,:)), 61, 1, 1), flip(obj.ay_neg,1));
                obj.az_neg_H_extended(11:num_grids_x-10,11:51,11:51) = cat(1, repmat(cat(1,flip(obj.az_neg,1), obj.az_neg(1:109,:,:)), 61, 1, 1), flip(obj.az_neg,1));
                obj.ax_neg_H_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.ax_neg_H_extended,'linear','linear');
                obj.ay_neg_H_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.ay_neg_H_extended,'linear','linear');
                obj.az_neg_H_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.az_neg_H_extended,'linear','linear');
                
                obj.ax_pos_extended = zeros(num_grids_x, num_grids_y, num_grids_z);
                obj.ay_pos_extended = zeros(num_grids_x, num_grids_y, num_grids_z);
                obj.az_pos_extended = zeros(num_grids_x, num_grids_y, num_grids_z);
                obj.ax_pos_extended(11:num_grids_x-10,11:51,11:51) = cat(1, repmat(cat(1,obj.ax_pos, flip(-obj.ax_pos(1:109,:,:),1)), 61, 1, 1), obj.ax_pos);
                obj.ay_pos_extended(11:num_grids_x-10,11:51,11:51) = cat(1, repmat(cat(1,obj.ay_pos, flip(obj.ay_pos(1:109,:,:),1)), 61, 1, 1), obj.ay_pos);
                obj.az_pos_extended(11:num_grids_x-10,11:51,11:51) = cat(1, repmat(cat(1,obj.az_pos, flip(obj.az_pos(1:109,:,:),1)), 61, 1, 1), obj.az_pos);
                obj.ax_pos_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.ax_pos_extended,'linear','linear');
                obj.ay_pos_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.ay_pos_extended,'linear','linear');
                obj.az_pos_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.az_pos_extended,'linear','linear');

                obj.ax_pos_H_extended = zeros(num_grids_x, num_grids_y, num_grids_z);
                obj.ay_pos_H_extended = zeros(num_grids_x, num_grids_y, num_grids_z);
                obj.az_pos_H_extended = zeros(num_grids_x, num_grids_y, num_grids_z);
                obj.ax_pos_H_extended(11:num_grids_x-10,11:51,11:51) = cat(1, repmat(cat(1,flip(-obj.ax_pos,1), obj.ax_pos(1:109,:,:)), 61, 1, 1), flip(-obj.ax_pos,1));
                obj.ay_pos_H_extended(11:num_grids_x-10,11:51,11:51) = cat(1, repmat(cat(1,flip(obj.ay_pos,1), obj.ay_pos(1:109,:,:)), 61, 1, 1), flip(obj.ay_pos,1));
                obj.az_pos_H_extended(11:num_grids_x-10,11:51,11:51) = cat(1, repmat(cat(1,flip(obj.az_pos,1), obj.az_pos(1:109,:,:)), 61, 1, 1), flip(obj.az_pos,1));
                obj.ax_pos_H_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.ax_pos_H_extended,'linear','linear');
                obj.ay_pos_H_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.ay_pos_H_extended,'linear','linear');
                obj.az_pos_H_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.az_pos_H_extended,'linear','linear');
            end
        end
        

        %% loadFortranTimeSequence
        % from appropriate folder and save it to the class variables
        % MANDATORY TO FIRST RUN THE SH CODE THAT CLEANS UP THE FILES...
        % matlab goes nut for the variable number of whitespaces. Python
        % should be fine with it. 
        function loadFortranTimeSequence(obj)
            if obj.params.FLY_focusing_mode_bool == false % folder path for normal mode
                timeSequenceFolder = ...
                    '../dec_Norm_' + strrep( string(obj.params.FLY_voltage_on_electrodes), '.', 'p') ...
                    + 'kV/timeseq/sequences/dec_' + num2str(obj.params.CALC_vel_synch_mol) + ...
                    '_' + num2str(obj.params.FLY_target_velocity) + '/T2jump_c.out' ;
            else % folder path for focusing files
                timeSequenceFolder = ...
                    '../dec_FM_' + strrep( string(obj.params.FLY_voltage_on_electrodes), '.', 'p') + ...
                    'kV/timeseq/sequences/dec_' + num2str(obj.params.CALC_vel_synch_mol) + ...
                    '_' + num2str(obj.params.FLY_target_velocity) + '/T2jump_c.out';
            end
            opt = detectImportOptions(timeSequenceFolder, 'FileType', 'text');
            opt.CommentStyle = '#'; opt.Delimiter = {'\t' ' '}; 
            opt.VariableTypes = {'double', 'char', 'char', 'uint8'};
            out = readtable(timeSequenceFolder, opt);
            obj.T2Jump_time_vec = table2array(out(:, 1));
            obj.T2Jump_trigger_pattern = cell2mat(table2array(out(:, 2)));
            obj.T2Jump_trigger_pattern = obj.T2Jump_trigger_pattern(:, (end-3):end); % gets rid of the b or 0x in front
            obj.T2Jump_stage_number = uint8(table2array(out(:, 4))); % also force uint8
            fprintf('Fortran time sequence loaded:\t\t' + timeSequenceFolder + '\n')
            clearvars out timeSequenceFolder opt;
        end
        
        %% loadFortranTimeSequence from T2jump.out
        % from appropriate folder and save it to the class variables
        % MANDATORY TO FIRST RUN THE SH CODE THAT CLEANS UP THE FILES...
        % matlab goes nut for the variable number of whitespaces. Python
        % should be fine with it. 
        function loadFortranTimeSequence2(obj)
            if obj.params.FLY_focusing_mode_bool == false % folder path for normal mode
                timeSequenceFolder = ...
                    '../dec_Norm_' + strrep( string(obj.params.FLY_voltage_on_electrodes), '.', 'p') ...
                    + 'kV/timeseq/sequences/dec_' + num2str(obj.params.CALC_vel_synch_mol) + ...
                    '_' + num2str(obj.params.FLY_target_velocity) + '/T2jump.out' ;
            else % folder path for focusing files
                timeSequenceFolder = ...
                    '../dec_FM_' + strrep( string(obj.params.FLY_voltage_on_electrodes), '.', 'p') + ...
                    'kV/timeseq/sequences/dec_' + num2str(obj.params.CALC_vel_synch_mol) + ...
                    '_' + num2str(obj.params.FLY_target_velocity) + '/T2jump.out';
            end
            fid = fopen(timeSequenceFolder,'r');
            out = textscan(fid,'%d %s %s %d','headerlines', 4); fclose(fid);
            obj.T2Jump_time_vec = out{1,1};
            obj.T2Jump_trigger_pattern = string(out{1,2});
%             obj.T2Jump_trigger_pattern = obj.T2Jump_trigger_pattern(:, (end-3):end); % gets rid of the b or 0x in front
            obj.T2Jump_stage_number = uint8(out{1,4}); % also force uint8
            fprintf('Fortran time sequence loaded:\t\t' + timeSequenceFolder + '\n')
            clearvars out timeSequenceFolder opt;
        end


        %% generateMatlabTimeSequence
        % This function integrates numerically the sequence, both in normal and focusing mode.
        % Returns - actually reassign - the variables:
        % M_time_vec, M_trigger_pattern, M_stage_number -- containing the
        % time vector of the sequence and the trigger pattern
        % and the variables
        % M_synch_position, M_synch_velocity, M_synch_time
        % that keeps the position, velocity and time on the x
        % (decelerator) coordinate of the synchronous molecules, which can
        % be quite usefull later on in the full simulation

        function [simulated_target_vel] = generateMatlabTimeSequence(obj, verbose)
            % verbose is optional argument to plot the result
            if nargin == 1 
                verbose = obj.verbose; % set to default
            end
            fprintf('Generating Matlab time sequence...');

            % 111 pops out from?
            %changed ax norm again since not cut anymore
            obj.ax_norm_1d_interpl = griddedInterpolant( linspace(0, obj.params.PHYS_distance_stages, 111), obj.ax_norm(:,21,21),'linear','none');          
            if obj.params.FLY_focusing_mode_bool
                obj.ax_neg_1d_interpl = griddedInterpolant(linspace(0,obj.params.PHYS_distance_stages,111), obj.ax_neg(:,21,21),'linear','none');
            end
            
            
            function [value, isterminal, direction] = EventsFcn(t, x) % This event function stops the ode solver once the molecule arrives at detection point
                value = x(1) < obj.params.PHYS_length_dec + obj.params.PHYS_exit_to_detection;
                isterminal = 1; 
                direction = 0;
            end

           % start integration
            use_ode_solver_bool = true;
            tic;
            if use_ode_solver_bool
                opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-8, 'Events', @(t,x) EventsFcn(t,x));
                [obj.M_synch_time, x_Vx_temp] = ode45( @(t,x) ...
                    obj.dxdt(t,x), [0, 5e-3], [0; obj.params.CALC_vel_synch_mol], opts); % ode23t seems a good solver
            else
                xx = [0; obj.params.CALC_vel_synch_mol];
                t_step = 10e-9;
                tt= 0:t_step:5e-3;
                x_Vx = zeros(length(tt),2);
                for i = 1:1:length(tt)
                    if xx(1) < obj.params.PHYS_length_dec + obj.params.PHYS_exit_to_detection
                        x_Vx(i,:) = xx;
                        xx = xx + obj.dxdt(tt(i), xx)*t_step;
                    else
                        break;
                    end
                end
                obj.M_synch_time = tt(1:i-1)';
                x_Vx_temp = x_Vx(1:i-1,:);
            end
            toc;
            % end of the numerical integration

            obj.M_synch_position = x_Vx_temp(:, 1); % reshape x and Vx
            obj.M_synch_velocity = x_Vx_temp(:, 2); clearvars x_Vx_temp; % ugly

            t_x_interpl = griddedInterpolant(obj.M_synch_position, obj.M_synch_time); % to obatain the exact field switching time based on the switching position
            if obj.params.FLY_focusing_mode_bool
                positions_to_jump = [0, union(obj.params.PHYS_distance_stages/2 + obj.params.CALC_phase_distance: obj.params.PHYS_distance_stages:...
                    obj.params.PHYS_length_dec-(obj.params.PHYS_distance_stages/2 - obj.params.CALC_phase_distance), obj.params.PHYS_distance_stages*3.0/2.0 -...
                    obj.params.CALC_phase_distance: obj.params.PHYS_distance_stages:obj.params.PHYS_length_dec - (obj.params.PHYS_distance_stages/2 + obj.params.CALC_phase_distance)), ...
                    obj.params.PHYS_length_dec + obj.params.PHYS_exit_to_detection]'; % list of positions of interest
                obj.M_time_vec = t_x_interpl(positions_to_jump);
                trigger_pattern = repmat(["b0011";"b0100";"b1100";"b0010";"b0011";"b1000";"b1100";"b0001"], obj.params.PHYS_number_of_electrodes,1);
                obj.M_trigger_pattern = trigger_pattern(1:length(obj.M_time_vec));
                obj.M_trigger_pattern(end-1:end)= "b0000";
            else
                positions_to_jump = [0, obj.params.PHYS_distance_stages/2 + obj.params.CALC_phase_distance: obj.params.PHYS_distance_stages:...
                    obj.params.PHYS_length_dec-(obj.params.PHYS_distance_stages/2 - ...
                    obj.params.CALC_phase_distance), obj.params.PHYS_length_dec + obj.params.PHYS_exit_to_detection]'; % list of positions of interest
                obj.M_time_vec = t_x_interpl(positions_to_jump);
                trigger_pattern = repmat(["0x0010";"0x0020"], obj.params.PHYS_number_of_electrodes,1);
                obj.M_trigger_pattern = trigger_pattern(1:length(obj.M_time_vec));
                obj.M_trigger_pattern(end-1:end)= "0x0000";
            end
            
            obj.params.FLY_simulated_target_vel = round( obj.M_synch_velocity(end), 0);
            % stores in a quite useless but safe variable the effective final velocity
            % Will be used to create the filename

            obj.saveMatlabSequence(); % saves both .out and .m files for the sequence just created
            % in the .out file we save M_time_vec and M_trigger_pattern
            % in the .mat file we save all the files with M_'something'

            %             if obj.params.CALC_save_sequence_bool
%                 seq_file = fopen(obj.M_sequence_path,'w'); % will fail if folder does not exists
%                 fprintf(seq_file, '%s\t%s\n', [string(round(obj.M_time_vec*1e9)), obj.M_trigger_pattern]');
%                 fclose(seq_file);
%             end

            if verbose % we plot the synch molecule
                figure('Name', 'Time seq with Matlab');

                subplot(2, 3, 1)    % x vs t
                plot(obj.M_synch_time .* 1e3, obj.M_synch_position);
                title('Position vs time'); ylabel('x (m)'); xlabel('t (ms)')

                % linear trendline, see below
                linear_velocity = obj.params.CALC_vel_synch_mol - ...
                        (obj.params.CALC_vel_synch_mol - ...21:131
                        obj.M_synch_velocity(end)) / ...
                        obj.M_synch_time(end) .* obj.M_synch_time;
                % THE FINAL TIME IS OVERESTIMATED, BECAUSE THAT IS THE TIME
                % AT WHICH IT EXITs
                % THE DECELERATOR; BUT THE DECELRATION
                % STOPS EARLIER; TODO TO BE FIXED

                subplot(2, 3, 2) % Vx vs t
                plot(obj.M_synch_time .* 1e3, obj.M_synch_velocity );            
                title('Velocity vs time'); ylabel('Vx (m/s)'); xlabel('t (ms)')

                subplot(2, 3, 4) % Vx vs x
                plot(obj.M_synch_position, obj.M_synch_velocity );            
                title('Velocity vs space'); ylabel('Vx (m/s)'); xlabel('x (m)')

                subplot(2, 3, 5) % Vx vs t, minus a trendline of linear deceleration
                % that follows V(t) =  V_0 - (V_0 - V_final) / t_final * t
                % whre V_0 = starting velocity, V_final t_final velocity
                % and time of the last instant of the simulation
                plot(obj.M_synch_time .* 1e3, obj.M_synch_velocity - linear_velocity);
                title('Vx (m/s) - linear decrease'); ylabel('V (m/s)'); xlabel('t (ms)');
                
                % plot the single timesteps, to see how they very in the
                % variable-timestep ODE solvers of MATLAB. 
                % get rid of first and few last ones as they screw up the
                % plotting
                subplot(2, 3, 3)
                time_step_difference = circshift(obj.M_synch_time*1e9, 1) - obj.M_synch_time*1e9;
                time_step_difference = - circshift(obj.M_synch_time, 1) + obj.M_synch_time;
                time_step_difference = time_step_difference(2:end-7);
                plot(time_step_difference * 1e6, '--o')
                xlabel('Index of vector'); ylabel('Intergation timestep (ns)'); title('Time steps of the integration')

                subplot(2, 3, 6)
                histogram(time_step_difference, 100)
                xlabel('Intergation timestep (ns)'); title('Histrogram of time steps')

            end
            fprintf('\tFinal Matlab velocity is %d\n', obj.params.FLY_simulated_target_vel)
            simulated_target_vel = obj.params.FLY_simulated_target_vel; 
            % I return the compute velocity, to be used if needed.
            % for some reason I cannot return the class variable but I must
            % declare a local one. Probably due to handle and class stuff.
        end
        
        function dxdt = dxdt(obj, t, x)
            if obj.params.FLY_focusing_mode_bool 
                if x(1) > obj.params.PHYS_length_dec - (obj.params.PHYS_distance_stages/2 - obj.params.CALC_phase_distance)
                    dxdt = [x(2); 0];
                elseif x(1) < obj.params.PHYS_distance_stages/2 + obj.params.CALC_phase_distance
                    dxdt = [x(2); obj.ax_norm_1d_interpl(x(1))];
                else
                    pos = mod(x(1), obj.params.PHYS_distance_stages);
                    if pos < obj.params.PHYS_distance_stages/2 - obj.params.CALC_phase_distance
                        dxdt = [x(2); obj.ax_neg_1d_interpl(pos)];
                    elseif pos > obj.params.PHYS_distance_stages/2 + obj.params.CALC_phase_distance
                        dxdt = [x(2); -obj.ax_neg_1d_interpl(obj.params.PHYS_distance_stages-pos)];
                    else
                        dxdt = [x(2); obj.ax_norm_1d_interpl(pos)];
                    end
                end
            else
                if x(1) > obj.params.PHYS_length_dec-(obj.params.PHYS_distance_stages/2 - obj.params.CALC_phase_distance)
                    dxdt = [x(2);0];
                else
                    pos = mod(x(1),obj.params.PHYS_distance_stages);
                    if pos <= obj.params.PHYS_distance_stages/2 + obj.params.CALC_phase_distance
                        dxdt = [x(2); obj.ax_norm_1d_interpl(pos)];
                    else
                        dxdt = [x(2); -obj.ax_norm_1d_interpl(obj.params.PHYS_distance_stages-pos)];
                    end
                end
            end
        end

        %% checkIfMatlabSequenceAlreadyExists
        % This function simply checks if the sequence currently loaded
        % already exists in the folder 
        % returns true (1) if it already exits, otherwise returns false (0).
        % Neat as forces re-loading of filepath, likely overkilled.
        % CHECK EXISTANCE OF BOTH .OUT AND .MAT FILES
        function [sequence_exist] = checkIfMatlabSequenceAlreadyExists(obj)
            obj.makeMatlabSequencePath(); % re-make the path for safety
            % make the .mat file extension to check existance of both
            [filepath, name, ~] = fileparts(obj.M_sequence_path); ext='.mat'; % remove .out
            filename_mat = fullfile(filepath, strcat(name,ext)); % glue toghether
            sequence_exist = isfile(obj.M_sequence_path) & isfile(filename_mat); % both .out and .mat files must exits
        end

        % This function checks if the Fortran sequence exits
        function [sequence_exist] = checkIfFortranSequenceExists(obj)
            if obj.params.FLY_focusing_mode_bool == false % folder path for normal mode
                timeSequenceFolder = ...
                    '../dec_Norm_' + strrep( string(obj.params.FLY_voltage_on_electrodes), '.', 'p') ...
                    + 'kV/timeseq/sequences/dec_' + num2str(obj.params.CALC_vel_synch_mol) + ...
                    '_' + num2str(obj.params.FLY_target_velocity) + '/T2jump_c.out' ;
            else % folder path for focusing files
                timeSequenceFolder = ...
                    '../dec_FM_' + strrep( string(obj.params.FLY_voltage_on_electrodes), '.', 'p') + ...
                    'kV/timeseq/sequences/dec_' + num2str(obj.params.CALC_vel_synch_mol) + ...
                    '_' + num2str(obj.params.FLY_target_velocity) + '/T2jump_c.out';
            end
            sequence_exist = isfile(timeSequenceFolder);
            if ~sequence_exist
                fprintf('Fortran sequence not found at filename %s\n', timeSequenceFolder)
            end
        end

        %% makeMatlabSequencePath
        % this function writes/overwrites in obj.M_sequence_path the path
        % of the matlab sequence. As extra argument you can specify the
        % final velocity you just simulated in generateMatlabSequence
        % function, which is what you need to create the new path
        % It is based on the .out file, but variables of synch molecule
        % should be saved in a .mat or .txt file too
        function makeMatlabSequencePath(obj, simulated_final_velocity)
            if nargin == 1 % target velocity if not specified
                simulated_final_velocity = obj.params.FLY_target_velocity;
            else 
                simulated_final_velocity = round(simulated_final_velocity, 1); % DANGEROUS AF
            end
            if obj.params.FLY_focusing_mode_bool
                obj.M_sequence_path = './sequences/focusing' + ...
                    strrep( string(obj.params.FLY_voltage_on_electrodes), '.', 'p') + ...
                    'kV/' + 'dec_' + string(obj.params.CALC_vel_synch_mol) ...
                    + '_' + string(simulated_final_velocity) + '.out';
            else
                obj.M_sequence_path = './sequences/norm' + ...
                    strrep( string(obj.params.FLY_voltage_on_electrodes), '.', 'p') + ...
                    'kV/' + 'dec_' + string(obj.params.CALC_vel_synch_mol) ...
                    + '_' + string(simulated_final_velocity) + '.out';
            end
        end
        
        %% saveMatlabSequence
        % in the .out file we save M_time_vec and M_trigger_pattern
        % in the .mat file we save all the files with M_'something'
        function saveMatlabSequence(obj)
            if isempty(obj.params.FLY_simulated_target_vel) % for safety
                fprintf('NOOB You are trying to save a sequence that was never simulated. Whaaat?\nBadly wrong, return without saving it.\n')
                return
            end
            obj.makeMatlabSequencePath(obj.params.FLY_simulated_target_vel); % rename/recreate the path filename of the sequence
            
            % save the .out file for the experiment
            seq_file = fopen(obj.M_sequence_path,'w'); % will fail if folder does not exists
            fprintf(seq_file, '%s\t%s\n', [ string( round( obj.M_time_vec*1e9) ), obj.M_trigger_pattern]' );
            fclose(seq_file);
            
            % rename filename to .mat to save Matlab variable
            [filepath, name, ~] = fileparts(obj.M_sequence_path); ext='.mat'; % remove .out
            filename_mat = fullfile(filepath, strcat(name,ext)); % glue toghether

            % follows ugly piece of code that creates local variables out
            % of "handle" of class variables to save them into a file. 

            M_time_vec = obj.M_time_vec;
            M_trigger_pattern = obj.M_trigger_pattern;
            M_stage_number = obj.M_stage_number;
            M_synch_position = obj.M_synch_position;
            M_synch_velocity = obj.M_synch_velocity;
            M_synch_time = obj.M_synch_time;
            M_sequence_path = obj.M_sequence_path;
            fprintf('Saving Matlab sequence in folder \t%s\n', obj.M_sequence_path)
            save(filename_mat, 'M_time_vec', "M_trigger_pattern", ...
                "M_stage_number", "M_synch_position", "M_synch_velocity", ...
                "M_synch_time", "M_sequence_path", '-mat')
            clearvars M_time_vec M_trigger_pattern M_stage_number ...
                M_synch_position M_synch_velocity M_synch_time M_sequence_path % important to delete them
        end

        %% loadMatlabSequence
        % this function loads the matlab sequence into the class' variables
        % named M_something form the .mat file
        function loadMatlabSequence(obj)
            fprintf('Loading Matlab time sequence...')
            obj.makeMatlabSequencePath()
            % get .mat filename 
            [filepath, name, ~] = fileparts(obj.M_sequence_path); ext='.mat'; % remove .out
            filename_mat = fullfile(filepath, strcat(name,ext)); % glue toghether
            load( filename_mat );
            obj.M_time_vec = M_time_vec;
            obj.M_trigger_pattern = M_trigger_pattern;
            obj.M_stage_number = M_stage_number;
            obj.M_synch_position = M_synch_position;
            obj.M_synch_velocity = M_synch_velocity;
            obj.M_synch_time = M_synch_time;
            obj.M_sequence_path = M_sequence_path;
            obj.params.FLY_simulated_target_vel = round( obj.M_synch_velocity(end) );
            clearvars M_time_vec M_trigger_pattern M_stage_number ...
                M_synch_position M_synch_velocity M_synch_time M_sequence_path % important to delete them            
            fprintf('\t\tdone\n')
        end

        %% changeFieldConfig
        % This function re-loads the fields and must be used whenever
        % you want to change voltage values or focusing/normal mode
        function changeFieldConfig(obj, new_voltage, new_focusing_mode_bool)
            if nargin ~= 3 
                fprintf('Wrong changeFieldConfig call.\nUsage changeFieldConfig( new_voltage, new_focusing_mode_bool)\n')
                return
            else
            % overwrite parameters
                obj.params.FLY_voltage_on_electrodes = new_voltage;
                obj.params.FLY_focusing_mode_bool = new_focusing_mode_bool;
            end
            % re-load acceleration fields
            obj.loadAccelerationFields();

            % re-load Fortran sequence, if fortran_seq_bool
            if obj.fortran_seq_bool
                obj.loadFortranTimeSequence();
            end 


            % re-load Matlab sequence if exists, otherwise re-generate it
            if obj.checkIfMatlabSequenceAlreadyExists % if the Matlab sequence already exists, just load it
                fprintf('Matlab sequence already exists, I will just load it.\n')
                obj.loadMatlabSequence(); % it loads all the M_something variables
            else 
                fprintf('Matlab sequence not found. Generating a new one.\n')
                obj.generateMatlabTimeSequence();
            end
            fprintf('TODO: to be tested properly')
        end

        %% changeVelocities
        % This function must be used to change the Velocities of the 
        % synchronous molecules and the target velocity
        % optional argument is the new phase
        function changeVelocities(obj, new_synch_mol_velocity, new_target_velocity, new_phase)
            if nargin ~= 3 && nargin ~=4
                fprintf('Wrong changeVelocities call.\nUsage changeVelocities( new_synch_mol_velocity, new_target_velocity, new_phase)\n(Phase is optional)\n')
                return
            else
            % overwrite parameters
                obj.params.CALC_vel_synch_mol = new_synch_mol_velocity;
                obj.params.FLY_target_velocity = new_target_velocity;
                if nargin == 4 % give new phase too
                    obj.setPhase( new_phase );
                end
            end
            
            % re-load Fortran sequence, if fortran_seq_bool
            if obj.fortran_seq_bool
                if obj.checkIfFortranSequenceExists()
                    obj.loadFortranTimeSequence();
                else
                    fprintf("The Fortran sequence you want to load does not exits. Set the flag obj.fortran_seq_bool to false form now on.\n")
                    obj.fortran_seq_bool = false;
                end
            end 


            % re-load Matlab sequence if exists, otherwise re-generate it
            if obj.checkIfMatlabSequenceAlreadyExists && ~obj.always_generate_M_seq % if the Matlab sequence already exists, just load it
                fprintf('Matlab sequence already exists, I will just load it.\n')
                obj.loadMatlabSequence(); % it loads all the M_something variables
            else 
                fprintf('Matlab sequence not found. Generating a new one.\n')
                obj.generateMatlabTimeSequence();
            end
            fprintf("Velocities changed\n")
            fprintf('TODO: to be tested properly')
        end

        %% compareSequences
        % compare Fortran w Matlab sequences
        % (Fortran is in ns...)
        function compareSequences(obj)
            if obj.fortran_seq_bool == false
                fprintf('Cannot compare cause Fortran sequence is not loaded\n')
                return
            end
            fprintf("Final velocity simulated with Matlab: %i\nFinal " + ...
                "velocity simulated with Fortran %i\nPhase of Matlab simulation: %d\nPrecision of +-1 m/s\n", ...
            obj.params.FLY_simulated_target_vel, obj.params.FLY_target_velocity, ...
            obj.params.CALC_phase_degrees);

            if obj.params.FLY_simulated_target_vel ~= obj.params.FLY_target_velocity
            fprintf("*** MATLAB and Fortran velocities are off! ***\nChange the phase of Matlab till matched\n")
            end
            
            % plot everything
            figure("Name", "Comparison Fortran - Matlab")
            subplot(2, 1, 1)
            plot(obj.T2Jump_time_vec, '-o', 'DisplayName', 'T2Jump_time_vec'); 
            hold on; plot(obj.M_time_vec .*1e9 , '-o', 'DisplayName', 'M_time_vec');
            xlabel('vector index'); ylabel('time (ns)'); legend()
            subplot(2, 1, 2)
            if obj.params.FLY_focusing_mode_bool
                plot(obj.T2Jump_time_vec(1:end-2) - obj.M_time_vec(1:end-1) *1e9, '-o', 'DisplayName', 'difference'); 
                fprintf('Removed one point in Matlab time vector and 2 points in Fortran time vector\n')
            else
                % there seems to be an offset of 1010 us due to some fishy
                % Fortran leftovers
                plot(obj.T2Jump_time_vec - 1010 - (obj.M_time_vec).*1e9, '-o', 'DisplayName', 'difference'); 
            end
            xlabel('vector index'); ylabel('time (ns)'); legend()

        end

        %% plotAccelerationFields
        % plot the fields
        function plotAccelerationFields(obj)
            fprintf('Plotting acceleration fields ...\t')

            xrange = [obj.params.SIMION_nbegin obj.params.SIMION_nend]; % modify this to modifi x rang of plotting
            ycut = ceil(obj.params.SIMION_nj/2); % modify these to change the y,z cut
            zcut = ceil(obj.params.SIMION_nk/2);
            ycut = 21; zcut=20;
            x = (xrange(1):xrange(2))./obj.params.SIMION_grid_units_p_meter.*1e3; % x axis in mm

            % plot slice along decelerator
            my_xlabel = 'mm (longitudinal axis)'; my_ylabel = 'acc (m/s^2)';
            figure('Name', 'Acceleration slice along decelerator axis') % for x-cuts of accelerations
            subplot(3, 1, 1)
            plot(x, obj.ax_norm(:, ycut, zcut), 'DisplayName', 'ax normal'); hold on;
            xlabel(my_xlabel); ylabel(my_ylabel);
            title('acc along x')
            subplot(3, 1, 2)
            plot(x, obj.ay_norm(:, ycut, zcut), 'DisplayName', 'ay normal'); hold on;
            xlabel(my_xlabel); ylabel(my_ylabel);
            title('acc along y')
            subplot(3, 1, 3)
            plot(x, obj.az_norm(:, ycut, zcut), 'DisplayName', 'ax normal'); hold on;
            xlabel(my_xlabel); ylabel(my_ylabel);
            title('acc along z')
            if obj.params.FLY_focusing_mode_bool == true % plot focusing mode too
                subplot(3, 1, 1); plot(x, obj.ax_pos(:, ycut, zcut), 'DisplayName', 'ax foc +');
                plot(x, obj.ax_neg(:, ycut, zcut), 'DisplayName', 'ax foc -'); legend();
                subplot(3, 1, 2); plot(x, obj.ay_pos(:, ycut, zcut), 'DisplayName', 'ay foc +');
                plot(x, obj.ay_neg(:, ycut, zcut), 'DisplayName', 'ay foc -'); legend();
                subplot(3, 1, 3); plot(x, obj.az_pos(:, ycut, zcut), 'DisplayName', 'az foc +');
                plot(x, abs(obj.az_neg(:, ycut, zcut)), 'DisplayName', 'az foc -'); legend();
            end % Again here ax etc. was adapted no need to adjust range

            % 2D plot
            figure()
            subplot(3, 3, 1); imagesc(obj.ax_norm(:, :, zcut))
            slice_ax_norm_a_y = permute(obj.ax_norm(:, ceil(41), :), [1, 3, 2]); % output is x, z axis
            slice_ay_norm_a_y = permute(obj.ay_norm(:, ceil(41), :), [1, 3, 2]); % output is x, z axis
            slice_az_norm_a_y = permute(obj.az_norm(:, ceil(41), :), [1, 3, 2]); % output is x, z axis
            figure()
            subplot(3, 1, 1)
            imagesc(slice_ax_norm_a_y); colorbar;
            xlabel('x (grid units)'); ylabel('z (grid units)'); title('acc on x (longitudinal)')
            subplot(3, 1, 2)
            imagesc(slice_ay_norm_a_y); colorbar;
            xlabel('x (grid units)'); ylabel('z (grid units)'); title('acc on y (vertical)')
            subplot(3, 1, 3)
            imagesc(slice_az_norm_a_y); colorbar;
            xlabel('x (grid units)'); ylabel('z (grid units)'); title('acc on z (horizonthal)')

            clearvars xrange ycut zcut x
            
            % 3 D vector plot
            figure('Name', '3D vector plot')
            x = 1:111; y = 1:41; z = 1:41; %here x was changed to 111 instead of 151 since our ax.norm is not full size anymore
            [X, Y, Z] = meshgrid(x, y, z); % but before here onlz point where whole ax_norm etc. was used
            X = permute(X, [2, 1, 3]);
            Y = permute(Y, [2, 1, 3]);
            Z = permute(Z, [2, 1, 3]);
            q = quiver3(X, Y, Z, obj.ax_norm, obj.ay_norm, obj.az_norm, 'r-', 'AutoScale', 'on', 'AutoScaleFactor', 30);
            xlabel('x (grid units)'); ylabel('z (grid units)'); title('acc on z (horizonthal)')

            %// Compute the magnitude of the vectors
            mags = sqrt(sum(cat(2, q.UData(:), q.VData(:), ...
                        reshape(q.WData, numel(q.UData), [])).^2, 2));
            
            %// Get the current colormap
            currentColormap = colormap(gca);
            
            %// Now determine the color to make each arrow using a colormap
            [~, ~, ind] = histcounts(mags, size(currentColormap, 1));
            
            %// Now map this to a colormap to get RGB
            cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
            cmap(:,:,4) = 255;
            cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);
            
            %// We repeat each color 3 times (using 1:3 below) because each arrow has 3 vertices
            set(q.Head, ...
                'ColorBinding', 'interpolated', ...
                'ColorData', reshape(cmap(1:3,:,:), [], 4).');   %'
            
            %// We repeat each color 2 times (using 1:2 below) because each tail has 2 vertices
            set(q.Tail, ...
                'ColorBinding', 'interpolated', ...
                'ColorData', reshape(cmap(1:2,:,:), [], 4).');


            fprintf('done\n')
        end


        
        
    end
end