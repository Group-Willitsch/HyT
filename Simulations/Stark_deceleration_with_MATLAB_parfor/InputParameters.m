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
        ax_norm_H_interpl, ay_norm_H_interpl, az_norm_H_interpl
        ax_neg_extended, ay_neg_extended, az_neg_extended
        ax_neg_interpl,ay_neg_interpl,az_neg_interpl
        ax_neg_H_interpl, ay_neg_H_interpl, az_neg_H_interpl
        ax_pos_interpl,ay_pos_interpl,az_pos_interpl
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
        electrode_sequences % variable to save the integer arrays containg the sequence for each of the electrode
        
        num_particles
        xyzVxyz_0                       % the vector of initial pos&vel, number_of_particles * 6
        xyzVxyz                         % the vector of computed pos&vel, number_of_particles * 6
        has_the_simulation_been_run     % boolean, default is False in constructor, will be switched to True at the end of the simulation
        num_trajectories_saved          % number of trajectories to be saved
        arrival_time                    % TODO; WTF? this variable is defined in some plotting function?? Maybe better just a lcoal variable
        output                          % save xyzVxyz , t and flag at certain times and display here
        traj_xyzVxyz                    % trajectories x,y,z, Vx, Vy, Vz
        traj_time                       % ?????
        ind_particles                 % index of particles created at begining saving the indexes of the ones that make it to the end
        TOF_xyzVxyz                   % Time of Flight: variable to save xyzVxyz at the time step/steps of detetcion
        TOF_save     
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
                options.FortranSeqBool      logical = false % default behaviour for Fortran sequences
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
            obj.params.BEAM_avg_velocity_beam = 480; % Average velocity incoming package (m/s) 
            obj.params.BEAM_radius_of_nozzle = 0.25e-3; % radius of the valve nozzle (m)
            obj.params.BEAM_long_pos_spread = 11.5e-3; % ????? longitudinal position spread (m) - along x aixs or beam propagation
            obj.params.BEAM_long_vel_spread = 0.13*obj.params.BEAM_avg_velocity_beam ; % relative velocity spread along beam axis 0.112 0.12
            obj.params.BEAM_trans_velocity_spread = 5; % will be 15 or 25 in the end, arbitrary was not measured velocity spread perp. to beam axis/average velocity 0.10

            % Parameters used in fly and in the simulation
            obj.params.FLY_incoupling_time = obj.params.PHYS_valve_to_dec/obj.params.CALC_vel_synch_mol; %710.2e-6; % valve - decelerator incoupling time (s)
            obj.params.FLY_detection_laser_diameter = 1e-3;
            obj.params.FLY_simulated_target_vel = []; % to be assigned while loading or generating the Matlab sequence.
            % Stores the rounded value of the last velocity vector. Ideally
            % equals obj.params.FLY_target_velocity

            % load the acceleration fields
            tic;obj.loadAccelerationFields();toc;
            obj.InterpolateAccField()
            
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
        end

        function setPhase(obj, new_phase) % mandatory to update the phase distance too!!
            obj.params.CALC_phase_degrees = new_phase;
            obj.params.CALC_phase_distance = obj.params.CALC_phase_degrees * obj.params.PHYS_distance_stages / 180;
        end
        
        %% Load the acceleration fields
        % (TODO: maybe the negative one can be obtained via a flip of the postive
        % one)
        % updated to load directly the .mat files
        function loadAccelerationFields(obj)           
            if ~obj.params.FLY_focusing_mode_bool % load normal fields
                if ~isempty(obj.ax_norm)
                    [obj.ax_norm, obj.ay_norm, obj.az_norm] = deal([],[],[]);
                end    
                fprintf('Loading normal mode fields ...\t')
                
                % set the focusing mode fields to empty. Necessary if you
                % first load focusing mode fields and after you re-load
                % normal mode fields.
                if not(isempty(obj.ax_pos))
                    [obj.ax_pos, obj.ax_neg, obj.ay_pos, obj.ay_neg, obj.az_pos, obj.az_neg] = deal([], [], [], [], [], []);
                end
                
                % load fields from the .mat files in the acc folder
                filename = './acc/dec_norm_' + strrep( string( obj.params.FLY_voltage_on_electrodes ), '.', 'p') + 'kV/';
                load(filename + 'a_norm') % load them in local workspace
                obj.ax_norm = ax_norm; obj.ay_norm = ay_norm; obj.az_norm = az_norm;
                clearvars ax_norm ay_norm az_norm % ugly: clear them from local workspace

                % symmetrize the y, z fields: has already been down when saving field mat files, see "generateMatAccFiles.m"
%                  obj.ay_norm = (obj.ay_norm - flip(obj.ay_norm, 2))/2;
%                  obj.az_norm = (obj.az_norm - flip(obj.az_norm, 3))/2;
                
            else % load focusing mode fields
                fprintf('Loading focusing mode accelerations...')

                if ~isempty(obj.ax_pos)
                    [obj.ax_pos, obj.ax_neg, obj.ay_pos, obj.ay_neg, obj.az_pos, obj.az_neg] = deal([], [], [], [], [], []);
                end
                filename = './acc/dec_foc_' + strrep( string(obj.params.FLY_voltage_on_electrodes), '.', 'p') + 'kV/';
                load(filename + 'a_norm'); load(filename + 'a_pos'); load(filename + 'a_neg'); 
                obj.ax_norm = ax_norm; obj.ay_norm = ay_norm; obj.az_norm = az_norm;
                obj.ax_pos = ax_pos; obj.ay_pos = ay_pos; obj.az_pos = az_pos;
                obj.ax_neg = ax_neg; obj.ay_neg = ay_neg; obj.az_neg = az_neg;

                clearvars ax_norm ax_pos ax_neg ay_norm ay_neg ay_pos az_norm az_pos az_neg
               
                
                % symmetrize the y, z fields: has already been down when saving field mat files, see "generateMatAccFiles.m"
%                  obj.ay_norm = (obj.ay_norm - flip(obj.ay_norm, 2))/2;
%                  obj.az_norm = (obj.az_norm - flip(obj.az_norm, 3))/2;
%                   % also for focusing mode
%                  obj.ay_pos = (obj.ay_pos + flip(obj.ay_neg, 3))/2;
%                  obj.ay_neg = flip( obj.ay_pos, 3);
%                  obj.az_pos = (obj.az_pos - flip(obj.az_neg, 3))/2;
%                  obj.az_neg = - flip(obj.az_pos, 3);
            end % when calling this function we have to also call interpolateAccField again in order 
            fprintf('\tloaded\n')
        end
        
        %% Interpolate acceleration field
        function InterpolateAccField(obj)
            num_grids_x = 111 + 110 * 122 + 4;
            num_grids_y = 41 + 4;
            num_grids_z = 41 + 4;
            gridded_x = linspace(-2/obj.params.SIMION_grid_units_p_meter, obj.params.PHYS_length_dec + 2/obj.params.SIMION_grid_units_p_meter, num_grids_x);
            gridded_y = linspace(-2/obj.params.SIMION_grid_units_p_meter-obj.params.PHYS_seperation_pins/2.0, obj.params.PHYS_seperation_pins/2.0 + 2/obj.params.SIMION_grid_units_p_meter, num_grids_y);
            gridded_z =	linspace(-2/obj.params.SIMION_grid_units_p_meter-obj.params.PHYS_seperation_pins/2.0, obj.params.PHYS_seperation_pins/2.0 + 2/obj.params.SIMION_grid_units_p_meter, num_grids_z);

            obj.ax_norm_extended = zeros(num_grids_x, num_grids_y, num_grids_z); %Vertical ones
            obj.ay_norm_extended = zeros(num_grids_x, num_grids_y, num_grids_z);
            obj.az_norm_extended = zeros(num_grids_x, num_grids_y, num_grids_z);
            obj.ax_norm_extended(3:num_grids_x-2,3:43,3:43) = cat(1, repmat(cat(1,obj.ax_norm, flip(-obj.ax_norm(2:110,:,:),1)), 61, 1, 1), obj.ax_norm);  % before  obj.ax_norm_extended(11:num_grids_x-10,11:51,11:51) = cat(1, repmat(cat(1,obj.ax_norm, flip(-obj.ax_norm(22:130,:,:),1)), 61, 1, 1), obj.ax_norm);
            obj.ay_norm_extended(3:num_grids_x-2,3:43,3:43) = cat(1, repmat(cat(1,obj.ay_norm, flip(obj.ay_norm(2:110,:,:),1)), 61, 1, 1), obj.ay_norm);
            obj.az_norm_extended(3:num_grids_x-2,3:43,3:43) = cat(1, repmat(cat(1,obj.az_norm, flip(obj.az_norm(2:110,:,:),1)), 61, 1, 1), obj.az_norm);
            
            obj.ax_norm_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.ax_norm_extended,'linear','linear');
            obj.ay_norm_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.ay_norm_extended,'linear','linear');
            obj.az_norm_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.az_norm_extended,'linear','linear');
            
            obj.ax_norm_H_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, -flip(obj.ax_norm_extended,1),'linear','linear');
            obj.ay_norm_H_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, flip(obj.ay_norm_extended,1),'linear','linear');
            obj.az_norm_H_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, flip(obj.az_norm_extended,1),'linear','linear');
            
            if obj.params.FLY_focusing_mode_bool
                obj.ax_neg_extended = zeros(num_grids_x, num_grids_y, num_grids_z);
                obj.ay_neg_extended = zeros(num_grids_x, num_grids_y, num_grids_z);
                obj.az_neg_extended = zeros(num_grids_x, num_grids_y, num_grids_z);
                obj.ax_neg_extended(3:num_grids_x-2,3:43,3:43) = cat(1, repmat(cat(1,obj.ax_neg, flip(-obj.ax_neg(2:110,:,:),1)), 61, 1, 1), obj.ax_neg);
                obj.ay_neg_extended(3:num_grids_x-2,3:43,3:43) = cat(1, repmat(cat(1,obj.ay_neg, flip(obj.ay_neg(2:110,:,:),1)), 61, 1, 1), obj.ay_neg);
                obj.az_neg_extended(3:num_grids_x-2,3:43,3:43) = cat(1, repmat(cat(1,obj.az_neg, flip(obj.az_neg(2:110,:,:),1)), 61, 1, 1), obj.az_neg);
                
                obj.ax_neg_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.ax_neg_extended,'linear','linear');
                obj.ay_neg_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.ay_neg_extended,'linear','linear');
                obj.az_neg_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.az_neg_extended,'linear','linear');

                obj.ax_neg_H_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, -flip(obj.ax_neg_extended,1),'linear','linear');
                obj.ay_neg_H_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, flip(obj.ay_neg_extended,1),'linear','linear');
                obj.az_neg_H_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, flip(obj.az_neg_extended,1),'linear','linear');
                
                obj.ax_pos_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, flip(obj.ax_neg_extended,3),'linear','linear');
                obj.ay_pos_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, flip(obj.ay_neg_extended,3),'linear','linear');
                obj.az_pos_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, -flip(obj.az_neg_extended,3),'linear','linear');

                obj.ax_pos_H_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, -flip(flip(obj.ax_neg_extended,3),1),'linear','linear');
                obj.ay_pos_H_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, flip(flip(obj.ay_neg_extended,3),1),'linear','linear');
                obj.az_pos_H_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, flip(-flip(obj.az_neg_extended,3),1),'linear','linear');
            end    
        end

        
        %% loadFortranTimeSequence from T2jump.out
        % from appropriate folder and save it to the class variables
        % MANDATORY TO FIRST RUN THE SH CODE THAT CLEANS UP THE FILES...
        % matlab goes nut for the variable number of whitespaces. Python
        % should be fine with it. 
        function loadFortranTimeSequence(obj)
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

        function [simulated_target_vel] = generateMatlabTimeSequence(obj)
            fprintf('Generating Matlab time sequence...');

            % 111 pops out from?
            %changed ax norm again since not cut anymore
            obj.ax_norm_1d_interpl = griddedInterpolant( linspace(0, obj.params.PHYS_distance_stages, 111), obj.ax_norm(:,21,21),'linear','none'); %Here we make new a_x values in the center of the trap at the trap and at the stages?electrodes?           
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
            if use_ode_solver_bool % use ode45 from Matlab
                opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-8, 'Events', @(t,x) EventsFcn(t,x));
                [obj.M_synch_time, x_Vx_temp] = ode45( @(t,x) ...
                    obj.dxdt(t,x), [0, 5e-3], [0; obj.params.CALC_vel_synch_mol], opts); % ode23t seems a good solver so why ode45 used?
            else % use ugly hand made Euler method
                xx = [0; obj.params.CALC_vel_synch_mol];
                t_step = 1e-8;
                tt= 0:t_step:5e-3;
                x_Vx = zeros(length(tt),2);
                for i = 1:1:length(tt)
                    if xx(1) < obj.params.PHYS_length_dec + obj.params.PHYS_exit_to_detection
                        x_Vx(i,:) = xx;
                        xx = xx + obj.dxdt(tt(i), xx)*t_step; %eulers method to update pos for each time step till pos is at end of decc.
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
            
            if obj.M_synch_velocity(end) < 0
                error("Sychronous molecule is bounced back, please lower the phase angle!")
            end

            t_x_interpl = griddedInterpolant(obj.M_synch_position, obj.M_synch_time); % to obatain the exact field switching time based on the switching position (pos synchronus molecule)
            if obj.params.FLY_focusing_mode_bool
                positions_to_jump = [0, union(obj.params.PHYS_distance_stages/2 + obj.params.CALC_phase_distance: obj.params.PHYS_distance_stages:... %union combines data into array
                    obj.params.PHYS_length_dec-(obj.params.PHYS_distance_stages/2 - obj.params.CALC_phase_distance), obj.params.PHYS_distance_stages*3.0/2.0 -... 
                    obj.params.CALC_phase_distance: obj.params.PHYS_distance_stages:obj.params.PHYS_length_dec - (obj.params.PHYS_distance_stages/2 + obj.params.CALC_phase_distance)), ...
                    obj.params.PHYS_length_dec + obj.params.PHYS_exit_to_detection]'; % list of positions of interest
                obj.M_time_vec = t_x_interpl(positions_to_jump);  %B = repmat(A,n) returns an array containing n copies of A in the row and column dimensions. The size of B is size(A)*n when A is a matrix.
                trigger_pattern = repmat(["b0011";"b0100";"b1100";"b0010";"b0011";"b1000";"b1100";"b0001"], obj.params.PHYS_number_of_electrodes,1);
                obj.M_trigger_pattern = trigger_pattern(1:length(obj.M_time_vec));
                obj.M_trigger_pattern(end-1:end)= "b0000"; %what is done here?
            else
                positions_to_jump = [0, obj.params.PHYS_distance_stages/2 + obj.params.CALC_phase_distance: obj.params.PHYS_distance_stages:...
                    obj.params.PHYS_length_dec-(obj.params.PHYS_distance_stages/2 - ...
                    obj.params.CALC_phase_distance), obj.params.PHYS_length_dec + obj.params.PHYS_exit_to_detection]'; % list of positions of interest
                obj.M_time_vec = t_x_interpl(positions_to_jump);
                trigger_pattern = repmat(["0x0010";"0x0020"], obj.params.PHYS_number_of_electrodes,1);
                obj.M_trigger_pattern = trigger_pattern(1:length(obj.M_time_vec));
                obj.M_trigger_pattern(end-1:end)= "0x0000";
            end
            fprintf("Precise final velocity of sequence is %s", num2str(obj.M_synch_velocity(end)))
            obj.params.FLY_simulated_target_vel = round( obj.M_synch_velocity(end), 0);
            % stores in a quite useless but safe variable the effective final velocity
            % Will be used to create the filename

            obj.saveMatlabSequence(); % saves both .out and .m files for the sequence just created
            % in the .out file we save M_time_vec and M_trigger_pattern
            % in the .mat file we save all the files with M_'something'

            if obj.verbose % we plot the synch molecule
                obj.plotTimeSequence()
            end

            fprintf('\tFinal Matlab velocity is %d\n', obj.params.FLY_simulated_target_vel)
            simulated_target_vel = obj.params.FLY_simulated_target_vel; 
%             I return the compute velocity, to be used if needed.
%             for some reason I cannot return the class variable but I must
%             declare a local one. Probably due to handle and class stuff.
        end
        
        function dxdt = dxdt(obj, t, x)
            if obj.params.FLY_focusing_mode_bool 
                if x(1) > obj.params.PHYS_length_dec -(obj.params.PHYS_distance_stages/2 - obj.params.CALC_phase_distance)
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
                    '_' + num2str(obj.params.FLY_target_velocity) + '/T2jump.out' ;
            else % folder path for focusing files
                timeSequenceFolder = ...
                    '../dec_FM_' + strrep( string(obj.params.FLY_voltage_on_electrodes), '.', 'p') + ...
                    'kV/timeseq/sequences/dec_' + num2str(obj.params.CALC_vel_synch_mol) + ...
                    '_' + num2str(obj.params.FLY_target_velocity) + '/T2jump.out';
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
            fprintf(seq_file,"#vx_i=%.1f m/s\n#vx_f=%.1f m/s\n#phase=%.2f deg\n#\n",obj.params.FLY_target_velocity, obj.M_synch_velocity(end), obj.params.CALC_phase_degrees);
            fprintf(seq_file, '%s\t%s\n', [ string( round( obj.M_time_vec*1e9)+1010 ), obj.M_trigger_pattern]' );
            fprintf(seq_file,"#\n#\n#\n");
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

            if obj.verbose
                obj.plotTimeSequence()
            end    


        end

        %% changeFieldConfig
        % This function re-loads the fields and must be used whenever
        % you want to change voltage values or focusing/normal mode
        function changeFieldConfig(obj, new_voltage, new_focusing_mode_bool) %There is a problem when chaning from NM to FM but not other way around
            if nargin ~= 3 
                fprintf('Wrong changeFieldConfig call.\nUsage changeFieldConfig( new_voltage, new_focusing_mode_bool)\n')
                return
            else
            % overwrite parameters
                obj.params.FLY_voltage_on_electrodes = new_voltage;
                obj.params.FLY_focusing_mode_bool = new_focusing_mode_bool;
            end
            if ~isempty(obj.ax_norm) % set fields to zero, for safety
                 [obj.ax_norm, obj.ay_norm, obj.az_norm] = deal([],[],[]);
            end
            if ~isempty(obj.ax_pos)
                 [obj.ax_pos, obj.ax_neg, obj.ay_pos, obj.ay_neg, obj.az_pos, obj.az_neg] = deal([], [], [], [], [], []);
            end
            % re-load acceleration fields
            obj.loadAccelerationFields();
            obj.InterpolateAccField()

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
%             fprintf('TODO: to be tested properly \n')
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
%             fprintf('TODO: to be tested properly \n')
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
                plot(double(obj.T2Jump_time_vec(1:end-2)) - obj.M_time_vec(1:end-1) *1e9, '-o', 'DisplayName', 'difference'); 
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

            xrange = [obj.params.SIMION_nbegin obj.params.SIMION_nend];
            xcut=ceil(111/2);
            ycut = ceil(obj.params.SIMION_nj/2); % modify these to change the y,z cut
            zcut = ceil(obj.params.SIMION_nk/2);

          
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
            plot(x, obj.az_norm(:, ycut, zcut), 'DisplayName', 'az normal'); hold on;
            xlabel(my_xlabel); ylabel(my_ylabel);
            title('acc along z')
            if obj.params.FLY_focusing_mode_bool == true % plot focusing mode too
                subplot(3, 1, 1); plot(x, obj.ax_pos(:, ycut, zcut), 'DisplayName', 'ax foc +');
                plot(x, obj.ax_neg(:, ycut, zcut), 'DisplayName', 'ax foc -'); legend();
                subplot(3, 1, 2); plot(x, obj.ay_pos(:, ycut, zcut), 'DisplayName', 'ay foc +');
                plot(x, obj.ay_neg(:, ycut, zcut), 'DisplayName', 'ay foc -'); legend();
                subplot(3, 1, 3); plot(x, obj.az_pos(:, ycut, zcut), 'DisplayName', 'az foc +');
                plot(x, obj.az_neg(:, ycut, zcut), 'DisplayName', 'az foc -'); legend();
            end % Again here ax etc. was adapted no need to adjust range

            % 2D plot

            slice_ax_norm_a_x= permute(obj.ax_norm(xcut, :, :), [3, 2, 1]);
            slice_ax_norm_a_y= permute(obj.ax_norm(:, ycut, :), [1, 3, 2]);
            slice_ax_norm_a_z= permute(obj.ax_norm(:, :, zcut), [1, 2, 3]);
            slice_ay_norm_a_x = permute(obj.ay_norm(xcut, :, :), [3, 2, 1]);
            slice_ay_norm_a_y = permute(obj.ay_norm(:, ycut, :), [1, 3, 2]);
            slice_ay_norm_a_z = permute(obj.ay_norm(:,:, zcut), [1, 2, 3]);
            slice_az_norm_a_x = permute(obj.az_norm(xcut, :, :), [3, 2, 1]);
            slice_az_norm_a_y = permute(obj.az_norm(:,ycut, :), [1, 3, 2]);
            slice_az_norm_a_z = permute(obj.az_norm(:, :, zcut), [1, 2, 3]);

            figure()
            subplot(3, 3, 1)
            imagesc(slice_ax_norm_a_x); colorbar; axis xy;
            xlabel('z (grid units)'); ylabel('y (grid units)'); title('a_x, z-y (longitudinal)')
            subplot(3, 3, 2)
            imagesc(slice_ax_norm_a_y); colorbar; axis xy;
            xlabel('x (grid units)'); ylabel('z (grid units)'); title('a_x, x-z (longitudinal)')
            subplot(3, 3, 3)
            imagesc(slice_ax_norm_a_z); colorbar; axis xy;
            xlabel('x(grid units)'); ylabel('y (grid units)'); title('a_x, x-y (longitudinal)')
            subplot(3, 3, 4)
            imagesc(slice_ay_norm_a_x); colorbar; axis xy;
            xlabel('z (grid units)'); ylabel('y (grid units)'); title('a_y, z-y (longitudinal)')
            subplot(3, 3, 5)
            imagesc(slice_ay_norm_a_y); colorbar; axis xy;
            xlabel('x (grid units)'); ylabel('z (grid units)'); title('a_y, x-z (longitudinal)')
            subplot(3, 3, 6)
            imagesc(slice_ay_norm_a_z); colorbar; axis xy;
            xlabel('x (grid units)'); ylabel('y (grid units)'); title('a_y, x-y (longitudinal)')
            subplot(3, 3, 7)
            imagesc(slice_az_norm_a_x); colorbar; axis xy;
            xlabel('z (grid units)'); ylabel('y (grid units)'); title('a_z, z-y (longitudinal)')
            subplot(3, 3, 8)
            imagesc(slice_az_norm_a_y); colorbar; axis xy;
            xlabel('x (grid units)'); ylabel('z (grid units)'); title('a_z, x-z (longitudinal)')
            subplot(3, 3, 9)
            imagesc(slice_az_norm_a_z); colorbar; axis xy;
            xlabel('x (grid units)'); ylabel('y (grid units)'); title('a_z, x-y (longitudinal)')

            if obj.params.FLY_focusing_mode_bool == true

                slice_ax_pos_a_x = permute(obj.ax_pos(xcut,:,:), [3, 2, 1]);
                slice_ax_pos_a_y = permute(obj.ax_pos(:, ycut,:), [3, 1, 2]);
                slice_ax_pos_a_z = permute(obj.ax_pos(:,:,zcut), [2, 1, 3]);
                
                slice_ax_neg_a_x = permute(obj.ax_neg(xcut,:,:), [3, 2, 1]);
                slice_ax_neg_a_y = permute(obj.ax_neg(:, ycut,:), [3, 1, 2]);
                slice_ax_neg_a_z = permute(obj.ax_neg(:,:,zcut), [2, 1, 3]);

                figure()
                subplot(2,3,1)
                imagesc(slice_ax_pos_a_x); colorbar; axis xy;
                xlabel('z (grid units)'); ylabel('y (grid units)'); title('acc ax_{pos}, z-y');

                subplot(2,3,2)
                imagesc(slice_ax_pos_a_y); colorbar; axis xy;
                xlabel('x (grid units)'); ylabel('z (grid units)'); title('acc ax_{pos}, x-z');

                subplot(2,3,3)
                imagesc(slice_ax_pos_a_z); colorbar; axis xy;
                xlabel('x (grid units)'); ylabel('y (grid units)'); title('acc ax_{pos}, x-y');

                subplot(2,3,4)
                imagesc(slice_ax_neg_a_x); colorbar; axis xy;
                xlabel('z (grid units)'); ylabel('x (grid units)'); title('acc ax_{neg}, z-y');

                subplot(2,3,5)
                imagesc(slice_ax_neg_a_y); colorbar; axis xy;
                xlabel('x (grid units)'); ylabel('z (grid units)'); title('acc ax_{neg}, x-z');

                subplot(2,3,6)
                imagesc(slice_ax_neg_a_z); colorbar; axis xy;
                xlabel('x (grid units)'); ylabel('y (grid units)'); title('acc ax_{neg}, x-y');


                slice_ay_pos_a_x = permute(obj.ay_pos(xcut,:,:), [3, 2, 1]);
                slice_ay_pos_a_y = permute(obj.ay_pos(:, ycut,:), [3, 1, 2]);
                slice_ay_pos_a_z = permute(obj.ay_pos(:,:,zcut), [2, 1, 3]);
                
                slice_ay_neg_a_x = permute(obj.ay_neg(xcut,:,:), [3, 2, 1]);
                slice_ay_neg_a_y = permute(obj.ay_neg(:, ycut,:), [3, 1, 2]);
                slice_ay_neg_a_z = permute(obj.ay_neg(:,:,zcut), [2, 1, 3]);

                figure()
                subplot(2,3,1)
                imagesc(slice_ay_pos_a_x()); colorbar; axis xy;
                xlabel('z (grid units)'); ylabel('y (grid units)'); title('acc ay_{pos}, z-y');

                subplot(2,3,2)
                imagesc(slice_ay_pos_a_y); colorbar; axis xy; 
                xlabel('x (grid units)'); ylabel('z (grid units)'); title('acc ay_{pos}, x-z');

                subplot(2,3,3)
                imagesc(slice_ay_pos_a_z); colorbar; axis xy;
                xlabel('x (grid units)'); ylabel('y (grid units)'); title('acc ay_{pos}, x-y');

                subplot(2,3,4)
                imagesc(slice_ay_neg_a_x); colorbar; axis xy;
                xlabel('z (grid units)'); ylabel('y (grid units)'); title('acc ay_{neg}, z-y');

                subplot(2,3,5)
                imagesc(slice_ay_neg_a_y); colorbar; axis xy; 
                xlabel('x (grid units)'); ylabel('z (grid units)'); title('acc ay_{neg}, x-z');

                subplot(2,3,6)
                imagesc(slice_ay_neg_a_z); colorbar; axis xy;
                xlabel('x (grid units)'); ylabel('y (grid units)'); title('acc ay_{neg}, x-y');


                slice_az_pos_a_x = squeeze(obj.az_pos(xcut,:,:));
                slice_az_pos_a_y = permute(obj.az_pos(:, ycut,:), [3, 1, 2]);
                slice_az_pos_a_z = permute(obj.az_pos(:,:,zcut), [2, 1, 3]);
                
                slice_az_neg_a_x = permute(obj.az_neg(xcut,:,:), [3, 2, 1]);
                slice_az_neg_a_y = permute(obj.az_neg(:, ycut,:), [3, 1, 2]);
                slice_az_neg_a_z = permute(obj.az_neg(:,:,zcut), [2, 1, 3]);


                figure()
                subplot(2,3,1)
                imagesc(slice_az_pos_a_x); colorbar; axis xy;
                xlabel('z (grid units)'); ylabel('x (grid units)'); title('acc az_{pos}, z-y');

                subplot(2,3,2)
                imagesc(slice_az_pos_a_y); colorbar; axis xy; caxis([0,1.7e5])
                xlabel('x (grid units)'); ylabel('z (grid units)'); title('acc az_{pos}, x-z');

                subplot(2,3,3)
                imagesc(slice_az_pos_a_z); colorbar; axis xy;
                xlabel('x (grid units)'); ylabel('y (grid units)'); title('acc az_{pos}, x-y');

                subplot(2,3,4)
                imagesc(slice_az_neg_a_x); colorbar; axis xy;
                xlabel('z (grid units)'); ylabel('y (grid units)'); title('acc az_{neg}, z-y');

                subplot(2,3,5)
                imagesc(slice_az_neg_a_y); colorbar; axis xy; caxis([-1.7e5,0])
                xlabel('x (grid units)'); ylabel('z (grid units)'); title('acc az_{neg}, x-z');

                subplot(2,3,6)
                imagesc(slice_az_neg_a_z); colorbar; axis xy;
                xlabel('x (grid units)'); ylabel('y (grid units)'); title('acc az_{neg}, x-y');
            end    


            clearvars xrange ycut zcut x
            
            % 3 D vector plot
            figure('Name', '3D vector plot')
            x = 1:80; y = 1:41; z = 1:41; %here x was changed to 111 instead of 151 since our ax.norm is not full size anymore
            [X, Y, Z] = meshgrid(x, y, z); % but before here only point where whole ax_norm till 151 etc. was used
            X = permute(X, [2, 1, 3]);
            Y = permute(Y, [2, 1, 3]);
            Z = permute(Z, [2, 1, 3]);
            q = quiver3(X, Y, Z, obj.ax_norm(1:80,:,:), obj.ay_norm(1:80,:,:), obj.az_norm(1:80,:,:),0.6, 'r-', 'ShowArrowHead', 'on');
            xlabel('x (grid units)'); ylabel('y (grid units)'); zlabel('z (grid units)'); title('Vector plot of accelerations ax,ay,az NM')

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

        function plotTimeSequence(obj)
            pattern=char(obj.M_trigger_pattern); % ugly code to turn M_pattern into array of doubles such that we can compare the pattern such that we get a string of ones and zeros           
            pattern=double(string(pattern(:,end-3:end))); % the same length as M_pattern where 1 means the correspondign electrode is on or off at this specific time     

            if  obj.params.FLY_focusing_mode_bool            
                rod1= double(pattern== 1000 | pattern==1100)+1; % numbering of electrtodes follows b1100 sample where electrode corresponds to the electrode represented by the first number of b1100 etc           
                rod2= double(pattern== 100 | pattern==1100)+3; % we add numbers to put the on same plot but still be able to see each individual sequence            
                rod3= double(pattern== 10 | pattern==11)+5; % to get meaning of the numbers b1000=1000, b1100=1100, b0011= 11, b0010=10, b0001=1            
                rod4= double(pattern== 1 | pattern==11)+7;


               
                title('Trigger sequence electrodes (Channel i labframe/other one)')
                hold on
                stairs(obj.M_time_vec*10^3,rod1) 
                ylim([0,9]); xlabel('time(ms)');
                stairs(obj.M_time_vec*10^3,rod2)
                stairs(obj.M_time_vec*10^3,rod3)
                stairs(obj.M_time_vec*10^3,rod4)
                legend('ch.8 (G+/H+)','ch.9 (E-/H-)', 'ch.6 (H+/V+)', 'ch.7 (F-/V-)');
                hold off
            else
                rod12= double(pattern== 10)+1;
                rod23= double(pattern== 20)+3;
                figure()
                title('Trigger sequence electrodes (labframe)')
                hold on
                stairs(obj.M_time_vec*10^3,rod12) 
                ylim([0,5]); xlabel('time(ms)');
                stairs(obj.M_time_vec*10^3,rod23)
                legend('vertical el. (H+,F-)','horizontal el. (E-,G+)');
                hold off
            end

             figure('Name', 'Time seq with Matlab');

                subplot(2, 3, 1)    % x vs t
                plot(obj.M_synch_time .* 1e3, obj.M_synch_position);
                title('Position vs time'); ylabel('x (m)'); xlabel('t (ms)')

                % linear trendline, see below
%                 linear_velocity = obj.params.CALC_vel_synch_mol - ...
%                         (obj.params.CALC_vel_synch_mol - ...21:131
%                         obj.M_synch_velocity(end)) / ...
%                         obj.M_synch_time(end) .* obj.M_synch_time;
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
                ax=diff(obj.M_synch_velocity)./diff(obj.M_synch_time);
                ax=[0;ax];
                plot(obj.M_synch_time .*1e3,ax, '-k.', 'LineWidth', 0.5);
                title('ax (m/s^2) vs t(ms)'); ylabel('ax (m/s62)'); xlabel('t (ms)');
                % that follows V(t) =  V_0 - (V_0 - V_final) / t_final * t
                % whre V_0 = starting velocity, V_final t_final velocity
                % and time of the last instant of the simulation
%                plot(obj.M_synch_time .* 1e3, obj.M_synch_velocity - linear_velocity);
%                 title('Vx (m/s) - linear decrease'); ylabel('Vx (m/s)'); xlabel('t (ms)');

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
                histogram(time_step_difference, 300)
                xlabel('Intergation timestep (ns)'); title('Histrogram of time steps')

        end
        
          %% create particles
%         
        function createParticles(obj)
            rng('default') % fix the seed of the random nunber generator to bremoved afterwards
            obj.xyzVxyz_0 = randn(obj.num_particles, 6)/sqrt(8*log(2)).*... %conversion factor to std from full width half max, below adjust normal dist. to represent particles
                [obj.params.BEAM_long_pos_spread, obj.params.BEAM_radius_of_nozzle, obj.params.BEAM_radius_of_nozzle,...                            
                 obj.params.BEAM_long_vel_spread, obj.params.BEAM_trans_velocity_spread, obj.params.BEAM_trans_velocity_spread] + ...
                [-obj.params.PHYS_valve_to_dec, 0, 0, obj.params.BEAM_avg_velocity_beam, 0, 0]; %set v_x to avergae and set all particles to valve pos. (decc x=0)

            obj.xyzVxyz_0(1,:) = [-obj.params.PHYS_valve_to_dec, 0, 0, 450, 0, 0]; % The first row is for a synchronous molecule
%             E_kin=((0.5*((obj.xyzVxyz_0(:,4)).^2 + (obj.xyzVxyz_0(:,5)).^2 + (obj.xyzVxyz_0(:,6)).^2))-0.5*(obj.xyzVxyz_0(1,4).^2 + obj.xyzVxyz_0(1,5).^2 + obj.xyzVxyz_0(1,6).^2))*1e-5;           
%             sort_E_kin= [obj.xyzVxyz_0,abs(E_kin)];
%             sort_E_kin= sortrows(sort_E_kin,7);
%             obj.xyzVxyz_0=sort_E_kin(:,1:6);   
        end

%% Propagate particles, intergate using time intervals in M_time_vector, propagating all the molecules together
        function propagateParticles_euler(obj)
%             obj.num_particles = 100000;
            obj.createParticles();
            
            ax_norm_interpl = obj.ax_norm_interpl;
            ay_norm_interpl = obj.ay_norm_interpl;
            az_norm_interpl = obj.az_norm_interpl;
            
            ax_norm_H_interpl =  obj.ax_norm_H_interpl;
            ay_norm_H_interpl = obj.ay_norm_H_interpl;
            az_norm_H_interpl = obj.az_norm_H_interpl;
            
            if obj.params.FLY_focusing_mode_bool
                 ax_neg_interpl = obj.ax_neg_interpl;
                 ay_neg_interpl = obj.ay_neg_interpl;
                 az_neg_interpl = obj.az_neg_interpl;
                 ax_neg_H_interpl = obj.ax_neg_H_interpl;
                 ay_neg_H_interpl = obj.ay_neg_H_interpl;
                 az_neg_H_interpl = obj.az_neg_H_interpl;
                 ax_pos_interpl = obj.ax_pos_interpl;
                 ay_pos_interpl = obj.ay_pos_interpl;
                 az_pos_interpl = obj.az_pos_interpl;
                 ax_pos_H_interpl = obj.ax_pos_H_interpl;
                 ay_pos_H_interpl = obj.ay_pos_H_interpl;
                 az_pos_H_interpl = obj.az_pos_H_interpl;
            end
            
            
%             a function that removes the lost molecules
            function [xyzVxyz,ind_particles] = removeHitParticles(xyzVxyz,ind_particles)
                hit_indices = xyzVxyz(:,1) < obj.params.PHYS_length_dec & (abs(xyzVxyz(:,2)) > obj.params.PHYS_seperation_pins/2 | abs(xyzVxyz(:,3)) > obj.params.PHYS_seperation_pins/2);
                xyzVxyz = xyzVxyz(~hit_indices,:);
                ind_particles = ind_particles(~hit_indices);
                hit_indices = (abs(xyzVxyz(:,1) - xyzVxyz(1,1)) < 2* 5.5e-3); % separated the cut of the fish in x into two lines in order to make indexing work and also to be sure its correct
                xyzVxyz = xyzVxyz(hit_indices, :);
                ind_particles = ind_particles(hit_indices);
            end
            
            % first propagate to entrance of decelerator
            obj.xyzVxyz = obj.xyzVxyz_0;
            obj.ind_particles = transpose(1:size(obj.xyzVxyz,1));
            obj.xyzVxyz(:,1:3) = obj.xyzVxyz(:,1:3) + obj.xyzVxyz(:,4:6) * obj.params.FLY_incoupling_time; % ekin at start of deacc.
            [obj.xyzVxyz,obj.ind_particles]= removeHitParticles(obj.xyzVxyz, obj.ind_particles );% select those that can enter dec
%             hit_indices = obj.xyzVxyz(:,1) < obj.params.PHYS_length_dec & (abs(obj.xyzVxyz(:,2)) > obj.params.PHYS_seperation_pins/2 | abs(obj.xyzVxyz(:,3)) > obj.params.PHYS_seperation_pins/2);
%             obj.xyzVxyz = obj.xyzVxyz(~hit_indices,:);
%             obj.xyzVxyz = obj.xyzVxyz((abs(obj.xyzVxyz(:,1) - obj.xyzVxyz(1,1)) < 2* 5.5e-3), :);

            % propagate inside dec
            
            dt = 4e-8;
%             fprintf("num/total switching\n");
            obj.Snapshot("start decelerator",obj.M_time_vec(1),obj.xyzVxyz, obj.ind_particles)
            if obj.params.FLY_focusing_mode_bool % focusing mode
                dxyzVxyz = {@(y) [y(:, 4:6), ax_norm_interpl(y(:, 1), y(:, 2), y(:, 3)),...
                                                 [0;ay_norm_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))],...
                                                 [0;az_norm_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))]]                   
                            @(y) [y(:, 4:6), ax_neg_H_interpl(y(:, 1), y(:, 3), -y(:, 2)),...
                                                 [0;-az_neg_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))],...
                                                 [0;ay_neg_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))]]
                            @(y) [y(:, 4:6), ax_norm_H_interpl(y(:, 1), y(:, 3), -y(:, 2)),...
                                                 [0;-az_norm_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))],...
                                                 [0;ay_norm_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))]]
                            @(y) [y(:, 4:6), ax_pos_interpl(y(:, 1), y(:, 2), y(:, 3)),...
                                                 [0;ay_pos_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))],...
                                                 [0;az_pos_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))]]
                            @(y) [y(:, 4:6), ax_norm_interpl(y(:, 1), y(:, 2), y(:, 3)),...
                                                 [0;ay_norm_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))],...
                                                 [0;az_norm_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))]]
                            @(y) [y(:, 4:6), ax_pos_H_interpl(y(:, 1), y(:, 3), -y(:, 2)),...
                                                 [0;-az_pos_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))],...
                                                 [0;ay_pos_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))]]
                            @(y) [y(:, 4:6), ax_norm_H_interpl(y(:, 1), y(:, 3), -y(:, 2)),...
                                                 [0;-az_norm_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))],...
                                                 [0;ay_norm_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))]]
                            @(y) [y(:, 4:6), ax_neg_interpl(y(:, 1), y(:, 2), y(:, 3)),...
                                                 [0;ay_neg_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))],...
                                                 [0;az_neg_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))]]};

                for i = 1:1: (length(obj.M_time_vec) - 2) %since free propagation is done with euler
                    if size(obj.xyzVxyz,1) == 0
                        error("All the molecules are lost after the %d switching!",i-1)
                    end
%                     fprintf("%d/%d\t%d\n",size(obj.xyzVxyz,1), obj.num_particles, i);
                    for t = obj.M_time_vec(i):dt: obj.M_time_vec(i+1)
                        obj.xyzVxyz = obj.xyzVxyz + dxyzVxyz{mod(i, 8)+8*(~mod(i, 8))}(obj.xyzVxyz) * dt;
                    end
                    [obj.xyzVxyz,obj.ind_particles]= removeHitParticles(obj.xyzVxyz, obj.ind_particles );
                end
                fprintf("%d out of %d particles left\n",size(obj.xyzVxyz,1), obj.num_particles);
            else    % normal mode
                dxyzVxyz = {@(y) [y(:, 4:6), ax_norm_interpl(y(:, 1), y(:, 2), y(:, 3)),...
                                                 [0;ay_norm_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))],...
                                                 [0;az_norm_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))]]                   
                            @(y) [y(:, 4:6), ax_norm_H_interpl(y(:, 1), y(:, 3), -y(:, 2)),...
                                                 [0;-az_norm_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))],...
                                                 [0;ay_norm_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))]]};
                
                for i = 1: 1: (length(obj.M_time_vec) - 2)
                    if size(obj.xyzVxyz,1) == 0
                        error("All the molecules are lost after the %d switching!", i-1)
                    end
%                     fprintf("%d/%d\t%d\n",size(obj.xyzVxyz,1), obj.num_particles,i);
                    for t = obj.M_time_vec(i):dt: obj.M_time_vec(i+1)
                        obj.xyzVxyz = obj.xyzVxyz + dxyzVxyz{mod(i, 2)+2*(~mod(i, 2))}(obj.xyzVxyz) * dt;
                    end
                    [obj.xyzVxyz,obj.ind_particles]= removeHitParticles(obj.xyzVxyz, obj.ind_particles );
                end
                fprintf("%d out of %d particles left\n",size(obj.xyzVxyz,1), obj.num_particles);
            end
            obj.Snapshot('end decelerator',t(end), obj.xyzVxyz, obj.ind_particles)
            obj.TOF()
            if obj.verbose         
            obj.Snapshotplot()
          
            obj.xyzVxyz(:,1:3)=obj.xyzVxyz(:,1:3)+obj.xyzVxyz(:,4:6)*(obj.M_time_vec(end)-obj.M_time_vec(end-1));                                                            
            figure();
            scatter(obj.xyzVxyz(:,1)*10^3, obj.xyzVxyz(:,4));  %scatter(x,y) creates a scatter plot with circular markers at the locations specified by the vectors x and y.
            xlabel('x-pos (mm)'); ylabel('v_x (m/s)');
            
            figure;
            obj.arrival_time = t(end) - (obj.xyzVxyz(:,1) - obj.params.PHYS_length_dec - obj.params.PHYS_exit_to_detection)./(obj.xyzVxyz(:,4)) + obj.params.FLY_incoupling_time;
            histogram(obj.arrival_time, 100);xlabel('time(s)')

            figure;
            scatter(obj.xyzVxyz(:,1), obj.xyzVxyz(:,4));%this compares and plots the phase spaces of snapshot (maybe make function sow that it could take any number of netries not jsut 2 as now)
            end

%           free propagation newton since velocity stays same thus
%           only update xyz by using x_new= x_old + v_x*dt for x,y,zphase

            
        end
        
        
%% Propagate particles, intergate using time intervals in M_time_vector, propagating all the molecules together
        function propagateParticles_euler_backup(obj)
            
            obj.createParticles();
            
%             a function that removes the lost molecules
            function [xyzVxyz,ind_particles] = removeHitParticles(xyzVxyz,ind_particles)
                hit_indices = xyzVxyz(:,1) < obj.params.PHYS_length_dec & (abs(xyzVxyz(:,2)) > obj.params.PHYS_seperation_pins/2 | abs(xyzVxyz(:,3)) > obj.params.PHYS_seperation_pins/2);
                xyzVxyz = xyzVxyz(~hit_indices,:);
                ind_particles = ind_particles(~hit_indices);
                hit_indices = (abs(xyzVxyz(:,1) - xyzVxyz(1,1)) < 2* 5.5e-3); % separated the cut of the fish in x into two lines in order to make indexing work and also to be sure its correct
                xyzVxyz = xyzVxyz(hit_indices, :);
                ind_particles = ind_particles(hit_indices);
            end
            
            % first propagate to entrance of decelerator
            obj.xyzVxyz = obj.xyzVxyz_0;
            obj.ind_particles = transpose(1:size(obj.xyzVxyz,1));
            obj.xyzVxyz(:,1:3) = obj.xyzVxyz(:,1:3) + obj.xyzVxyz(:,4:6) * obj.params.FLY_incoupling_time; % ekin at start of deacc.
            [obj.xyzVxyz,obj.ind_particles]= removeHitParticles(obj.xyzVxyz, obj.ind_particles );% select those that can enter dec
%             hit_indices = obj.xyzVxyz(:,1) < obj.params.PHYS_length_dec & (abs(obj.xyzVxyz(:,2)) > obj.params.PHYS_seperation_pins/2 | abs(obj.xyzVxyz(:,3)) > obj.params.PHYS_seperation_pins/2);
%             obj.xyzVxyz = obj.xyzVxyz(~hit_indices,:);
%             obj.xyzVxyz = obj.xyzVxyz((abs(obj.xyzVxyz(:,1) - obj.xyzVxyz(1,1)) < 2* 5.5e-3), :);

            % propagate inside dec
            
            dt = 4e-8;
%             fprintf("num/total switching\n");
            obj.Snapshot("start decelerator",obj.M_time_vec(1),obj.xyzVxyz, obj.ind_particles)
            if obj.params.FLY_focusing_mode_bool % focusing mode
                dxyzVxyz = {@(y) [y(:, 4:6), obj.ax_norm_interpl(y(:, 1), y(:, 2), y(:, 3)),...
                                                 [0;obj.ay_norm_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))],...
                                                 [0;obj.az_norm_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))]]                   
                            @(y) [y(:, 4:6), obj.ax_neg_H_interpl(y(:, 1), y(:, 3), -y(:, 2)),...
                                                 [0;-obj.az_neg_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))],...
                                                 [0;obj.ay_neg_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))]]
                            @(y) [y(:, 4:6), obj.ax_norm_H_interpl(y(:, 1), y(:, 3), -y(:, 2)),...
                                                 [0;-obj.az_norm_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))],...
                                                 [0;obj.ay_norm_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))]]
                            @(y) [y(:, 4:6), obj.ax_pos_interpl(y(:, 1), y(:, 2), y(:, 3)),...
                                                 [0;obj.ay_pos_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))],...
                                                 [0;obj.az_pos_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))]]
                            @(y) [y(:, 4:6), obj.ax_norm_interpl(y(:, 1), y(:, 2), y(:, 3)),...
                                                 [0;obj.ay_norm_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))],...
                                                 [0;obj.az_norm_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))]]
                            @(y) [y(:, 4:6), obj.ax_pos_H_interpl(y(:, 1), y(:, 3), -y(:, 2)),...
                                                 [0;-obj.az_pos_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))],...
                                                 [0;obj.ay_pos_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))]]
                            @(y) [y(:, 4:6), obj.ax_norm_H_interpl(y(:, 1), y(:, 3), -y(:, 2)),...
                                                 [0;-obj.az_norm_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))],...
                                                 [0;obj.ay_norm_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))]]
                            @(y) [y(:, 4:6), obj.ax_neg_interpl(y(:, 1), y(:, 2), y(:, 3)),...
                                                 [0;obj.ay_neg_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))],...
                                                 [0;obj.az_neg_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))]]};

                for i = 1:1: (length(obj.M_time_vec) - 2) %since free propagation is done with euler
                    if size(obj.xyzVxyz,1) == 0
                        error("All the molecules are lost after the %d switching!",i-1)
                    end
%                     fprintf("%d/%d\t%d\n",size(obj.xyzVxyz,1), obj.num_particles, i);
                    for t = obj.M_time_vec(i):dt: obj.M_time_vec(i+1)
                        obj.xyzVxyz = obj.xyzVxyz + dxyzVxyz{mod(i, 8)+8*(~mod(i, 8))}(obj.xyzVxyz) * dt;
                    end
                    [obj.xyzVxyz,obj.ind_particles]= removeHitParticles(obj.xyzVxyz, obj.ind_particles );
                end
                fprintf("%d out of %d particles left\n",size(obj.xyzVxyz,1), obj.num_particles);
            else    % normal mode
                dxyzVxyz = {@(y) [y(:, 4:6), obj.ax_norm_interpl(y(:, 1), y(:, 2), y(:, 3)),...
                                                 [0;obj.ay_norm_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))],...
                                                 [0;obj.az_norm_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))]]                   
                            @(y) [y(:, 4:6), obj.ax_norm_H_interpl(y(:, 1), y(:, 3), -y(:, 2)),...
                                                 [0;-obj.az_norm_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))],...
                                                 [0;obj.ay_norm_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))]]};
                
                for i = 1: 1: (length(obj.M_time_vec) - 2)
                    if size(obj.xyzVxyz,1) == 0
                        error("All the molecules are lost after the %d switching!", i-1)
                    end
%                     fprintf("%d/%d\t%d\n",size(obj.xyzVxyz,1), obj.num_particles,i);
                    for t = obj.M_time_vec(i):dt: obj.M_time_vec(i+1)
                        obj.xyzVxyz = obj.xyzVxyz + dxyzVxyz{mod(i, 2)+2*(~mod(i, 2))}(obj.xyzVxyz) * dt;
                    end
                    [obj.xyzVxyz,obj.ind_particles]= removeHitParticles(obj.xyzVxyz, obj.ind_particles );
                end
                fprintf("%d out of %d particles left\n",size(obj.xyzVxyz,1), obj.num_particles);
            end
            obj.Snapshot('end decelerator',t(end), obj.xyzVxyz, obj.ind_particles)
            obj.TOF()
            if obj.verbose
              
            obj.Snapshotplot()
          
            obj.xyzVxyz(:,1:3)=obj.xyzVxyz(:,1:3)+obj.xyzVxyz(:,4:6)*(obj.M_time_vec(end)-obj.M_time_vec(end-1));                                                            
            figure();
            scatter(obj.xyzVxyz(:,1)*10^3, obj.xyzVxyz(:,4));  %scatter(x,y) creates a scatter plot with circular markers at the locations specified by the vectors x and y.
            xlabel('x-pos (mm)'); ylabel('v_x (m/s)');
            
            figure;
            obj.arrival_time = t(end) - (obj.xyzVxyz(:,1) - obj.params.PHYS_length_dec - obj.params.PHYS_exit_to_detection)./(obj.xyzVxyz(:,4)) + obj.params.FLY_incoupling_time;
            histogram(obj.arrival_time, 100);xlabel('time(s)')

            figure;
            scatter(obj.xyzVxyz(:,1), obj.xyzVxyz(:,4));%this compares and plots the phase spaces of snapshot (maybe make function sow that it could take any number of netries not jsut 2 as now)
            end

%           free propagation newton since velocity stays same thus
%           only update xyz by using x_new= x_old + v_x*dt for x,y,zphase

            
        end
        
         %% Propagate particles, intergate using time intervals in M_time_vector, propagating all the molecules together
         function propagateParticles_ode45(obj,in_l)
             
             obj.createParticles();
             
             % make local variables of interpolations improves running speed
             %a  bit

             ax_norm_interpl= in_l.ax_norm_interpl;
             ay_norm_interpl= in_l.ay_norm_interpl;
             az_norm_interpl= in_l.az_norm_interpl;

             ax_norm_H_interpl = in_l.ax_norm_H_interpl;
             ay_norm_H_interpl = in_l.ay_norm_H_interpl;
             az_norm_H_interpl = in_l.az_norm_H_interpl;

             ax_neg_interpl = in_l.ax_neg_interpl;
             ay_neg_interpl = in_l.ay_neg_interpl;
             az_neg_interpl = in_l.az_neg_interpl;

             ax_neg_H_interpl = in_l.ax_neg_H_interpl;
             ay_neg_H_interpl = in_l.ay_neg_H_interpl;
             az_neg_H_interpl = in_l.az_neg_H_interpl;

             ax_pos_interpl = in_l.ax_pos_interpl;
             ay_pos_interpl = in_l.ay_pos_interpl;
             az_pos_interpl = in_l.az_pos_interpl;
             
             ax_pos_H_interpl = in_l.ax_pos_H_interpl;
             ay_pos_H_interpl = in_l.ay_pos_H_interpl;
             az_pos_H_interpl = in_l.az_pos_H_interpl;

             M_time_vec_l = in_l.M_time_vec;



            % THIS FUNCTION IS DEFINED TWICE!!! DOES IT MAKES ANY SENSE?
            % a function that remove the lost molecules
        function [xyzVxyz,ind_particles] = removeHitParticles(xyzVxyz,ind_particles)
                hit_indices = xyzVxyz(:,1) < in_l.params.PHYS_length_dec & (abs(xyzVxyz(:,2)) > in_l.params.PHYS_seperation_pins/2 | abs(xyzVxyz(:,3)) > in_l.params.PHYS_seperation_pins/2);
                xyzVxyz = xyzVxyz(~hit_indices,:);
                ind_particles = ind_particles(~hit_indices);
                hit_indices = (abs(xyzVxyz(:,1) - xyzVxyz(1,1)) < 2* 5.5e-3); % separated the cut of the fish in x into two lines in order to make indexing work and also to be sure its correct
                xyzVxyz = xyzVxyz(hit_indices, :);
                ind_particles = ind_particles(hit_indices);
        end        
            
            % first propagate to entrance of decelerator
            obj.xyzVxyz = obj.xyzVxyz_0;
            obj.ind_particles = transpose(1:size(obj.xyzVxyz,1));
            obj.xyzVxyz(:,1:3) = obj.xyzVxyz(:,1:3) + obj.xyzVxyz(:,4:6) * in_l.params.FLY_incoupling_time;
            [obj.xyzVxyz,obj.ind_particles]= removeHitParticles(obj.xyzVxyz, obj.ind_particles );% select those that can enter dec
           

            % propagate inside dec
%             fprintf("num/total switching\n");
            opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
            times = [];   % used for commented plots below which tell you plot you the time steps of the integration
            num_par = [];
            obj.output={};
            obj.Snapshot("start decelerator", M_time_vec_l(1),obj.xyzVxyz, obj.ind_particles)
            if in_l.params.FLY_focusing_mode_bool % focusing mode
                dydt_array = {@(t,y) dydtNormVerticalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) dydtNegHorizontalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) dydtNormHorizontalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) dydtPosVerticalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) dydtNormVerticalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) dydtPosHorizontalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) dydtNormHorizontalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) dydtNegVerticalOn(t,y,size(obj.xyzVxyz,1))};
                tic
                for i = 1:1: (length( M_time_vec_l) - 2)
                    if size(obj.xyzVxyz,1) == 0
                        error("All the molecules are lost after the %d switching!",i-1)
                    end
                    [t, y] = ode45(dydt_array{mod(i, 8)+8*(~mod(i, 8))}, [ M_time_vec_l(i), ( M_time_vec_l(i) +  M_time_vec_l(i+1))/2,  M_time_vec_l(i+1)], reshape(obj.xyzVxyz, [], 1), opts);
                    obj.xyzVxyz = reshape(y(end,:), [], 6);
                    [obj.xyzVxyz,obj.ind_particles]= removeHitParticles(obj.xyzVxyz,obj.ind_particles);
                end
                fprintf("%d out of %d particles left\n",size(obj.xyzVxyz,1), obj.num_particles);
            else    % normal mode
                dydt_array = {@(t,y) obj.dydtNormVerticalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtNormHorizontalOn(t,y,size(obj.xyzVxyz,1))};
                for i = 1: 1: (length( M_time_vec_l) - 2)
                    if size(obj.xyzVxyz,1) == 0
                        error("All the molecules are lost after the %d switching!", i-1)
                    end
%                     fprintf("%d/%d\t%d\n",size(obj.xyzVxyz,1), obj.num_particles,i);
                    [t, y] = ode45(dydt_array{mod(i, 2)+2*(~mod(i, 2))}, [ M_time_vec_l(i), ( M_time_vec_l(i) +  M_time_vec_l(i+1))/2,  M_time_vec_l(i+1)], reshape(obj.xyzVxyz, [], 1), opts);             
%                   times = [times toc];
%                   num_par = [num_par, size(obj.xyzVxyz,1)]; 
                    obj.xyzVxyz = reshape(y(end,:), [], 6);
                    [obj.xyzVxyz,obj.ind_particles]= removeHitParticles(obj.xyzVxyz,obj.ind_particles);
                end
                fprintf("%d out of %d particles left\n",size(obj.xyzVxyz,1), obj.num_particles);
                toc
            end
            obj.Snapshot('end decelerator',t(end), obj.xyzVxyz, obj.ind_particles)

            obj.TOF()  
            if in_l.verbose

%                 figure()
%                 subplot(2, 2, 1)
%                 title('Time steps')
%                 plot(times,'.')
%                 xlabel('integration step'); ylabel('integration time step')
%                 subplot(2, 2, 2)
%                 title('number particles')
%                 plot(num_par,'.')
%                 xlabel('integration step'); ylabel('# particles')
%                 subplot(2, 2, 3)
%                 plot(num_par, times,'.')
%                 xlabel('# partciles'); ylabel('integration time of step')
% 

            obj.Snapshotplot() %this compares and plots the phase spaces of snapshot (maybe make function sow that it could take any number of netries not jsut 2 as now)
            obj.xyzVxyz(:,1:3)=obj.xyzVxyz(:,1:3)+obj.xyzVxyz(:,4:6)*( M_time_vec_l(end)- M_time_vec_l(end-1));                                                            
            figure();
            scatter(obj.xyzVxyz(:,1)*10^3, obj.xyzVxyz(:,4));  %scatter(x,y) creates a scatter plot with circular markers at the locations specified by the vectors x and y.
            xlabel('x-pos (mm)'); ylabel('v_x (m/s)');
            
            figure;
            obj.arrival_time = t(end) - (obj.xyzVxyz(:,1) - in_l.params.PHYS_length_dec - in_l.params.PHYS_exit_to_detection)./(obj.xyzVxyz(:,4)) + in_l.params.FLY_incoupling_time;
            histogram(obj.arrival_time, 100);xlabel('time(s)')
            end

             function dydt = dydtNormVerticalOn(t,y,n)
                 dydt = [y(3*n+1:end); ax_norm_interpl(y(1:n), y(n+1:2*n), y(2*n+1:3*n)); [0;ay_norm_interpl(y(2:n), y(n+2:2*n), y(2*n+2:3*n))]; [0; az_norm_interpl(y(2:n), y(n+2:2*n), y(2*n+2:3*n))]];
             end

             function dydt = dydtNormHorizontalOn(t,y,n)
                 dydt = [y(3*n+1:end); ax_norm_H_interpl(y(1:n), y(2*n+1:3*n), -y(n+1:2*n)); [0;-az_norm_H_interpl(y(2:n), y(2*n+2:3*n), -y(n+2:2*n))]; [0;ay_norm_H_interpl(y(2:n), y(2*n+2:3*n), -y(n+2:2*n))]];
             end

             function dydt = dydtNegVerticalOn(t,y,n)
                 dydt = [y(3*n+1:end); ax_neg_interpl(y(1:n), y(n+1:2*n), y(2*n+1:3*n)); [0; ay_neg_interpl(y(2:n), y(n+2:2*n), y(2*n+2:3*n))]; [0; az_neg_interpl(y(2:n), y(n+2:2*n), y(2*n+2:3*n))]];
             end

             function dydt = dydtNegHorizontalOn(t,y,n)
                 dydt = [y(3*n+1:end); ax_neg_H_interpl(y(1:n), y(2*n+1:3*n), -y(n+1:2*n)); [0;-az_neg_H_interpl(y(2:n), y(2*n+2:3*n), -y(n+2:2*n))]; [0;ay_neg_H_interpl(y(2:n), y(2*n+2:3*n), -y(n+2:2*n))]];
             end

             function dydt = dydtPosVerticalOn(t,y,n)
                 dydt = [y(3*n+1:end); ax_pos_interpl(y(1:n), y(n+1:2*n), y(2*n+1:3*n)); [0; ay_pos_interpl(y(2:n), y(n+2:2*n), y(2*n+2:3*n))]; [0; az_pos_interpl(y(2:n), y(n+2:2*n), y(2*n+2:3*n))]];
             end
             function dydt = dydtPosHorizontalOn(t,y,n)
                 dydt = [y(3*n+1:end); ax_pos_H_interpl(y(1:n), y(2*n+1:3*n), -y(n+1:2*n)); [0;-az_pos_H_interpl(y(2:n), y(2*n+2:3*n), -y(n+2:2*n))]; [0; ay_pos_H_interpl(y(2:n), y(2*n+2:3*n), -y(n+2:2*n))]];
             end                                             
            
        end
        
        %% Obtain the acceleration by interpolation
        function dydt = dydtNormVerticalOn(obj,t,y,n)
           dydt = [y(3*n+1:end); obj.in.ax_norm_interpl(y(1:n), y(n+1:2*n), y(2*n+1:3*n)); [0;obj.in.ay_norm_interpl(y(2:n), y(n+2:2*n), y(2*n+2:3*n))]; [0; obj.in.az_norm_interpl(y(2:n), y(n+2:2*n), y(2*n+2:3*n))]]; 
        end
        
        function dydt = dydtNormHorizontalOn(obj,t,y,n)
           dydt = [y(3*n+1:end); obj.in.ax_norm_H_interpl(y(1:n), y(2*n+1:3*n), -y(n+1:2*n)); [0;-obj.in.az_norm_H_interpl(y(2:n), y(2*n+2:3*n), -y(n+2:2*n))]; [0;obj.in.ay_norm_H_interpl(y(2:n), y(2*n+2:3*n), -y(n+2:2*n))]]; 
        end
        
        function dydt = dydtNegVerticalOn(obj,t,y,n)
            dydt = [y(3*n+1:end); obj.in.ax_neg_interpl(y(1:n), y(n+1:2*n), y(2*n+1:3*n)); [0; obj.in.ay_neg_interpl(y(2:n), y(n+2:2*n), y(2*n+2:3*n))]; [0; obj.in.az_neg_interpl(y(2:n), y(n+2:2*n), y(2*n+2:3*n))]];
        end
        
        function dydt = dydtNegHorizontalOn(obj,t,y,n)
            dydt = [y(3*n+1:end); obj.in.ax_neg_H_interpl(y(1:n), y(2*n+1:3*n), -y(n+1:2*n)); [0;-obj.in.az_neg_H_interpl(y(2:n), y(2*n+2:3*n), -y(n+2:2*n))]; [0;obj.in.ay_neg_H_interpl(y(2:n), y(2*n+2:3*n), -y(n+2:2*n))]]; 
        end

        function dydt = dydtPosVerticalOn(obj,t,y,n)
            dydt = [y(3*n+1:end); obj.in.ax_pos_interpl(y(1:n), y(n+1:2*n), y(2*n+1:3*n)); [0; obj.in.ay_pos_interpl(y(2:n), y(n+2:2*n), y(2*n+2:3*n))]; [0; obj.in.az_pos_interpl(y(2:n), y(n+2:2*n), y(2*n+2:3*n))]];
        end 
        function dydt = dydtPosHorizontalOn(obj,t,y,n)
            dydt = [y(3*n+1:end); obj.in.ax_pos_H_interpl(y(1:n), y(2*n+1:3*n), -y(n+1:2*n)); [0;-obj.in.az_pos_H_interpl(y(2:n), y(2*n+2:3*n), -y(n+2:2*n))]; [0;obj.in.ay_pos_H_interpl(y(2:n), y(2*n+2:3*n), -y(n+2:2*n))]]; 
        end

        
        function propagateParticlesAndSaveTrajectories(obj)
            
            % a function that remove the lost molecules
            function xyzVxyz = removeHitParticles(y)
                xyzVxyz = reshape(y(end,:), [], 6);
                remaining_indices = ((xyzVxyz(:,1) < obj.in.params.PHYS_length_dec & abs(xyzVxyz(:,2)) < obj.in.params.PHYS_seperation_pins/2 & abs(xyzVxyz(:,3)) < obj.in.params.PHYS_seperation_pins/2)...
                                    | xyzVxyz(:,1) >= obj.in.params.PHYS_length_dec) & abs(xyzVxyz(:,1) - xyzVxyz(1,1)) < 2*obj.in.params.PHYS_distance_stages;
                xyzVxyz = xyzVxyz(remaining_indices,:);
%                 size(remaining_indices)
%                 size(y)
                if obj.num_trajectories_saved > 0
%                     y = permute(reshape(y', 6, [], size(y,1)),[2,1,3]);
                      y = reshape(y', [], 6, size(y,1));
%                     size(y)
                    y = y(remaining_indices,:,:);
%                     size(y)
                    obj.traj_xyzVxyz = obj.traj_xyzVxyz(remaining_indices,:,:);
                    obj.traj_xyzVxyz = cat(3, obj.traj_xyzVxyz, y);
                end
            end

            obj.traj_xyzVxyz = zeros(obj.num_particles, 6, 1);
            if obj.num_trajectories_saved > 0
                obj.traj_time = 0.0;
                obj.traj_xyzVxyz(:,:,1) = obj.xyzVxyz_0;% xlabel('time ( /mu s)'); ylabel('detecetd moelcuels'); legend('euler','ode45','synch. mol.')
            else
                fprintf("obj.num_trajectories == 0, no trajectories will be saved!");
            end
            
            % first propagate to entrance of decelerator
            obj.xyzVxyz = obj.xyzVxyz_0;
            obj.xyzVxyz(:,1:3) = obj.xyzVxyz(:,1:3) + obj.xyzVxyz(:,4:6) * obj.in.params.FLY_incoupling_time;
%             obj.xyzVxyz= removeHitParticles(obj.xyzVxyz);% select those that can enter dec

            
            % propagate inside dec
            fprintf("num/total switching\n");
%             line_temp = 0;
            
            opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
            if obj.in.params.FLY_focusing_mode_bool % focusing mode
                dydt_array = {@(t,y) obj.dydtNormVerticalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtNegHorizontalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtNormHorizontalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtPosVerticalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtNormVerticalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtPosHorizontalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtNormHorizontalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtNegVerticalOn(t,y,size(obj.xyzVxyz,1))};
                 
                for i = 1:1:(length(obj.in.M_time_vec) - 1)
                    if size(obj.xyzVxyz,1) == 0
                        error("All the molecules are lost after the %d switching!",i-1)
                    end
%                     fprintf("%d/%d\t%d\n",size(obj.xyzVxyz,1), obj.num_particles, i);
%                     [t, y] = ode45(dydt_array{mod(i, 8)+8*(~mod(i, 8))}, obj.in.M_time_vec(i):1e-6:obj.in.M_time_vec(i+1), reshape(obj.xyzVxyz, [], 1), opts);
                    [t, y] = ode45(dydt_array{mod(i, 8)+8*(~mod(i, 8))}, [obj.in.M_time_vec(i), (obj.in.M_time_vec(i) + obj.in.M_time_vec(i+1))/2, obj.in.M_time_vec(i+1)], reshape(obj.xyzVxyz, [], 1), opts);

                    obj.traj_time = [obj.traj_time; t + obj.in.params.FLY_incoupling_time];
                    obj.xyzVxyz= removeHitParticles(y);
                end
                
            else    % normal mode
                dydt_array = {@(t,y) obj.dydtNormVerticalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtNormHorizontalOn(t,y,size(obj.xyzVxyz,1))};
                for i = 1: 1: (length(obj.in.M_time_vec) - 1)
                    if size(obj.xyzVxyz,1) == 0
                        error("All the molecules are lost after the %d switching!", i-1)
                    end
                    fprintf("%d/%d\t%d\n",size(obj.xyzVxyz,1), obj.num_particles,i);
                    [t, y] = ode45(dydt_array{mod(i, 2)+2*(~mod(i, 2))}, [obj.in.M_time_vec(i), obj.in.M_time_vec(i+1)], reshape(obj.xyzVxyz, [], 1), opts);
%                     [t, y] = ode45(dydt_array{mod(i, 2)+2*(~mod(i, 2))}, [obj.in.M_time_vec(i), (obj.in.M_time_vec(i) + obj.in.M_time_vec(i+1))/2, obj.in.M_time_vec(i+1)], reshape(obj.xyzVxyz, [], 1), opts);
                    
                    obj.traj_time = [obj.traj_time; t];
                    obj.xyzVxyz= removeHitParticles(y);
                end
            end

            obj.arrival_time = t(end) - (obj.xyzVxyz(:,1) - obj.in.params.PHYS_length_dec - obj.in.params.PHYS_exit_to_detection)./(obj.xyzVxyz(:,4)) + obj.in.params.FLY_incoupling_time;
            
        end
        
        function plotTrajectories(obj)
            if size(obj.traj_xyzVxyz, 1) > obj.num_trajectories_saved
                obj.traj_xyzVxyz = obj.traj_xyzVxyz(1:obj.num_trajectories_saved,:,:);
            end
            
            subplot(2,2,1);
            histogram(obj.arrival_time, 100);
            xlabel('arrival time(s)'); ylabel('arb.u.'); title('TOF at detection')
            
            subplot(2,2,2);
            scatter(obj.xyzVxyz(:,1), obj.xyzVxyz(:,4));
            xlabel('longitudinal position (m)'); ylabel('longitudinal velocity(m/s)'); title('Phase space at detection')
            
            subplot(2,2,3);
            for i=1:size(obj.traj_xyzVxyz,1)
                plot(squeeze(obj.traj_xyzVxyz(i,1,:)), squeeze(obj.traj_xyzVxyz(i,4,:)))
                hold on;
            end
            xlabel('longitudinal position(m)'); ylabel('longitudinal velocity(m/s)'); title('longitudinal velocity vs position')
            
            subplot(2,2,4);
            for i=1:size(obj.traj_xyzVxyz,1)
                plot(obj.traj_time, squeeze(obj.traj_xyzVxyz(i,4,:)))
                hold on;
            end
            xlabel('time(s)'); ylabel('longitudinal velocity(m/s)'); title('longitudial velocity evolution')
            
            
            figure;
            subplot(2,2,1)
            for i=1:size(obj.traj_xyzVxyz,1)
                plot(squeeze(obj.traj_xyzVxyz(i,1,:)), squeeze(obj.traj_xyzVxyz(i,2,:)))
                hold on;
            end
            xlabel('longitudinal position (m)'); ylabel('transverse position(m)'); title('2d trajectories (z vs x)')
            
            subplot(2,2,2)
            for i=1:size(obj.traj_xyzVxyz,1)
                plot(squeeze(obj.traj_xyzVxyz(i,1,:)), squeeze(obj.traj_xyzVxyz(i,3,:)))
                hold on;
            end
            xlabel('longitudinal position (m)'); ylabel('transverse position(m)'); title('2d trajectories (y vs x)')
            
            subplot(2,2,3)
            for i=1:size(obj.traj_xyzVxyz,1)
                plot3(squeeze(obj.traj_xyzVxyz(i,1,:)), squeeze(obj.traj_xyzVxyz(i,2,:)), squeeze(obj.traj_xyzVxyz(i,3,:)))
                hold on;
            end
            xlabel('longitudinal position (m)'); ylabel('position z(m)'); ylabel('position y(m)'); title('3d trajectories')
             
            subplot(2,2,4);
            for i=1:size(obj.traj_xyzVxyz,1)
                plot(obj.traj_time, squeeze(obj.traj_xyzVxyz(i,4,:)))
                hold on;
            end
            xlabel('time(s)'); ylabel('longitudinal velocity(m/s)'); title('longitudial velocity evolution')
            
        end
        
        %% Plot time-of-flight signal
        function plotTOF(obj)
            min_time = 0; max_time = 10e-3;
            num_bins = 6000;
            binsize = (max_time - min_time)/num_bins;
            tof_profile = zeros(num_bins,1);
            
            
            obj.arrival_time = obj.arrival_time(abs(obj.xyzVxyz(:,2)) < obj.params.FLY_detection_laser_diameter/2 & abs(obj.xyzVxyz(:,3)) < 2e-3);
            bin_begin_each_particle = int16((obj.arrival_time - obj.params.FLY_detection_laser_diameter/2/obj.xyzVxyz(:,4) - min_time)/binsize + 0.5);
            bin_end_each_particle =int16((obj.arrival_time + obj.params.FLY_detection_laser_diameter/2/obj.xyzVxyz(:,4) - min_time)/binsize + 0.5);
            
            time = min_time:binsize:max_time;
            for i = 1: length(bin_begin_each_particle)
                tof_profile(bin_begin_each_particle(i):bin_end_each_particle(i)) = tof_profile(bin_begin_each_particle(i):bin_end_each_particle(i)) + 1;
            end
            
%             disp(size(time));
%             disp(size(tof_profile));
            figure;
            plot(time(1:end-1), tof_profile);
            xlabel('arrival time (s)')
            ylabel('signal (arb. u)')
        end

         %% Plot TOF with laser volume
        % redone with laser volume and snapshot of xyzVxyz when snych. molecule is at detection point
        function TOF(obj) 
            t_steps_TOF = 1e-6; %time srteps for gfree propagation one smaller than t_profile
            % t_profile = (obj.M_time_vec(end-1) + t_steps_TOF ):t_steps_TOF:(obj.M_time_vec(end)+1e-3); % start 1us after, go on for
            t_profile = t_steps_TOF : t_steps_TOF : (obj.M_time_vec(end) - obj.M_time_vec(end-1) + 1e-3); % start 1us after, go on for 1 ms
            
            x_laser = obj.params.PHYS_length_dec + obj.params.PHYS_exit_to_detection;  % x coordinate laser center
            h_laser = 2e-3; %half of height laser volume to make check since height is 4 mm but goes from y=-2mm to y=2mm
            r_laser = obj.params.FLY_detection_laser_diameter/2; % radius of laser we model as cylinder
            
            obj.TOF_save = zeros(length(t_profile), 2); % matrix wher number of molecules in laser volume will be saved with the corresponding time
            obj.TOF_save(:, 1) = t_profile + obj.M_time_vec(end-1);
            output_xyzVxyz = {};

            % propagate till synch molecules
            xyzVxyz_TOF = obj.xyzVxyz; % local copy not to modify class variable
            xyzVxyz_TOF(:, 1:3) = obj.xyzVxyz(:, 1:3) + obj.xyzVxyz(:,4:6)*(obj.M_time_vec(end) - obj.M_time_vec(end-1)); 
            in_volume =  abs( xyzVxyz_TOF(:,3) ) <= h_laser & sqrt((xyzVxyz_TOF(:,1)-x_laser).^2 + (xyzVxyz_TOF(:,2)).^2) <= r_laser;
            xyzVxyz_TOF = xyzVxyz_TOF(in_volume, :);
            obj.Snapshot('synch. molecule detection', obj.M_time_vec(end), xyzVxyz_TOF,[]);


            for i=1:length(t_profile)

                xyzVxyz_TOF = obj.xyzVxyz; % local copy not to modify class variable
 %               xyzVxyz_TOF(:, 1:3) = xyzVxyz_TOF(:, 1:3) + xyzVxyz_TOF(:,4:6)*t_steps_TOF; % let particles propagate 1 mu s after decc           
                xyzVxyz_TOF(:, 1:3) = obj.xyzVxyz(:, 1:3) + obj.xyzVxyz(:,4:6)*t_profile(i); % let particles propagate 1 mu s after decc           
                in_volume =  abs( xyzVxyz_TOF(:,3) ) <= h_laser & sqrt((xyzVxyz_TOF(:,1)-x_laser).^2 + (xyzVxyz_TOF(:,2)).^2) <= r_laser;   %check if a particle is in y range of cylinder (array of ones and zeros) one if inside range
                xyzVxyz_TOF = xyzVxyz_TOF(in_volume, :);
                number_part_det = sum(in_volume);

                if number_part_det<0
                    fprintf('Ahiahiahi')
                end
                obj.TOF_save(i,2) = number_part_det; %save the number of detected particles at that time
                if ~isempty(xyzVxyz_TOF)     
                    output_xyzVxyz = [output_xyzVxyz, {xyzVxyz_TOF} ];
                end
                clearvars xyzVxyz_TOF
            end

            xyzVxyz_TOF = obj.xyzVxyz;         % initalize again local variable to propagate each particle by itself and see if it gets detected or not
            xyzVxyz_in_volume = zeros(1,length(obj.xyzVxyz(:,1)));  % save it here 1 if particle was detectted 0 if not

            for j = 1:length(obj.xyzVxyz(:,1))
                for i = 1:length(t_profile)
                    xyzVxyz_TOF(j, 1:3) = obj.xyzVxyz(j, 1:3) + obj.xyzVxyz(j,4:6)*t_profile(i); % let particles propagate 1 mu s after decc
                    if abs( xyzVxyz_TOF(j,3) ) <= h_laser & sqrt((xyzVxyz_TOF(j,1)-x_laser).^2 + (xyzVxyz_TOF(j,2)).^2) <= r_laser
                        xyzVxyz_in_volume(j)=1;
                        break
                    end
                end
            end

            obj.Snapshot('end decc. detected (safe)',obj.M_time_vec(end-1), obj.xyzVxyz(logical(xyzVxyz_in_volume),:), []) % here time is chosen at end decc. since pos. and vel. xyzVxyz correspond to this time
           
            compare=output_xyzVxyz{1,1}; % idea below is to compare velocities of all particles detected around the synch molecule se above if statement such that
                                         % we get a good estiamte for the
                                        % area of the peak
            for j =1: length(output_xyzVxyz)
                row_a= ismember(output_xyzVxyz{1,j}(:,4:end),compare(:,4:end),'rows'); %compare vel. of partciles to distinguish them compare all rows at same time
                compare= [compare;output_xyzVxyz{1,j}(~row_a,:)];  % add only the rows which were not yet present in compare to not count them twice
            end
            row_end = ismember(obj.xyzVxyz(:,4:end),compare(:,4:end),'rows');
            obj.Snapshot('end decc, detected', obj.M_time_vec(end-1),obj.xyzVxyz(row_end,:),[])
%% now there is fucntion we can call that does that for us plot_TOF_laser
%             if obj.verbose
%                 figure('Name','TOF')
%                 plot(obj.TOF_save(:,1)*1e6,obj.TOF_save(:,2))     
%                 xlabel('time (/mu s'); ylabel('# particles in laser volume');            
%                 ylim([0,max(obj.TOF_save(:,2)+2)]);            
%                 xlim([obj.TOF_save(1,1)*1e6,obj.TOF_save(end,1)*1e6]);
%             end
        end

        %fucntion taking a snapshot of time, xyzVxyz and the corresponding
        %indexes of the particles as well a flag that tells you the point
        function Snapshot(obj, flag, t, pos_vel, index_p)
            obj.output= [obj.output; {flag,t,pos_vel,index_p}];
        end

        function Snapshotplot(obj)
            [leng,~]=size(obj.output);
            output_xyzVxyz={}; 
            for j=1:leng   %Phase space plots of each time instance saved in output file
                 figure()
                 subplot(3,1,1)
                 title('Phase space x,v_x', obj.output{j,1})
                 hold on 
                 plot(obj.output{j,3}(:,1)*10^3,obj.output{j,3}(:,4)*10^3, '.')
                 plot(obj.output{j,3}(1,1)*10^3,obj.output{j,3}(1,4)*10^3,'.r')
                 hold off
                 xlabel('x_{pos.} (mm)'), ylabel('v_x (m/s)')
                 subplot(3,1,2)
                 title('Phase space y,v_y',obj.output{j,1})
                 hold on 
                 plot(obj.output{j,3}(:,2)*10^3,obj.output{j,3}(:,5)*10^3, '.')
                 plot(obj.output{j,3}(1,2)*10^3,obj.output{j,3}(1,5)*10^3,'.r')
                 hold off
                 xlabel('y_{pos.} (mm)'), ylabel('v_y (m/s)')
                 subplot(3,1,3)
                 title('Phase space z, v_z',obj.output{j,1} )
                 hold on 
                 plot(obj.output{j,3}(:,3)*10^3,obj.output{j,3}(:,6)*10^3, '.')
                 plot(obj.output{j,3}(1,3)*10^3,obj.output{j,3}(1,6)*10^3,'.r')
                 hold off
                 xlabel('z_{pos.} (mm)'), ylabel('v_z (m/s)')

%                  figure() %3D vector plots not very usefull but noce to look at
%                  quiver3(obj.output{j,3}(:,1)*10^3,obj.output{j,3}(:,3)*10^3,obj.output{j,3}(:,2)*10^3,obj.output{j,3}(:,4),obj.output{j,3}(:,6),obj.output{j,3}(:,5))
%                  xlabel('x (m)'); ylabel('z (m)'); zlabel('y(m)')

            end

            figure() %compare end and start by histograms plotting rho= |x-x_synch|*|v_x-v_x(synch)| for the x,vx y,vy, z,vz pahse spaces
            subplot(3,1,1)
            hold on
            histogram(abs(obj.output{1,3}(:,1)-obj.output{1,3}(1,1)).*abs(obj.output{1,3}(:,4)-obj.output{1,3}(1,4)),100)
            histogram(abs(obj.output{2,3}(:,1)-obj.output{2,3}(1,1)).*abs(obj.output{2,3}(:,4)-obj.output{2,3}(1,4)),100)
            legend(obj.output{1,1}, obj.output{2,1})
            hold off
            xlabel("/rho")
            subplot(3,1,2)
            hold on
            histogram(abs(obj.output{1,3}(:,2)-obj.output{1,3}(1,2)).*abs(obj.output{1,3}(:,5)-obj.output{1,3}(1,5)),100)
            histogram(abs(obj.output{2,3}(:,2)-obj.output{2,3}(1,2)).*abs(obj.output{2,3}(:,5)-obj.output{2,3}(1,5)),100)
            legend(obj.output{1,1}, obj.output{2,1})
            hold off
            xlabel("/rho")
            subplot(3,1,3)
            hold on
            histogram(abs(obj.output{1,3}(:,3)-obj.output{1,3}(1,3)).*abs(obj.output{1,3}(:,6)-obj.output{1,3}(1,6)),100)
            histogram(abs(obj.output{2,3}(:,3)-obj.output{2,3}(1,3)).*abs(obj.output{2,3}(:,6)-obj.output{2,3}(1,6)),100)
            legend(obj.output{1,1}, obj.output{2,1})
            hold off
            xlabel("/rho")

           

            E_kin_s= 0.5 * (obj.output{1,3}(:,4).^2 + obj.output{1,3}(:,6).^2 + obj.output{1,3}(:,5).^2);
            E_kin_e= 0.5* (obj.output{2,3}(:,4).^2 +obj.output{2,3}(:,6).^2 + obj.output{2,3}(:,5).^2);
            E_kin_sy_s= 0.5*(obj.output{1,3}(1,4).^2 + obj.output{1,3}(1,6).^2 + obj.output{1,3}(1,5).^2);
            E_kin_sy_e= 0.5*(obj.output{2,3}(1,4).^2 + obj.output{2,3}(1,6).^2 + obj.output{2,3}(1,5).^2);

            figure("Name",'kinetic enery particles with respect to synch. mol.(subtracted E_{kin,synch}')
            hold on
            histogram(E_kin_s-E_kin_sy_s,100)
            histogram(E_kin_s(~ismember(obj.output{1,4},obj.output{2,4}))-E_kin_sy_s,100);
            xlabel('E_{kin} particles (m=1) (arb.units)')
            legend('start', 'start-end')
            hold off
            figure("Name",'E_{kin} particles with respect to synch. mol. (divided by E_{kin,synch})')
            hold on
            histogram(E_kin_s./E_kin_sy_s,100)
            histogram(E_kin_s(~ismember(obj.output{1,4},obj.output{2,4}))./E_kin_sy_s,100);
            xlabel('E_{kin} particles (m=1) (arb.units)')
            legend('start', 'start-end')
            hold off

            % Compare the phase spaces of each axis x,y,z between the start and end of the deccelerator
            figure('Name', 'Compare Phase spaces x,v_x/y,v_y/z,v_z  start decelerator. ,end decelerator.')
            subplot(3,1,1)
            title('Phase space x, vx')
            hold on 
            plot(obj.output{1,3}(:,1)-obj.output{1,3}(1,1),obj.output{1,3}(:,4)-obj.output{1,3}(1,4), '.')
            plot(obj.output{2,3}(:,1)-obj.output{2,3}(1,1),obj.output{2,3}(:,4)-obj.output{2,3}(1,4),'.r')
            legend('start', 'end')
            hold off
            xlabel('x_{pos.} (m)'), ylabel('v_x (m/s)')
            subplot(3,1,2)
            title('y,v_y')
            hold on 
            plot(obj.output{1,3}(:,2)-obj.output{1,3}(1,2),obj.output{1,3}(:,5)-obj.output{1,3}(1,5), '.')
            plot(obj.output{2,3}(:,2)-obj.output{2,3}(1,2),obj.output{2,3}(:,5)-obj.output{2,3}(1,5),'.r')
            legend('start', 'end')
            hold off
            xlabel('y_{pos.} (m)'), ylabel('v_y (m/s)')
            subplot(3,1,3)
            title('z,v_z')
            hold on 
            plot(obj.output{1,3}(:,3)-obj.output{1,3}(1,3),obj.output{1,3}(:,6)-obj.output{1,3}(1,6), '.')
            plot(obj.output{2,3}(:,3)-obj.output{2,3}(1,3),obj.output{2,3}(:,6)-obj.output{2,3}(1,6),'.r')
            legend('start', 'end')
            hold off
            xlabel('z_{pos.} (m)'), ylabel('v_z (m/s)')
        end

        function [h,area,TOF_cut]=gain_TOF(obj) % save max TOF as well as area to compare gfain factor to experiment
            h=max(obj.TOF_save(:,2));
            particle_entries=obj.TOF_save(:,2)>0;  % make array with logical one where particles are detected 0 otherwise
            TOF_cut=obj.TOF_save(particle_entries,:); % make a variable of cut TOF where entries are non zero also for compare plots
            [area,~]= size(obj.output{end-1,3});   % integrate area of TOF to compare to experiments and check performance        
        end

        function plot_TOF_laser(obj)
            plot(obj.TOF_save(:,1)*1e6+obj.params.FLY_incoupling_time*1e6,obj.TOF_save(:,2))
            xlabel('time (/mu s'); ylabel('# particles in laser volume');
            ylim([0,max(obj.TOF_save(:,2)+2)]);
            xlim([obj.TOF_save(1,1)*1e6+obj.params.FLY_incoupling_time*1e6,obj.TOF_save(end,1)*1e6+obj.params.FLY_incoupling_time*1e6]);
        end

        
            

    end
end