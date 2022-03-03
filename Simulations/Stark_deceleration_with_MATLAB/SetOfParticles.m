classdef SetOfParticles < handle  % handle makes it works more like a Python class than a C++ class, but dangerous!!
    %This class holds the 1e6 or so different particles you shoot across
    %the decelerator, for a given FIXED set of input and general settings.
    %This class will:
%     - initialize a set of number_of_particles
%     - call the numerical solver (with possible parallelization here)
%     - plot phase space at the beggining and at the end
%     - plot trayectories and 3d plots of general interest
%     - plot everything else you may think of
%     - slim down the trayectories you may want to save, or keep the relevant results of the simulation (like TOF profile) in internal and .mat variables.

    properties
        myInput
        num_particles
        xyzVxyz_0               % the vector of initial pos&vel, number_of_particles * 6
        xyzVxyz                 % the vector of computed pos&vel, number_of_particles * 6
        has_the_simulation_been_run     % boolean, default is False in constructor, will be switched to True at the end of the simulation
        num_trajectories_saved
        arrival_time
        
        traj_xyzVxyz % trajectories x,y,z, Vx, Vy, Vz
        traj_time
    end
    
    methods
        %% CONSTRUCTOR
        function obj = SetOfParticles(InputParameters)
            if nargin >= 1
                obj.myInput = InputParameters;
            end
            
            obj.num_particles = 1;
            obj.has_the_simulation_been_run = false;
            obj.num_trajectories_saved = 100;
            obj.createParticles();
        end
        
        %% create particles
        
        function createParticles(obj)
            obj.xyzVxyz_0 = randn(obj.num_particles, 6)/sqrt(8*log(2)).*...
                [obj.myInput.params.BEAM_long_pos_spread, obj.myInput.params.BEAM_radius_of_nozzle, obj.myInput.params.BEAM_radius_of_nozzle,...                            
                 obj.myInput.params.BEAM_long_vel_spread, obj.myInput.params.BEAM_trans_velocity_spread, obj.myInput.params.BEAM_trans_velocity_spread] + ...
                [-obj.myInput.params.PHYS_valve_to_dec, 0, 0, obj.myInput.params.BEAM_avg_velocity_beam, 0, 0];
            obj.xyzVxyz_0(1,:) = [-obj.myInput.params.PHYS_valve_to_dec, 0, 0, 450, 0, 0]; % The first row is for a synchronous molecule
              
        end
        
         %% Propagate particles, intergate using time intervals in M_time_vector, propagating all the molecules together
        function propagateParticles(obj)
            
            % a function that remove the lost molecules
            function xyzVxyz = removeHitParticles(xyzVxyz)
                hit_indices = xyzVxyz(:,1) < obj.myInput.params.PHYS_length_dec & (abs(xyzVxyz(:,2)) > obj.myInput.params.PHYS_seperation_pins/2 | abs(xyzVxyz(:,3)) > obj.myInput.params.PHYS_seperation_pins/2);
                xyzVxyz = xyzVxyz(~hit_indices,:);
                xyzVxyz = xyzVxyz((abs(xyzVxyz(:,1) - xyzVxyz(1,1)) < 2* 5.5e-3), :);
            end
            
            % first propagate to entrance of decelerator
            obj.xyzVxyz = obj.xyzVxyz_0;
            obj.xyzVxyz(:,1:3) = obj.xyzVxyz(:,1:3) + obj.xyzVxyz(:,4:6) * obj.myInput.params.FLY_incoupling_time;
            obj.xyzVxyz= removeHitParticles(obj.xyzVxyz);% select those that can enter dec

            % propagate inside dec
            fprintf("num/total switching\n");
            opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
            if obj.myInput.params.FLY_focusing_mode_bool % focusing mode
                dydt_array = {@(t,y) obj.dydtNormVerticalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtNegHorizontalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtNormHorizontalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtPosVerticalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtNormVerticalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtPosHorizontalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtNormHorizontalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtNegVerticalOn(t,y,size(obj.xyzVxyz,1))};
                 
                for i = 1:1: (length(obj.myInput.M_time_vec) - 1)
                    if size(obj.xyzVxyz,1) == 0
                        error("All the molecules are lost after the %d switching!",i-1)
                    end
                    fprintf("%d/%d\t%d\n",size(obj.xyzVxyz,1), obj.num_particles, i);
                    [t, y] = ode45(dydt_array{mod(i, 8)+8*(~mod(i, 8))}, [obj.myInput.M_time_vec(i), (obj.myInput.M_time_vec(i) + obj.myInput.M_time_vec(i+1))/2, obj.myInput.M_time_vec(i+1)], reshape(obj.xyzVxyz, [], 1), opts);
                    obj.traj_time = [obj.traj_time; t]; 
                    obj.xyzVxyz = reshape(y(end,:), [], 6);
                    obj.xyzVxyz= removeHitParticles(obj.xyzVxyz);
                end
                
            else    % normal mode
                dydt_array = {@(t,y) obj.dydtNormVerticalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtNormHorizontalOn(t,y,size(obj.xyzVxyz,1))};
                for i = 1: 1: (length(obj.myInput.M_time_vec) - 1)
                    if size(obj.xyzVxyz,1) == 0
                        error("All the molecules are lost after the %d switching!", i-1)
                    end
                    fprintf("%d/%d\t%d\n",size(obj.xyzVxyz,1), obj.num_particles,i);
                    [t, y] = ode45(dydt_array{mod(i, 2)+2*(~mod(i, 2))}, [obj.myInput.M_time_vec(i), (obj.myInput.M_time_vec(i) + obj.myInput.M_time_vec(i+1))/2, obj.myInput.M_time_vec(i+1)], reshape(obj.xyzVxyz, [], 1), opts);             
                    obj.xyzVxyz = reshape(y(end,:), [], 6);
                    obj.xyzVxyz= removeHitParticles(obj.xyzVxyz);
                end
            end
            
            figure;
            obj.arrival_time = t(end) - (obj.xyzVxyz(:,1) - obj.myInput.params.PHYS_length_dec - obj.myInput.params.PHYS_exit_to_detection)./(obj.xyzVxyz(:,4)) + obj.myInput.params.FLY_incoupling_time;
            histogram(obj.arrival_time, 100);
            
            
            figure;
            scatter(obj.xyzVxyz(:,1), obj.xyzVxyz(:,4));
            
        end
        
        %% Obtain the acceleration by interpolation
        function dydt = dydtNormVerticalOn(obj,t,y,n)
            dydt = [y(3*n+1:end); obj.myInput.ax_norm_interpl(y(1:n), y(n+1:2*n), y(2*n+1:3*n)); [0;obj.myInput.ay_norm_interpl(y(2:n), y(n+2:2*n), y(2*n+2:3*n))]; [0; obj.myInput.az_norm_interpl(y(2:n), y(n+2:2*n), y(2*n+2:3*n))]]; 
        end
        
        function dydt = dydtNormHorizontalOn(obj,t,y,n)
            dydt = [y(3*n+1:end); obj.myInput.ax_norm_H_interpl(y(1:n), y(2*n+1:3*n), -y(n+1:2*n)); [0;-obj.myInput.az_norm_H_interpl(y(2:n), y(2*n+2:3*n), -y(n+2:2*n))]; [0;obj.myInput.ay_norm_H_interpl(y(2:n), y(2*n+2:3*n), -y(n+2:2*n))]]; 
        end
        
        function dydt = dydtNegVerticalOn(obj,t,y,n)
            dydt = [y(3*n+1:end); obj.myInput.ax_neg_interpl(y(1:n), y(n+1:2*n), y(2*n+1:3*n)); [0; obj.myInput.ay_neg_interpl(y(2:n), y(n+2:2*n), y(2*n+2:3*n))]; [0; obj.myInput.az_neg_interpl(y(2:n), y(n+2:2*n), y(2*n+2:3*n))]];
        end
        
        function dydt = dydtNegHorizontalOn(obj,t,y,n)
            dydt = [y(3*n+1:end); obj.myInput.ax_neg_H_interpl(y(1:n), y(2*n+1:3*n), -y(n+1:2*n)); [0;-obj.myInput.az_neg_H_interpl(y(2:n), y(2*n+2:3*n), -y(n+2:2*n))]; [0;obj.myInput.ay_neg_H_interpl(y(2:n), y(2*n+2:3*n), -y(n+2:2*n))]]; 
        end

        function dydt = dydtPosVerticalOn(obj,t,y,n)
            dydt = [y(3*n+1:end); obj.myInput.ax_pos_interpl(y(1:n), y(n+1:2*n), y(2*n+1:3*n)); [0; obj.myInput.ay_pos_interpl(y(2:n), y(n+2:2*n), y(2*n+2:3*n))]; [0; obj.myInput.az_pos_interpl(y(2:n), y(n+2:2*n), y(2*n+2:3*n))]];
        end
        
        function dydt = dydtPosHorizontalOn(obj,t,y,n)
            dydt = [y(3*n+1:end); obj.myInput.ax_pos_H_interpl(y(1:n), y(2*n+1:3*n), -y(n+1:2*n)); [0;-obj.myInput.az_pos_H_interpl(y(2:n), y(2*n+2:3*n), -y(n+2:2*n))]; [0;obj.myInput.ay_pos_H_interpl(y(2:n), y(2*n+2:3*n), -y(n+2:2*n))]]; 
        end

        
        function propagateParticlesAndSaveTrajectories(obj)
            
            % a function that remove the lost molecules
            function xyzVxyz = removeHitParticles(y)
                xyzVxyz = reshape(y(end,:), [], 6);
                remaining_indices = ((xyzVxyz(:,1) < obj.myInput.params.PHYS_length_dec & abs(xyzVxyz(:,2)) < obj.myInput.params.PHYS_seperation_pins/2 & abs(xyzVxyz(:,3)) < obj.myInput.params.PHYS_seperation_pins/2)...
                                    | xyzVxyz(:,1) >= obj.myInput.params.PHYS_length_dec) & abs(xyzVxyz(:,1) - xyzVxyz(1,1)) < 2*obj.myInput.params.PHYS_distance_stages;
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
                obj.traj_xyzVxyz(:,:,1) = obj.xyzVxyz_0;
            else
                fprintf("obj.num_trajectories == 0, no trajectories will be saved!");
            end
            
            % first propagate to entrance of decelerator
            obj.xyzVxyz = obj.xyzVxyz_0;
            obj.xyzVxyz(:,1:3) = obj.xyzVxyz(:,1:3) + obj.xyzVxyz(:,4:6) * obj.myInput.params.FLY_incoupling_time;
%             obj.xyzVxyz= removeHitParticles(obj.xyzVxyz);% select those that can enter dec

            
            % propagate inside dec
            fprintf("num/total switching\n");
%             line_temp = 0;
            
            opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
            if obj.myInput.params.FLY_focusing_mode_bool % focusing mode
                dydt_array = {@(t,y) obj.dydtNormVerticalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtNegHorizontalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtNormHorizontalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtPosVerticalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtNormVerticalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtPosHorizontalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtNormHorizontalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtNegVerticalOn(t,y,size(obj.xyzVxyz,1))};
                 
                for i = 1:1:(length(obj.myInput.M_time_vec) - 1)
                    if size(obj.xyzVxyz,1) == 0
                        error("All the molecules are lost after the %d switching!",i-1)
                    end
                    fprintf("%d/%d\t%d\n",size(obj.xyzVxyz,1), obj.num_particles, i);
%                     [t, y] = ode45(dydt_array{mod(i, 8)+8*(~mod(i, 8))}, obj.myInput.M_time_vec(i):1e-6:obj.myInput.M_time_vec(i+1), reshape(obj.xyzVxyz, [], 1), opts);
                    [t, y] = ode45(dydt_array{mod(i, 8)+8*(~mod(i, 8))}, [obj.myInput.M_time_vec(i), (obj.myInput.M_time_vec(i) + obj.myInput.M_time_vec(i+1))/2, obj.myInput.M_time_vec(i+1)], reshape(obj.xyzVxyz, [], 1), opts);

                    obj.traj_time = [obj.traj_time; t + obj.myInput.params.FLY_incoupling_time];
                    obj.xyzVxyz= removeHitParticles(y);
                end
                
            else    % normal mode
                dydt_array = {@(t,y) obj.dydtNormVerticalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtNormHorizontalOn(t,y,size(obj.xyzVxyz,1))};
                for i = 1: 1: (length(obj.myInput.M_time_vec) - 1)
                    if size(obj.xyzVxyz,1) == 0
                        error("All the molecules are lost after the %d switching!", i-1)
                    end
                    fprintf("%d/%d\t%d\n",size(obj.xyzVxyz,1), obj.num_particles,i);
                    [t, y] = ode45(dydt_array{mod(i, 2)+2*(~mod(i, 2))}, [obj.myInput.M_time_vec(i), obj.myInput.M_time_vec(i+1)], reshape(obj.xyzVxyz, [], 1), opts);
%                     [t, y] = ode45(dydt_array{mod(i, 2)+2*(~mod(i, 2))}, [obj.myInput.M_time_vec(i), (obj.myInput.M_time_vec(i) + obj.myInput.M_time_vec(i+1))/2, obj.myInput.M_time_vec(i+1)], reshape(obj.xyzVxyz, [], 1), opts);
                    
                    obj.traj_time = [obj.traj_time; t];
                    obj.xyzVxyz= removeHitParticles(y);
                end
            end

            obj.arrival_time = t(end) - (obj.xyzVxyz(:,1) - obj.myInput.params.PHYS_length_dec - obj.myInput.params.PHYS_exit_to_detection)./(obj.xyzVxyz(:,4)) + obj.myInput.params.FLY_incoupling_time;
            
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
            
            bin_begin_each_particle = int16((obj.arrival_time - obj.myInput.params.FLY_detection_laser_diameter/2/obj.xyzVxyz(:,4) - min_time)/binsize + 0.5);
            bin_end_each_particle =int16((obj.arrival_time + obj.myInput.params.FLY_detection_laser_diameter/2/obj.xyzVxyz(:,4) - min_time)/binsize + 0.5);
            
            time = min_time:binsize:max_time;
            for i = 1: length(bin_begin_each_particle)
                tof_profile(bin_begin_each_particle(i):bin_end_each_particle(i)) = tof_profile(bin_begin_each_particle(i):bin_end_each_particle(i)) + 1;
            end
            
            disp(size(time));
            disp(size(tof_profile));
            plot(time(1:end-1), tof_profile);
        end
    end
end

