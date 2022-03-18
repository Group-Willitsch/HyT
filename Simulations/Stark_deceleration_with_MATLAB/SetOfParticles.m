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
        output                     % save xyzVxyz , t and flag at certain times and display here
        
        traj_xyzVxyz % trajectories x,y,z, Vx, Vy, Vz
        traj_time
        ind_particles                 % index of particles created at begining saving the indexes of the ones that make it to the end
        TOF_xyzVxyz                   % variable to save xyzVxyz at the time step/steps of detetcion
        TOF_save     
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
            rng('default') % fix the seed of the random nunber generator to bremoved afterwards
            obj.xyzVxyz_0 = randn(obj.num_particles, 6)/sqrt(8*log(2)).*... %conversion factor to std from full width half max, below adjust normal dist. to represent particles
                [obj.myInput.params.BEAM_long_pos_spread, obj.myInput.params.BEAM_radius_of_nozzle, obj.myInput.params.BEAM_radius_of_nozzle,...                            
                 obj.myInput.params.BEAM_long_vel_spread, obj.myInput.params.BEAM_trans_velocity_spread, obj.myInput.params.BEAM_trans_velocity_spread] + ...
                [-obj.myInput.params.PHYS_valve_to_dec, 0, 0, obj.myInput.params.BEAM_avg_velocity_beam, 0, 0]; %set v_x to avergae and set all particles to valve pos. (decc x=0)
            obj.xyzVxyz_0(1,:) = [-obj.myInput.params.PHYS_valve_to_dec, 0, 0, 450, 0, 0]; % The first row is for a synchronous molecule
              
        end

        
         %% Propagate particles, intergate using time intervals in M_time_vector, propagating all the molecules together
        function propagateParticles_euler(obj)
            
            % a function that removes the lost molecules
            function [xyzVxyz,ind_particles] = removeHitParticles(xyzVxyz,ind_particles)
                hit_indices = xyzVxyz(:,1) < obj.myInput.params.PHYS_length_dec & (abs(xyzVxyz(:,2)) > obj.myInput.params.PHYS_seperation_pins/2 | abs(xyzVxyz(:,3)) > obj.myInput.params.PHYS_seperation_pins/2);
                xyzVxyz = xyzVxyz(~hit_indices,:);
                ind_particles = ind_particles(~hit_indices);
                hit_indices = (abs(xyzVxyz(:,1) - xyzVxyz(1,1)) < 2* 5.5e-3); % separated the cut of the fish in x into two lines in order to make indexing work and also to be sure its correct
                xyzVxyz = xyzVxyz(hit_indices, :);
                ind_particles = ind_particles(hit_indices);
            end
            
            
            % first propagate to entrance of decelerator
            obj.xyzVxyz = obj.xyzVxyz_0;
            obj.ind_particles = transpose(1:size(obj.xyzVxyz,1));
            obj.xyzVxyz(:,1:3) = obj.xyzVxyz(:,1:3) + obj.xyzVxyz(:,4:6) * obj.myInput.params.FLY_incoupling_time;
             [obj.xyzVxyz,obj.ind_particles]= removeHitParticles(obj.xyzVxyz, obj.ind_particles );% select those that can enter dec
%             hit_indices = obj.xyzVxyz(:,1) < obj.myInput.params.PHYS_length_dec & (abs(obj.xyzVxyz(:,2)) > obj.myInput.params.PHYS_seperation_pins/2 | abs(obj.xyzVxyz(:,3)) > obj.myInput.params.PHYS_seperation_pins/2);
%             obj.xyzVxyz = obj.xyzVxyz(~hit_indices,:);
%             obj.xyzVxyz = obj.xyzVxyz((abs(obj.xyzVxyz(:,1) - obj.xyzVxyz(1,1)) < 2* 5.5e-3), :);

            % propagate inside dec
            
            dt = 5e-8;
            obj.output={};
%             fprintf("num/total switching\n");
            obj.Snapshot("start decelerator",obj.myInput.M_time_vec(1),obj.xyzVxyz, obj.ind_particles)
            if obj.myInput.params.FLY_focusing_mode_bool % focusing mode
                dxyzVxyz = {@(y) [y(:, 4:6), obj.myInput.ax_norm_interpl(y(:, 1), y(:, 2), y(:, 3)),...
                                                 [0;obj.myInput.ay_norm_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))],...
                                                 [0;obj.myInput.az_norm_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))]]                   
                            @(y) [y(:, 4:6), obj.myInput.ax_neg_H_interpl(y(:, 1), y(:, 3), -y(:, 2)),...
                                                 [0;-obj.myInput.az_neg_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))],...
                                                 [0;obj.myInput.ay_neg_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))]]
                            @(y) [y(:, 4:6), obj.myInput.ax_norm_H_interpl(y(:, 1), y(:, 3), -y(:, 2)),...
                                                 [0;-obj.myInput.az_norm_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))],...
                                                 [0;obj.myInput.ay_norm_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))]]
                            @(y) [y(:, 4:6), obj.myInput.ax_pos_interpl(y(:, 1), y(:, 2), y(:, 3)),...
                                                 [0;obj.myInput.ay_pos_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))],...
                                                 [0;obj.myInput.az_pos_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))]]
                            @(y) [y(:, 4:6), obj.myInput.ax_norm_interpl(y(:, 1), y(:, 2), y(:, 3)),...
                                                 [0;obj.myInput.ay_norm_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))],...
                                                 [0;obj.myInput.az_norm_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))]]
                            @(y) [y(:, 4:6), obj.myInput.ax_pos_H_interpl(y(:, 1), y(:, 3), -y(:, 2)),...
                                                 [0;-obj.myInput.az_pos_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))],...
                                                 [0;obj.myInput.ay_pos_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))]]
                            @(y) [y(:, 4:6), obj.myInput.ax_norm_H_interpl(y(:, 1), y(:, 3), -y(:, 2)),...
                                                 [0;-obj.myInput.az_norm_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))],...
                                                 [0;obj.myInput.ay_norm_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))]]
                            @(y) [y(:, 4:6), obj.myInput.ax_neg_interpl(y(:, 1), y(:, 2), y(:, 3)),...
                                                 [0;obj.myInput.ay_neg_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))],...
                                                 [0;obj.myInput.az_neg_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))]]};

                for i = 1:1: (length(obj.myInput.M_time_vec) - 2) %since free propagation is done with euler
                    if size(obj.xyzVxyz,1) == 0
                        error("All the molecules are lost after the %d switching!",i-1)
                    end
%                     fprintf("%d/%d\t%d\n",size(obj.xyzVxyz,1), obj.num_particles, i);
                    
                    for t = obj.myInput.M_time_vec(i):dt: obj.myInput.M_time_vec(i+1)
                        obj.xyzVxyz = obj.xyzVxyz + dxyzVxyz{mod(i, 8)+8*(~mod(i, 8))}(obj.xyzVxyz) * dt;
                    end
                    [obj.xyzVxyz,obj.ind_particles]= removeHitParticles(obj.xyzVxyz, obj.ind_particles );
                end
                fprintf("%d out of %d particles left\n",size(obj.xyzVxyz,1), obj.num_particles);
            else    % normal mode
                dxyzVxyz = {@(y) [y(:, 4:6), obj.myInput.ax_norm_interpl(y(:, 1), y(:, 2), y(:, 3)),...
                                                 [0;obj.myInput.ay_norm_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))],...
                                                 [0;obj.myInput.az_norm_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))]]                   
                            @(y) [y(:, 4:6), obj.myInput.ax_norm_H_interpl(y(:, 1), y(:, 3), -y(:, 2)),...
                                                 [0;-obj.myInput.az_norm_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))],...
                                                 [0;obj.myInput.ay_norm_H_interpl(y(2:end, 1), y(2:end, 3), -y(2:end, 2))]]};
                
                for i = 1: 1: (length(obj.myInput.M_time_vec) - 2)
                    if size(obj.xyzVxyz,1) == 0
                        error("All the molecules are lost after the %d switching!", i-1)
                    end
%                     fprintf("%d/%d\t%d\n",size(obj.xyzVxyz,1), obj.num_particles,i);
                    for t = obj.myInput.M_time_vec(i):dt: obj.myInput.M_time_vec(i+1)
                        obj.xyzVxyz = obj.xyzVxyz + dxyzVxyz{mod(i, 2)+2*(~mod(i, 2))}(obj.xyzVxyz) * dt;
                    end
                    [obj.xyzVxyz,obj.ind_particles]= removeHitParticles(obj.xyzVxyz, obj.ind_particles );
                end
            end
            obj.Snapshot('end decelerator',t(end), obj.xyzVxyz, obj.ind_particles)
            if obj.myInput.verbose

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
            end
            obj.TOF()                                                   
                        
            
%           free propagation newton since velocity stays same thus
%           only update xyz by using x_new= x_old + v_x*dt for x,y,z
            obj.xyzVxyz(:,1:3)=obj.xyzVxyz(:,1:3)+obj.xyzVxyz(:,4:6)*(obj.myInput.M_time_vec(end)-obj.myInput.M_time_vec(end-1));                                                            
            figure();
            scatter(obj.xyzVxyz(:,1)*10^3, obj.xyzVxyz(:,4));  %scatter(x,y) creates a scatter plot with circular markers at the locations specified by the vectors x and y.
            xlabel('x-pos (mm)'); ylabel('v_x (m/s)');
            
            figure;
            obj.arrival_time = t(end) - (obj.xyzVxyz(:,1) - obj.myInput.params.PHYS_length_dec - obj.myInput.params.PHYS_exit_to_detection)./(obj.xyzVxyz(:,4)) + obj.myInput.params.FLY_incoupling_time;
            histogram(obj.arrival_time, 100);xlabel('time(s)')

            figure;
            scatter(obj.xyzVxyz(:,1), obj.xyzVxyz(:,4));
            
        end
        
         %% Propagate particles, intergate using time intervals in M_time_vector, propagating all the molecules together
        function propagateParticles_ode45(obj)
            
            % a function that remove the lost molecules
        function [xyzVxyz,ind_particles] = removeHitParticles(xyzVxyz,ind_particles)
                hit_indices = xyzVxyz(:,1) < obj.myInput.params.PHYS_length_dec & (abs(xyzVxyz(:,2)) > obj.myInput.params.PHYS_seperation_pins/2 | abs(xyzVxyz(:,3)) > obj.myInput.params.PHYS_seperation_pins/2);
                xyzVxyz = xyzVxyz(~hit_indices,:);
                ind_particles = ind_particles(~hit_indices);
                hit_indices = (abs(xyzVxyz(:,1) - xyzVxyz(1,1)) < 2* 5.5e-3); % separated the cut of the fish in x into two lines in order to make indexing work and also to be sure its correct
                xyzVxyz = xyzVxyz(hit_indices, :);
                ind_particles = ind_particles(hit_indices);
        end        
            
            % first propagate to entrance of decelerator
            obj.xyzVxyz = obj.xyzVxyz_0;
            obj.ind_particles = transpose(1:size(obj.xyzVxyz,1));
            obj.xyzVxyz(:,1:3) = obj.xyzVxyz(:,1:3) + obj.xyzVxyz(:,4:6) * obj.myInput.params.FLY_incoupling_time;
            [obj.xyzVxyz,obj.ind_particles]= removeHitParticles(obj.xyzVxyz, obj.ind_particles );% select those that can enter dec
           

            % propagate inside dec
            fprintf("num/total switching\n");
            opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
            times = [];   % used for commented plots below which tell you plot you the time steps of the integration
            num_par = [];
            obj.output={};
            obj.Snapshot("start decelerator",obj.myInput.M_time_vec(1),obj.xyzVxyz, obj.ind_particles)
            if obj.myInput.params.FLY_focusing_mode_bool % focusing mode
                dydt_array = {@(t,y) obj.dydtNormVerticalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtNegHorizontalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtNormHorizontalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtPosVerticalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtNormVerticalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtPosHorizontalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtNormHorizontalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtNegVerticalOn(t,y,size(obj.xyzVxyz,1))};
                tic
                for i = 1:1: (length(obj.myInput.M_time_vec) - 2)
                    if size(obj.xyzVxyz,1) == 0
                        error("All the molecules are lost after the %d switching!",i-1)
                    end
                    fprintf("%d/%d\t%d\n",size(obj.xyzVxyz,1), obj.num_particles, i);
                    [t, y] = ode45(dydt_array{mod(i, 8)+8*(~mod(i, 8))}, [obj.myInput.M_time_vec(i), (obj.myInput.M_time_vec(i) + obj.myInput.M_time_vec(i+1))/2, obj.myInput.M_time_vec(i+1)], reshape(obj.xyzVxyz, [], 1), opts);
                    obj.xyzVxyz = reshape(y(end,:), [], 6);
                    [obj.xyzVxyz,obj.ind_particles]= removeHitParticles(obj.xyzVxyz,obj.ind_particles);
                    fprintf("%d/%d\t%d\n",size(obj.xyzVxyz,1), obj.num_particles,i);
                end
                fprintf("%d out of %d particles left\n",size(obj.xyzVxyz,1), obj.num_particles);
            else    % normal mode
                dydt_array = {@(t,y) obj.dydtNormVerticalOn(t,y,size(obj.xyzVxyz,1)),...
                              @(t,y) obj.dydtNormHorizontalOn(t,y,size(obj.xyzVxyz,1))};
                for i = 1: 1: (length(obj.myInput.M_time_vec) - 2)
                    if size(obj.xyzVxyz,1) == 0
                        error("All the molecules are lost after the %d switching!", i-1)
                    end
%                     fprintf("%d/%d\t%d\n",size(obj.xyzVxyz,1), obj.num_particles,i);
                    [t, y] = ode45(dydt_array{mod(i, 2)+2*(~mod(i, 2))}, [obj.myInput.M_time_vec(i), (obj.myInput.M_time_vec(i) + obj.myInput.M_time_vec(i+1))/2, obj.myInput.M_time_vec(i+1)], reshape(obj.xyzVxyz, [], 1), opts);             
%                   times = [times toc];
%                   num_par = [num_par, size(obj.xyzVxyz,1)]; 
                    obj.xyzVxyz = reshape(y(end,:), [], 6);
                    [obj.xyzVxyz,obj.ind_particles]= removeHitParticles(obj.xyzVxyz,obj.ind_particles);
                    fprintf("%d/%d\t%d\n",size(obj.xyzVxyz,1), obj.num_particles,i);
                end
                toc
            end
            obj.Snapshot('end decelerator',t(end), obj.xyzVxyz, obj.ind_particles)
            
            if obj.myInput.verbose

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
            end
            obj.TOF()                                                   
                        
            
%           free propagation newton since velocity stays same thus
%           only update xyz by using x_new= x_old + v_x*dt for x,y,z
            obj.xyzVxyz(:,1:3)=obj.xyzVxyz(:,1:3)+obj.xyzVxyz(:,4:6)*(obj.myInput.M_time_vec(end)-obj.myInput.M_time_vec(end-1));                                                            
            figure();
            scatter(obj.xyzVxyz(:,1)*10^3, obj.xyzVxyz(:,4));  %scatter(x,y) creates a scatter plot with circular markers at the locations specified by the vectors x and y.
            xlabel('x-pos (mm)'); ylabel('v_x (m/s)');
            
            figure;
            obj.arrival_time = t(end) - (obj.xyzVxyz(:,1) - obj.myInput.params.PHYS_length_dec - obj.myInput.params.PHYS_exit_to_detection)./(obj.xyzVxyz(:,4)) + obj.myInput.params.FLY_incoupling_time;
            histogram(obj.arrival_time, 100);xlabel('time(s)')
            
            
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
            
            
            obj.arrival_time = obj.arrival_time(abs(obj.xyzVxyz(:,2)) < obj.myInput.params.FLY_detection_laser_diameter/2 & abs(obj.xyzVxyz(:,3)) < 2e-3);
            bin_begin_each_particle = int16((obj.arrival_time - obj.myInput.params.FLY_detection_laser_diameter/2/obj.xyzVxyz(:,4) - min_time)/binsize + 0.5);
            bin_end_each_particle =int16((obj.arrival_time + obj.myInput.params.FLY_detection_laser_diameter/2/obj.xyzVxyz(:,4) - min_time)/binsize + 0.5);
            
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

         %% Plot TOF
        % redone with laser volume and snapshot of xyzVxyz when snych. molecule is at detection point
        function TOF(obj) 
            t_profile=obj.myInput.M_time_vec(end-1):1e-6:obj.myInput.M_time_vec(end)+1e-3;
            t_steps=diff(t_profile); %time srteps for gfree propagation one smaller than t_profile
            x_laser= obj.myInput.params.PHYS_length_dec + obj.myInput.params.PHYS_exit_to_detection;  % x coordinate laser
            h_laser= 2*1e-3; %half of height laser volume to make check since height is 4 mm but goes from y=-2mm to y=2mm
            r_laser=obj.myInput.params.FLY_detection_laser_diameter/2; % radius of laser we model as cylinder
            obj.TOF_save= zeros(length(t_steps), 2); % matrix wher number of molecules in laser volume will be saved with the corresponding time
            obj.TOF_save(:,1)=t_profile(2:end);
            xyzVxyz_TOF=obj.xyzVxyz;    % make local variable since we use then obj.xyzVxyz to take pahse space picture at M_time_vec(end) (time synch reaches detection spot middle laser?)


            
            for i=1:length(t_steps)
                xyzVxyz_TOF(:,1:3)=xyzVxyz_TOF(:,1:3)+xyzVxyz_TOF(:,4:6)*t_steps(i); % let particles propagate 1 mu s after decc.
                in_volume= abs(xyzVxyz_TOF(:,2))<= h_laser & sqrt((xyzVxyz_TOF(:,1)-x_laser).^2 + (xyzVxyz_TOF(:,3)).^2) <= r_laser;   %check if a particle is in y range of cylinder (array of ones and zeros) one if inside range
               number_part_det=sum(in_volume);
               if number_part_det>0    % save the time at which we checked if they are inside
                  obj.TOF_save(i,2)= number_part_det; %save the number of detected particles at that time
                  obj.TOF_save(i,1)=t_profile(i+1);
               else
                  obj.TOF_save(i,1)=t_profile(i+1);
               end
%                would be possible to take snap shots of partciles here at
%                wanted detecion stages
%                if round(t_profile(i),5) == round(obj.myInput.M_time_vec(end),5)
%                    obj.Snapshot('end laser',t_profile(i+1), xyzVxyz_TOF(in_volume,:),1); % index saving not implemented for this fucntion
%                end

            end
            
            figure('Name','TOF')
            plot(obj.TOF_save(:,1)*1e6,obj.TOF_save(:,2))
            xlabel('time (/mu s'); ylabel('# particles in laser volume');
            ylim([0,max(obj.TOF_save(:,2)+2)]);
            xlim([obj.TOF_save(1,1)*1e6,obj.TOF_save(end,1)*1e6]);
        end

        %fucntion taking a snapshot of time, xyzVxyz and the corresponding
        %indexes of the particles as well a flag that tells you the point
        function Snapshot(obj, flag, t, pos_vel, index_p)
            obj.output= [obj.output; {flag,t,pos_vel,index_p}];
        end

        function Snapshotplot(obj)
            [leng,~]=size(obj.output);

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
        
    
%             figure() % if one wants to see a bar plot of the indexes showing which molecules made it into the decelerator
%             bar(obj.output{1,4},ismember(obj.output{1,4},obj.output{2,4}))
%             xlabel('Index'); ylabel('present at end deacc.'); ylim([0,2])

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


            

            
%             figure()
%             subplot(2,3,1)
%             title('Phase space x,v_x start')
%             hold on 
%             plot(obj.output{1,3}(:,1)*10^3,obj.output{1,3}(:,4)*10^3, '.')
%             plot(obj.output{1,3}(1,1)*10^3,obj.output{1,3}(1,4)*10^3,'.r')
%             hold off
%             xlabel('x_{pos.} (mm)'), ylabel('v_x (m/s)')
%             subplot(2,3,2)
%             title('Phase space y,v_y start')
%             hold on 
%             plot(obj.output{1,3}(:,2)*10^3,obj.output{1,3}(:,5)*10^3, '.')
%             plot(obj.output{1,3}(1,2)*10^3,obj.output{1,3}(1,5)*10^3,'.r')
%             hold off
%             xlabel('y_{pos.} (mm)'), ylabel('v_y (m/s)')
%             subplot(2,3,3)
%             title('Phase space z, v_z start')
%             hold on 
%             plot(obj.output{1,3}(:,3)*10^3,obj.output{1,3}(:,6)*10^3, '.')
%             plot(obj.output{1,3}(1,3)*10^3,obj.output{1,3}(1,6)*10^3,'.r')
%             hold off
%             xlabel('z_{pos.} (mm)'), ylabel('v_z (m/s)')
%             subplot(2,3,4)
%             title('Phase space x, v_x end')
%             hold on 
%             plot(obj.output{2,3}(:,1)*10^3,obj.output{2,3}(:,4)*10^3, '.')
%             plot(obj.output{2,3}(1,1)*10^3,obj.output{2,3}(1,4)*10^3,'.r')
%             hold off
%             xlabel('x_{pos.} (mm)'), ylabel('v_x (m/s)')
%             subplot(2,3,5)
%             title('Phase space y, v_y end')
%             hold on 
%             plot(obj.output{2,3}(:,2)*10^3,obj.output{2,3}(:,5)*10^3, '.')
%             plot(obj.output{2,3}(1,2)*10^3,obj.output{2,3}(1,5)*10^3,'.r')
%             hold off
%             xlabel('y_{pos.} (mm)'), ylabel('v_y (m/s)')
%             subplot(2,3,6)
%             title('Phase space z, v_z end')
%             hold on 
%             plot(obj.output{2,3}(:,3)*10^3,obj.output{2,3}(:,6)*10^3, '.')
%             plot(obj.output{2,3}(1,3)*10^3,obj.output{2,3}(1,6)*10^3,'.r')
%             hold off
%             xlabel('z_{pos.} (mm)'), ylabel('v_z (m/s)')


            

        end  


    end
end


