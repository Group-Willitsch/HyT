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
        in                              % InputParams class, to be removed maybe?
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
        %% CONSTRUCTOR
        function obj = SetOfParticles(in, num_particles, options)
            arguments % list of mandary and optional arguments
            in                          InputParameters            
            num_particles               int64 
            options.NumTrajToBeSaved    int64 = 100 % default number of trajectors to be saved, if you want to svae them
            end

            if nargin < 2 % (if ... will probably fail earlier...)
                 fprintf('Usage: SetOfParticles( InputParameters, num_particles, options)\n(Options are: \t NumTrajToBeSaved\)\n')
            end
            obj.in = in;
            obj.num_particles = num_particles;            
            obj.num_trajectories_saved = options.NumTrajToBeSaved;
            obj.has_the_simulation_been_run = false;
        end
        
        %% create particles
%         
        function createParticles(obj)
            rng('default') % fix the seed of the random nunber generator to bremoved afterwards
            obj.xyzVxyz_0 = randn(obj.num_particles, 6)/sqrt(8*log(2)).*... %conversion factor to std from full width half max, below adjust normal dist. to represent particles
                [obj.in.params.BEAM_long_pos_spread, obj.in.params.BEAM_radius_of_nozzle, obj.in.params.BEAM_radius_of_nozzle,...                            
                 obj.in.params.BEAM_long_vel_spread, obj.in.params.BEAM_trans_velocity_spread, obj.in.params.BEAM_trans_velocity_spread] + ...
                [-obj.in.params.PHYS_valve_to_dec, 0, 0, obj.in.params.BEAM_avg_velocity_beam, 0, 0]; %set v_x to avergae and set all particles to valve pos. (decc x=0)

            obj.xyzVxyz_0(1,:) = [-obj.in.params.PHYS_valve_to_dec, 0, 0, 450, 0, 0]; % The first row is for a synchronous molecule
            E_kin=((0.5*((obj.xyzVxyz_0(:,4)).^2 + (obj.xyzVxyz_0(:,5)).^2 + (obj.xyzVxyz_0(:,6)).^2))-0.5*(obj.xyzVxyz_0(1,4).^2 + obj.xyzVxyz_0(1,5).^2 + obj.xyzVxyz_0(1,6).^2))*1e-5;           
            sort_E_kin= [obj.xyzVxyz_0,abs(E_kin)];
            sort_E_kin= sortrows(sort_E_kin,7);
            obj.xyzVxyz_0=sort_E_kin(:,1:6);   
        end

%% Propagate particles, intergate using time intervals in M_time_vector, propagating all the molecules together
        function propagateParticles_euler(obj)
            ax_norm_interpl = obj.in.ax_norm_interpl;
            ay_norm_interpl = obj.in.ay_norm_interpl;
            az_norm_interpl = obj.in.az_norm_interpl;
            ax_norm_H_interpl =  obj.in.ax_norm_H_interpl;
            ay_norm_H_interpl = obj.in.ay_norm_H_interpl;
            az_norm_H_interpl = obj.in.az_norm_H_interpl;
            
            if obj.in.params.FLY_focusing_mode_bool
                 ax_neg_interpl = obj.in.ax_neg_interpl;
                 ay_neg_interpl = obj.in.ay_neg_interpl;
                 az_neg_interpl = obj.in.az_neg_interpl;
                 ax_neg_H_interpl = obj.in.ax_neg_H_interpl;
                 ay_neg_H_interpl = obj.in.ay_neg_H_interpl;
                 az_neg_H_interpl = obj.in.az_neg_H_interpl;
                 ax_pos_interpl = obj.in.ax_pos_interpl;
                 ay_pos_interpl = obj.in.ay_pos_interpl;
                 az_pos_interpl = obj.in.az_pos_interpl;
                 ax_pos_H_interpl = obj.in.ax_pos_H_interpl;
                 ay_pos_H_interpl = obj.in.ay_pos_H_interpl;
                 az_pos_H_interpl = obj.in.az_pos_H_interpl;
            end
            
%             a function that removes the lost molecules
            function [xyzVxyz,ind_particles] = removeHitParticles(xyzVxyz,ind_particles)
                hit_indices = xyzVxyz(:,1) < obj.in.params.PHYS_length_dec & (abs(xyzVxyz(:,2)) > obj.in.params.PHYS_seperation_pins/2 | abs(xyzVxyz(:,3)) > obj.in.params.PHYS_seperation_pins/2);
                xyzVxyz = xyzVxyz(~hit_indices,:);
                ind_particles = ind_particles(~hit_indices);
                hit_indices = (abs(xyzVxyz(:,1) - xyzVxyz(1,1)) < 2* 5.5e-3); % separated the cut of the fish in x into two lines in order to make indexing work and also to be sure its correct
                xyzVxyz = xyzVxyz(hit_indices, :);
                ind_particles = ind_particles(hit_indices);
            end
            
            % first propagate to entrance of decelerator
            obj.xyzVxyz = obj.xyzVxyz_0;
            obj.ind_particles = transpose(1:size(obj.xyzVxyz,1));
            obj.xyzVxyz(:,1:3) = obj.xyzVxyz(:,1:3) + obj.xyzVxyz(:,4:6) * obj.in.params.FLY_incoupling_time; % ekin at start of deacc.
            [obj.xyzVxyz,obj.ind_particles]= removeHitParticles(obj.xyzVxyz, obj.ind_particles );% select those that can enter dec
%             hit_indices = obj.xyzVxyz(:,1) < obj.params.PHYS_length_dec & (abs(obj.xyzVxyz(:,2)) > obj.params.PHYS_seperation_pins/2 | abs(obj.xyzVxyz(:,3)) > obj.params.PHYS_seperation_pins/2);
%             obj.xyzVxyz = obj.xyzVxyz(~hit_indices,:);
%             obj.xyzVxyz = obj.xyzVxyz((abs(obj.xyzVxyz(:,1) - obj.xyzVxyz(1,1)) < 2* 5.5e-3), :);

            % propagate inside dec
            dt = 4e-8;
%             fprintf("num/total switching\n");
            obj.Snapshot("start decelerator",obj.in.M_time_vec(1),obj.xyzVxyz, obj.ind_particles)
            if obj.in.params.FLY_focusing_mode_bool % focusing mode
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

                for i = 1:1: (length(obj.in.M_time_vec) - 2) %since free propagation is done with euler
                    if size(obj.xyzVxyz,1) == 0
                        error("All the molecules are lost after the %d switching!",i-1)
                    end
%                     fprintf("%d/%d\t%d\n",size(obj.xyzVxyz,1), obj.num_particles, i);
                    for t = obj.in.M_time_vec(i):dt: obj.in.M_time_vec(i+1)
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
                
                for i = 1: 1: (length(obj.in.M_time_vec) - 2)
                    if size(obj.xyzVxyz,1) == 0
                        error("All the molecules are lost after the %d switching!", i-1)
                    end
%                     fprintf("%d/%d\t%d\n",size(obj.xyzVxyz,1), obj.num_particles,i);
                    for t = obj.in.M_time_vec(i):dt: obj.in.M_time_vec(i+1)
                        obj.xyzVxyz = obj.xyzVxyz + dxyzVxyz{mod(i, 2)+2*(~mod(i, 2))}(obj.xyzVxyz) * dt;
                    end
                    [obj.xyzVxyz,obj.ind_particles]= removeHitParticles(obj.xyzVxyz, obj.ind_particles );
                end
                fprintf("%d out of %d particles left\n",size(obj.xyzVxyz,1), obj.num_particles);
            end
            obj.Snapshot('end decelerator',t(end), obj.xyzVxyz, obj.ind_particles)
            obj.TOF()
            if obj.in.verbose         
            obj.Snapshotplot()
          
            obj.xyzVxyz(:,1:3)=obj.xyzVxyz(:,1:3)+obj.xyzVxyz(:,4:6)*(obj.in.M_time_vec(end)-obj.in.M_time_vec(end-1));                                                            
            figure();
            scatter(obj.xyzVxyz(:,1)*10^3, obj.xyzVxyz(:,4));  %scatter(x,y) creates a scatter plot with circular markers at the locations specified by the vectors x and y.
            xlabel('x-pos (mm)'); ylabel('v_x (m/s)');
            
            figure;
            obj.arrival_time = t(end) - (obj.xyzVxyz(:,1) - obj.in.params.PHYS_length_dec - obj.in.params.PHYS_exit_to_detection)./(obj.xyzVxyz(:,4)) + obj.in.params.FLY_incoupling_time;
            histogram(obj.arrival_time, 100);xlabel('time(s)')

            figure;
            scatter(obj.xyzVxyz(:,1), obj.xyzVxyz(:,4));%this compares and plots the phase spaces of snapshot (maybe make function sow that it could take any number of netries not jsut 2 as now)
            end

%           free propagation newton since velocity stays same thus
%           only update xyz by using x_new= x_old + v_x*dt for x,y,zphase
            
        end
        
        
         %% Propagate particles, intergate using time intervals in M_time_vector, propagating all the molecules together
         function propagateParticles_ode45(obj,in_l)
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
            
            
            obj.arrival_time = obj.arrival_time(abs(obj.xyzVxyz(:,2)) < obj.in.params.FLY_detection_laser_diameter/2 & abs(obj.xyzVxyz(:,3)) < 2e-3);
            bin_begin_each_particle = int16((obj.arrival_time - obj.in.params.FLY_detection_laser_diameter/2/obj.xyzVxyz(:,4) - min_time)/binsize + 0.5);
            bin_end_each_particle =int16((obj.arrival_time + obj.in.params.FLY_detection_laser_diameter/2/obj.xyzVxyz(:,4) - min_time)/binsize + 0.5);
            
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
            % t_profile = (obj.in.M_time_vec(end-1) + t_steps_TOF ):t_steps_TOF:(obj.in.M_time_vec(end)+1e-3); % start 1us after, go on for
            t_profile = t_steps_TOF : t_steps_TOF : (obj.in.M_time_vec(end) - obj.in.M_time_vec(end-1) + 1e-3); % start 1us after, go on for 1 ms
            
            x_laser = obj.in.params.PHYS_length_dec + obj.in.params.PHYS_exit_to_detection;  % x coordinate laser center
            h_laser = 2e-3; %half of height laser volume to make check since height is 4 mm but goes from y=-2mm to y=2mm
            r_laser = obj.in.params.FLY_detection_laser_diameter/2; % radius of laser we model as cylinder
            
            obj.TOF_save = zeros(length(t_profile), 2); % matrix wher number of molecules in laser volume will be saved with the corresponding time
            obj.TOF_save(:, 1) = t_profile + obj.in.M_time_vec(end-1);
            output_xyzVxyz = {};

            % propagate till synch molecules
            xyzVxyz_TOF = obj.xyzVxyz; % local copy not to modify class variable
            xyzVxyz_TOF(:, 1:3) = obj.xyzVxyz(:, 1:3) + obj.xyzVxyz(:,4:6)*(obj.in.M_time_vec(end) - obj.in.M_time_vec(end-1)); 
            in_volume =  abs( xyzVxyz_TOF(:,3) ) <= h_laser & sqrt((xyzVxyz_TOF(:,1)-x_laser).^2 + (xyzVxyz_TOF(:,2)).^2) <= r_laser;
            xyzVxyz_TOF = xyzVxyz_TOF(in_volume, :);
            obj.Snapshot('synch. molecule detection', obj.in.M_time_vec(end), xyzVxyz_TOF,[]);


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

            obj.Snapshot('end decc. detected (safe)',obj.in.M_time_vec(end-1), obj.xyzVxyz(logical(xyzVxyz_in_volume),:), []) % here time is chosen at end decc. since pos. and vel. xyzVxyz correspond to this time
           
            compare=output_xyzVxyz{1,1}; % idea below is to compare velocities of all particles detected around the synch molecule se above if statement such that
                                         % we get a good estiamte for the
                                        % area of the peak
            for j =1: length(output_xyzVxyz)
                row_a= ismember(output_xyzVxyz{1,j}(:,4:end),compare(:,4:end),'rows'); %compare vel. of partciles to distinguish them compare all rows at same time
                compare= [compare;output_xyzVxyz{1,j}(~row_a,:)];  % add only the rows which were not yet present in compare to not count them twice
            end
            row_end = ismember(obj.xyzVxyz(:,4:end),compare(:,4:end),'rows');
            obj.Snapshot('end decc, detected', obj.in.M_time_vec(end-1),obj.xyzVxyz(row_end,:),[])
%% now there is fucntion we can call that does that for us plot_TOF_laser
%             if obj.in.verbose
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
        end

        function [h,area,TOF_cut]=gain_TOF(obj) % save max TOF as well as area to compare gfain factor to experiment
            h=max(obj.TOF_save(:,2));
            particle_entries=obj.TOF_save(:,2)>0;  % make array with logical one where particles are detected 0 otherwise
            TOF_cut=obj.TOF_save(particle_entries,:); % make a variable of cut TOF where entries are non zero also for compare plots
            [area,~]= size(obj.output{end-1,3});   % integrate area of TOF to compare to experiments and check performance        
        end

        function plot_TOF_laser(obj)
            plot(obj.TOF_save(:,1)*1e6+obj.in.params.FLY_incoupling_time*1e6,obj.TOF_save(:,2))
            xlabel('time (/mu s'); ylabel('# particles in laser volume');
            ylim([0,max(obj.TOF_save(:,2)+2)]);
            xlim([obj.TOF_save(1,1)*1e6+obj.in.params.FLY_incoupling_time*1e6,obj.TOF_save(end,1)*1e6+obj.in.params.FLY_incoupling_time*1e6]);
        end

        function saveWorkspace(obj)
%             arguments
%                 obj   double
%                 options.file_Name  string = dtaeTime('now')
%             end
            params = obj.in.params;  % create local varibale to safe otherwise matlab complains
            TOF_save_l = obj.TOF_save;
            output_l= obj.output;
            electrode_sequences = obj.in.electrode_sequences;
            filename = append('./output/',regexprep(datestr(datetime('now')),' ','_'),'_',num2str(obj.in.params.FLY_voltage_on_electrodes),'kV','_','FM',num2str(obj.in.params.FLY_focusing_mode_bool),'_',num2str(obj.in.params.FLY_target_velocity));
            fprintf(append('saving output, TOF_sa  ve, electrode sequences, in.params to',filename))
            save(filename,'TOF_save_l', 'output_l','electrode_sequences', 'params');
            
        end

    end
end


