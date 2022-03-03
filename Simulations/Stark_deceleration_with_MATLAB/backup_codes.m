%10/2/2022
%parfor?

opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
figure;

parfor i = 2:1:obj.params.FLY_n_molecules
    [t, y] = ode45(@(t,y) obj.dydt_beam_norm(t,y,2), [0, 5e-3], reshape([obj.xyzVxyz_0(1,:);obj.xyzVxyz_0(i,:)], 12, 1), opts);
    size(y)
    plot(t,y(:,1));
end


%    08/02/2022
        function generateMatlabTimeSequence(obj)
            fprintf('Generating Matlab time sequence...');
            obj.params.CALC_phase_distance = obj.params.CALC_phase_degrees*obj.params.PHYS_distance_stages/180;
            obj.ax_norm_interpl = griddedInterpolant(linspace(0,obj.params.PHYS_distance_stages,111), obj.ax_norm(21:131,21,21));
            
            % start integration
            opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
            [t, x] = ode45(@(t,x) obj.dxdt(t,x), [0,5e-3], [0; obj.params.CALC_vel_synch_mol],opts);
            
            figure;
            plot(t,x(:,2));
            title('Velocity vs time');
            
            t_x_interpl = griddedInterpolant(x(:,1),t);
            positions_to_jump = [0, obj.params.PHYS_distance_stages/2 + obj.params.CALC_phase_distance: obj.params.PHYS_distance_stages: obj.params.PHYS_length_dec-(obj.params.PHYS_distance_stages/2 - obj.params.CALC_phase_distance)];
            obj.M_time_vec = t_x_interpl(positions_to_jump);
            obj.M_trigger_pattern = repmat(["0x0010","0x0020"], 1, round((length(obj.M_time_vec)-1)/2));
            obj.M_trigger_pattern(end)= "0x0000";
            
            if obj.params.CALC_save_sequence_bool
                seq_file=fopen(sprintf('./sequences/dec_%d_%d.out', obj.params.CALC_vel_synch_mol, round(x(end,2))),'w');
                fprintf(seq_file, '%8s\t%s\n',[string(round(obj.M_time_vec*1e9)); obj.M_trigger_pattern]);
                fclose(seq_file);
            end
            
            fprintf('done\n')
        end
        
        function dxdt = dxdt(obj,t,x)                    
            if x(2) < 0
                fprintf('ERROR: Velocity less than zero in Stark decelerator');
                return
            end
            
            if x(1) > obj.params.PHYS_length_dec-(obj.params.PHYS_distance_stages/2 - obj.params.CALC_phase_distance)
                dxdt = [x(2);0];
            else
                pos = mod(x(1),obj.params.PHYS_distance_stages);
                if pos <= obj.params.PHYS_distance_stages/2 + obj.params.CALC_phase_distance
                    dxdt = [x(2); obj.ax_norm_interpl(pos)];
                else
                    dxdt = [x(2); -obj.ax_norm_interpl(obj.params.PHYS_distance_stages-pos)];
                end
            end
        end



% 07/02/2022

function generateMatlabTimeSequence(obj)
        fprintf('Generating Matlab time sequence...')
        sequence_file=fopen('T2jump.dat','w');
        x = 0;
        vx = obj.params.CALC_vel_synch_mol;
        tof = 0;
        t_step = 5e-9;
        obj.ax_norm_interp = griddedInterpolant(linspace(0,obj.params.PHYS_distance_stages,111), obj.ax_norm(21:131,21,21));
        for i = 1:1:obj.params.PHYS_number_of_switching_times
            i
            fprintf(sequence_file, '%.2f\t%d\n', tof*1e9,mod(i,2));
            while x-i*obj.params.PHYS_distance_stages < obj.params.CALC_phase_in_distance
                [t, y] = ode45(@(t,y) obj.dydt_seq(t,y,i), [tof, tof+t_step], [x,vx]);
                tof = t(end);
                x = y(end,1);                
                vx = y(end,2);
                if vx <= 0
                    fprintf('ERROR: Velocity less than zero in Stark decelerator');
                    break;
                end
            end
        end
        fprintf(sequence_file, '%.2f\t0000\n', tof*1e9);
        tof = tof+(obj.params.PHYS_exit_to_detection-obj.params.CALC_phase_in_distance)/vx;
        fprintf(sequence_file, '%.2f\t0000\n', tof*1e9,mod(i,2));
        fprintf('\tdone\n')
        end
        
        function dydt = dydt_seq(obj,t,y,i)
            
            if ~mod(i,2)
                y(1)=y(1)+obj.params.PHYS_distance_stages;
            end
            if mod(y(1),obj.params.PHYS_distance_stages*2) < obj.params.PHYS_distance_stages
                dydt = [y(2); obj.ax_norm_interp(mod(y(1),obj.params.PHYS_distance_stages*2))];
            else
                dydt = [y(2); -obj.ax_norm_interp(obj.params.PHYS_distance_stages*2-mod(y(1),obj.params.PHYS_distance_stages*2))];
            end
        end
        
        
        function constructAccelerationFields(obj)
            %This function constructs position dependent acceleration fields for
            %synchronous molecules
            ax_norm_interp = griddedInterpolant(linspace(0,obj.params.PHYS_distance_stages,111), obj.ax_norm(21:131,21,21));
            obj.ax_norm_constructed = repmat([ax_norm_interp(0:1e-5:obj.params.PHYS_distance_stages/2 + obj.params.CALC_phase_distance),...
                -ax_norm_interp(obj.params.PHYS_distance_stages/2-obj.params.CALC_phase_distance:-1e-5:0)],1,obj.params.PHYS_number_of_electrodes);
%             plot(linspace(0,5.5*123,67773), obj.ax_norm_constructed)
        end