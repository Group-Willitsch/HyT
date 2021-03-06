function optimizer_speed_test(in,options)
arguments
in
options.Amplitude   double =0
options.x_tilde    double =0
options.r          double =3
options.a  double =0
options.b  double =0
options.c  double =0
options.d  double =0
end

% load fields loacally

ax_norm_interpl= in.ax_norm_interpl;
ax_norm_H_interpl = in.ax_norm_H_interpl;
ax_neg_interpl = in.ax_neg_interpl;
ax_neg_H_interpl = in.ax_neg_H_interpl;
ax_pos_interpl = in.ax_pos_interpl;
ax_pos_H_interpl = in.ax_pos_H_interpl;




fprintf('making new M_time_vec and trigger pattern...');

switch_pos_NM = in.params.PHYS_distance_stages/2 + in.params.CALC_phase_distance: in.params.PHYS_distance_stages : in.params.PHYS_length_dec-in.params.PHYS_distance_stages/2 + in.params.CALC_phase_distance;
switch_pos_NM = switch_pos_NM + phase_distribution_atan(switch_pos_NM)* in.params.PHYS_distance_stages / 180;

switch_pos_Foc_pulse = in.params.PHYS_distance_stages*options.r/2.0 -in.params.CALC_phase_distance: in.params.PHYS_distance_stages:in.params.PHYS_length_dec - in.params.PHYS_distance_stages/2 - in.params.CALC_phase_distance;
switch_pos_Foc_pulse = switch_pos_Foc_pulse + phase_distribution_atan(switch_pos_Foc_pulse)* in.params.PHYS_distance_stages / 180;

switch_pos_FM = sort([switch_pos_NM,switch_pos_Foc_pulse]);


    function [value, isterminal, direction] = EventsFcn_NM(t, x,i) % This event function stops the ode solver once the molecule arrives at a switching point
        value = x(1) < switch_pos_NM(i);
        isterminal = 1;
        direction = 0;
    end

    function [value, isterminal, direction] = EventsFcn_FM(t, x,i) % This event function stops the ode solver once the molecule arrives at a switching point
        value = x(1) < switch_pos_FM(i) ;
        isterminal = 1;
        direction = 0;
    
    end
tic;
if in.params.FLY_focusing_mode_bool



    dydt_array = {@(t,x) dydtNormVerticalOn(t,x),...
        @(t,x) dydtNegHorizontalOn(t,x),...
        @(t,x) dydtNormHorizontalOn(t,x),...
        @(t,x) dydtPosVerticalOn(t,x),...
        @(t,x) dydtNormVerticalOn(t,x),...
        @(t,x) dydtPosHorizontalOn(t,x),...
        @(t,x) dydtNormHorizontalOn(t,x),...
        @(t,x) dydtNegVerticalOn(t,x)};
    safer = cell (length(switch_pos_FM)+1, 3); % idea is to safe all times, psoitions and velcoites from each ode call in the for loop
    safer{1,1} = [0];
    safer{1,2} = [0];
    safer{1,3} = [in.params.CALC_vel_synch_mol,in.params.CALC_vel_synch_mol];

    for i = 1 : 1 : length(switch_pos_FM)

        opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-8, 'Events', @(t,x) EventsFcn_FM(t,x,i));
        [M_synch_time_temp, x_Vx_temp] = ode45( dydt_array{mod(i, 8)+8*(~mod(i, 8))}, [safer{i,1}(end),5e-3], [safer{i,2}(end); safer{i,3}(end)], opts);
        safer{i+1,1} = M_synch_time_temp;
        safer{i+1,2} = x_Vx_temp(:,1);
        safer{i+1,3} = x_Vx_temp(:,2);


    end

    t_end = safer{end,1}(end) + (in.params.PHYS_length_dec + in.params.PHYS_exit_to_detection - safer{end,2}(end))/safer{end,3}(end);
    v_free = safer{end,3}(end);
    x_free = safer{end,2}(end) + v_free *(t_end-safer{end,1}(end));

else
    
    dydt_array = {@(t,x) dydtNormVerticalOn(t,x),...
        @(t,x) dydtNormHorizontalOn(t,x)};

    safer = cell (length(switch_pos_NM)+1, 3); % idea is to safe all times, psoitions and velcoites from each ode call in the for loop
    safer{1,1} = [0];
    safer{1,2} = [0];
    safer{1,3} = [in.params.CALC_vel_synch_mol,in.params.CALC_vel_synch_mol];


    for i = 1 : 1 : length(switch_pos_NM)


        opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-8, 'Events', @(t,x) EventsFcn_NM(t,x,i));
        [M_synch_time_temp, x_Vx_temp] = ode45( dydt_array{mod(i, 2)+2*(~mod(i, 2))}, [safer{i,1}(end), 5e-3], [safer{i,2}(end); safer{i,3}(end)], opts);
        safer{i+1,1} = M_synch_time_temp;
        safer{i+1,2} = x_Vx_temp(:,1);
        safer{i+1,3} = x_Vx_temp(:,2);

    end

    t_end = safer{end,1}(end) + (in.params.PHYS_length_dec + in.params.PHYS_exit_to_detection - safer{end,2}(end))/safer{end,3}(end);
    v_free = safer{end,3}(end);
    x_free = safer{end,2}(end) + v_free *(t_end-safer{end,1}(end));
    

end
    in.M_synch_time = [];
    in.M_synch_position =[]; 
    in.M_synch_velocity =[];

for j =1 : length(safer)

    in.M_synch_time = union(in.M_synch_time,safer{j,1});
    in.M_synch_position = union(in.M_synch_position,safer{j,2}); 
    in.M_synch_velocity =cat(1,in.M_synch_velocity,safer{j,3}(2:end));
end

in.M_synch_time = cat(1,in.M_synch_time,t_end);
in.M_synch_position = cat(1,in.M_synch_position,x_free);
in.M_synch_velocity = cat(1,in.M_synch_velocity,v_free);

toc;
% end of the numerical integration


% This error message is skipped and isnteadt we give NaN in return to the fitness fucntion
if safer{1,3}(end) < 0
    %     error("Sychronous molecule is bounced back, please lower the phase angle!") %
    return
end



t_x_interpl = griddedInterpolant(in.M_synch_position, in.M_synch_time); % inetrpolate pos and times of integration above to obatain the exact field switching time based on the switching position (pos synchronus molecule)
if in.params.FLY_focusing_mode_bool

    positions_to_jump = [0,switch_pos_FM,in.params.PHYS_length_dec + in.params.PHYS_exit_to_detection];

    in.M_time_vec = t_x_interpl(positions_to_jump);  %using the interpolate we get switcvhing times at pos of interest
    trigger_pattern = repmat(["b0011";"b0100";"b1100";"b0010";"b0011";"b1000";"b1100";"b0001"], in.params.PHYS_number_of_electrodes,1); %sequence of focusing mode is what is in [], % B = repmat(A,n) returns an array containing n copies of A in the row and column dimensions. The size of B is size(A)*n when A is a matrix.
    in.M_trigger_pattern = trigger_pattern(1:length(in.M_time_vec));
    in.M_trigger_pattern(end-1:end)= "b0000";
else
    %in NM we only add the phase dist. thus can be done like below
positions_to_jump = [0,switch_pos_NM,in.params.PHYS_length_dec + in.params.PHYS_exit_to_detection];
   
    in.M_time_vec = t_x_interpl(positions_to_jump);
    trigger_pattern = repmat(["0x0010";"0x0020"], in.params.PHYS_number_of_electrodes,1);
    in.M_trigger_pattern = trigger_pattern(1:length(in.M_time_vec));
    in.M_trigger_pattern(end-1:end)= "0x0000";
end



fprintf("Precise final velocity of sequence is %s", num2str(in.M_synch_velocity(end)))
in.params.FLY_simulated_target_vel = round( in.M_synch_velocity(end), 0);
% stores in a quite useless but safe variable the effective final velocity
% Will be used to create the filename

%           Uncomment to see the distribution of the phase change along the
%           deccc. usefull to set starting values
%             figure()
%             plot(positions_to_jump(2:end-1),Phase_changes_dist)


fprintf('\tFinal Matlab velocity is %d\n', in.params.FLY_simulated_target_vel)
simulated_target_vel = in.params.FLY_simulated_target_vel;
%             I return the compute velocity, to be used if needed.
%             for some reason I cannot return the class variable but I must
%             declare a local one. Probably due to handle and class stuff.


%% change phase according to distribution defined here, we start with arctan dist (lower phase beginning, higher towards the end) but can be changed
    function phase_change = phase_distribution_atan(x)
        phase_change = options.Amplitude*atan((x-in.params.PHYS_length_dec/2)/options.x_tilde);
    end
    function phase_change = phase_distribution_pol(x)
        x = x - in.params.PHYS_length_dec/2;
        phase_change = options.a + options.b*x + options.c*x.^2 + options.d*x.^3;
    end
    function phase_change = phase_distribution_sin(x)
        phase_change = options.Amplitude*sin((x - in.params.PHYS_length_dec/2)/options.x_tilde);
    end


    function dxdt = dydtNormVerticalOn(t,x)
        dxdt = [x(2); ax_norm_interpl(x(1),0,0)];
       
    end

    function dxdt = dydtNormHorizontalOn(t,x)
         dxdt = [x(2); ax_norm_H_interpl(x(1),0,0)];
    end

    function dxdt = dydtNegVerticalOn(t,x)
        dxdt = [x(2); ax_neg_interpl(x(1),0,0)];     
    end

    function dxdt = dydtNegHorizontalOn(t,x)
        dxdt = [x(2); ax_neg_H_interpl(x(1),0,0)]; 
       
    end

    function dxdt = dydtPosVerticalOn(t,x)
         dxdt = [x(2); ax_pos_interpl(x(1),0,0)]; 

    end
        
    function dxdt = dydtPosHorizontalOn(t,x)
         dxdt = [x(2); ax_pos_H_interpl(x(1),0,0)]; 
        
    end
                                                             
                        



end
