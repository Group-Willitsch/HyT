function phase_changer_poly_5(in,param_vec)
arguments
in
param_vec
end
a = param_vec(1);
b = param_vec(2);
c = param_vec(3);
d = param_vec(4);
e = param_vec(5);







% save old version of sequence to be able to compare old and new version of
% sequence generation
M_trigger_pattern_old = in.M_trigger_pattern;
M_time_old = in.M_time_vec;
% define fields loacally for spped up
ax_norm_interpl= in.ax_norm_interpl;
ax_neg_interpl = in.ax_neg_interpl;
% if ~isempty(in.M_time_vec) || ~ismepty(in.M_synch_velocity) || ~ismepty(in.M_synch_position) || ~ismepty(in.M_synch_time)
%     in.M_time_vec = [];
%     in.M_synch_velocity = [];
%     in.M_synch_position = [];
%     in.M_trigger_pattern = [];
%     in.M_synch_time = [];
% 
% end




fprintf('making new switching sequence using polynomial paramters given....');
% Important to note here is that we can decide if we want to change
% switch_pos_NM and switch_pos_FM_pulse  in the same directions we either subtarct or add 
% the change in pahse distance in both cases and if we want to change them in different directions we make a sign change in one of the two places

% define switching postions where we switch off normal mode fields/turn on focusing pulses in the FM or where we switch from V to H and vice versa in NM
switch_pos_NM = in.params.PHYS_distance_stages/2 + in.params.CALC_phase_distance: in.params.PHYS_distance_stages : in.params.PHYS_length_dec-in.params.PHYS_distance_stages/2 + in.params.CALC_phase_distance;
switch_pos_NM = switch_pos_NM + phase_distribution_pol(switch_pos_NM)* in.params.PHYS_distance_stages / 180; % change the pahse distance using an analytical fucntion f(x_pos)

% define switching positions where we turn off the focusing pulses and turn on the normal mode fileds again (only used in focusing mode)
switch_pos_Foc_pulse = in.params.PHYS_distance_stages*3/2.0 -in.params.CALC_phase_distance: in.params.PHYS_distance_stages:in.params.PHYS_length_dec - in.params.PHYS_distance_stages/2 - in.params.CALC_phase_distance;
switch_pos_Foc_pulse = switch_pos_Foc_pulse - phase_distribution_pol(switch_pos_Foc_pulse)* in.params.PHYS_distance_stages / 180;

%combine both postions to get switching psotions of the focusing mode
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


% fucntion array containing fucntions for all field configurations
    dydt_array = {@(t,x) dydtNormVerticalOn(t,x),...
        @(t,x) dydtNegHorizontalOn(t,x),...
        @(t,x) dydtNormHorizontalOn(t,x),...
        @(t,x) dydtPosVerticalOn(t,x),...
        @(t,x) dydtNormVerticalOn(t,x),...
        @(t,x) dydtPosHorizontalOn(t,x),...
        @(t,x) dydtNormHorizontalOn(t,x),...
        @(t,x) dydtNegVerticalOn(t,x)};
    saver = cell (length(switch_pos_FM)+1, 3); % variable to save all times, psoitions and velcoites from each ode call in the for loop to make interpolation
    saver{1,1} = [0]; % position {i,1} is the time
    saver{1,2} = [0]; % {i,2} is the position
    saver{1,3} = [in.params.CALC_vel_synch_mol,in.params.CALC_vel_synch_mol]; %{i,3} is the velocity


    % loop 
    for i = 1 : 1 : length(switch_pos_FM)
        
        % check if velocity smaller than zero if so return but safe
        % velocities such that fintess sees them and also sees that vel was
        % negative and can give NaN for revalution of parameters
        if  any(saver{i,3} <0)
            fprintf('synch. molecule bounced back get new parameters \n')
            in.M_synch_velocity = saver{i,3};
            return
        end
        % Not sure why this was needed issue occured where molecule was flying longer than 5e-3 which gives error in ode but had not a negative velocity
        % same as above return it with corresponding velocites to fitness
        % and if v < 20 m/s fitness gives NaN thus reevaluation takes place
        if  saver{i,1}(end) >= 5e-3
            fprintf('synch moleucle travaled longer than 5 ms get new parameters \n ')
            in.M_synch_velocity = saver{i,3};
            return
        end
        % actual integration till switching point is reached which leads to
        % a stop of the ode and thus in the next ode call the next field configuration is used
        opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-8, 'Events', @(t,x) EventsFcn_FM(t,x,i),InitialStep=2e-9);
        [M_synch_time_temp, x_Vx_temp] = ode45( dydt_array{mod(i, 8)+8*(~mod(i, 8))}, [saver{i,1}(end),5e-3], [saver{i,2}(end); saver{i,3}(end)], opts);
        % safe time, vel and pos given by ode solver
        saver{i+1,1} = M_synch_time_temp;
        saver{i+1,2} = x_Vx_temp(:,1);
        saver{i+1,3} = x_Vx_temp(:,2);
    end
    % make free propagation untill detection point
    t_end = saver{end,1}(end) + (in.params.PHYS_length_dec + in.params.PHYS_exit_to_detection - saver{end,2}(end))/saver{end,3}(end);
    v_free = saver{end,3}(end);
    x_free = saver{end,2}(end) + v_free *(t_end-saver{end,1}(end));

else
    % same principle for the normal mode
    dydt_array = {@(t,x) dydtNormVerticalOn(t,x),...
        @(t,x) dydtNormHorizontalOn(t,x)};

    saver = cell (length(switch_pos_NM)+1, 3); % idea is to safe all times, psoitions and velcoites from each ode call in the for loop
    saver{1,1} = [0];
    saver{1,2} = [0];
    saver{1,3} = [in.params.CALC_vel_synch_mol,in.params.CALC_vel_synch_mol];


    for i = 1 : 1 : length(switch_pos_NM)


        opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-8, 'Events', @(t,x) EventsFcn_NM(t,x,i));
        [M_synch_time_temp, x_Vx_temp] = ode45( dydt_array{mod(i, 2)+2*(~mod(i, 2))}, [saver{i,1}(end), 5e-3], [saver{i,2}(end); saver{i,3}(end)], opts);
        saver{i+1,1} = M_synch_time_temp;
        saver{i+1,2} = x_Vx_temp(:,1);
        saver{i+1,3} = x_Vx_temp(:,2);

    end

    t_end = saver{end,1}(end) + (in.params.PHYS_length_dec + in.params.PHYS_exit_to_detection - saver{end,2}(end))/saver{end,3}(end);
    v_free = saver{end,3}(end);
    x_free = saver{end,2}(end) + v_free *(t_end-saver{end,1}(end));
    

end 
% empty ou the old variables
    in.M_synch_time = [];
    in.M_synch_position =[]; 
    in.M_synch_velocity =[];

    % fuse together all the times, pos, vel given by the ode solver above
for j =1 : length(saver)
    in.M_synch_time = union(in.M_synch_time,saver{j,1});
    in.M_synch_position = union(in.M_synch_position,saver{j,2}); 
    in.M_synch_velocity =cat(1,in.M_synch_velocity,saver{j,3}(2:end));
end

 % fuse them to
in.M_synch_time = union(in.M_synch_time,t_end);
in.M_synch_position = union(in.M_synch_position,x_free);
in.M_synch_velocity = cat(1,in.M_synch_velocity,v_free);

toc;
% end of the numerical integration


%again make check if synch molecule vel is smaller zero (here artefact of old code)
if saver{1,3}(end) < 0
    %     error("Sychronous molecule is bounced back, please lower the phase angle!") %
    return
end



t_x_interpl = griddedInterpolant(in.M_synch_position, in.M_synch_time); % inetrpolate pos and times of integration above to obatain the exact field switching time based on the switching position (pos synchronus molecule)
if in.params.FLY_focusing_mode_bool
    % define all positions to get the correct M_time_vec 0 and detection
    % point have to be added to the switch positions
    positions_to_jump = [0,switch_pos_FM,in.params.PHYS_length_dec + in.params.PHYS_exit_to_detection];

    in.M_time_vec = t_x_interpl(positions_to_jump);  % using the interpolate we get switching times of all positions of interest
    trigger_pattern = repmat(["b0011";"b0100";"b1100";"b0010";"b0011";"b1000";"b1100";"b0001"], in.params.PHYS_number_of_electrodes,1); % intialize the pattern of the FM mode
    in.M_trigger_pattern = trigger_pattern(1:length(in.M_time_vec)); % Make trigger pattern same length as M_time_vec
    in.M_trigger_pattern(end-1:end)= "b0000"; % set electrodes to zero for last entry (detection)
else
%     same for NM mode as above
positions_to_jump = [0,switch_pos_NM,in.params.PHYS_length_dec + in.params.PHYS_exit_to_detection];
   
    in.M_time_vec = t_x_interpl(positions_to_jump);
    trigger_pattern = repmat(["0x0010";"0x0020"], in.params.PHYS_number_of_electrodes,1);
    in.M_trigger_pattern = trigger_pattern(1:length(in.M_time_vec));
    in.M_trigger_pattern(end-1:end)= "0x0000";
end



fprintf("Precise final velocity of sequence is %s", num2str(in.M_synch_velocity(end)))
in.params.FLY_simulated_target_vel = round( in.M_synch_velocity(end), 0);
% stores in a quite useless but safe variable the effective final velocity

fprintf('\tFinal Matlab velocity is %d\n', in.params.FLY_simulated_target_vel)
simulated_target_vel = in.params.FLY_simulated_target_vel;
%             I return the compute velocity, to be used if needed.
%             for some reason I cannot return the class variable but I must
%             declare a local one. Probably due to handle and class stuff.

%% Uncomment to have a comparsion between the time sequences from the old version in Inputparamters and the one using the extended fields here

% pattern_old=char(M_trigger_pattern_old); % ugly code to turn M_pattern into array of doubles such that we can compare the pattern such that we get a string of ones and zeros
% pattern_old=double(string(pattern_old(:,end-3:end))); % the same length as M_pattern where 1 means the correspondign electrode is on or off at this specific time
% 
% 
% rod1_old = double(pattern_old== 1000 | pattern_old==1100)+1; % numbering of electrtodes follows b1100 sample where electrode corresponds to the electrode represented by the first number of b1100 etc
% rod2_old = double(pattern_old== 100 | pattern_old==1100)+3; % we add numbers to put the on same plot but still be able to see each individual sequence
% rod3_old = double(pattern_old== 10 | pattern_old==11)+5; % to get meaning of the numbers b1000=1000, b1100=1100, b0011= 11, b0010=10, b0001=1
% rod4_old = double(pattern_old== 1 | pattern_old==11)+7;
% 
% pattern=char(in.M_trigger_pattern); % ugly code to turn M_pattern into array of doubles such that we can compare the pattern such that we get a string of ones and zeros
% pattern=double(string(pattern(:,end-3:end))); % the same length as M_pattern where 1 means the correspondign electrode is on or off at this specific time
% 
% 
% rod1= double(pattern== 1000 | pattern==1100)+1; % numbering of electrtodes follows b1100 sample where electrode corresponds to the electrode represented by the first number of b1100 etc
% rod2= double(pattern== 100 | pattern==1100)+3; % we add numbers to put the on same plot but still be able to see each individual sequence
% rod3= double(pattern== 10 | pattern==11)+5; % to get meaning of the numbers b1000=1000, b1100=1100, b0011= 11, b0010=10, b0001=1
% rod4= double(pattern== 1 | pattern==11)+7;

%% Plot of time sequence using stairs fucntion to get the rigth pattern also used in the experiment
% figure()
% title('Trigger sequence electrodes (Channel i labframe/other one)')
% hold on
% stairs(in.M_time_vec*10^3,rod1,'k')
% stairs(in.M_time_vec*10^3,rod2,'k')
% stairs(in.M_time_vec*10^3,rod3,'k')
% stairs(in.M_time_vec*10^3,rod4,'k')
% 
% stairs(M_time_old*10^3,rod1_old,'r')
% stairs(M_time_old*10^3,rod2_old,'r')
% stairs(M_time_old*10^3,rod3_old,'r')
% stairs(M_time_old*10^3,rod4_old,'r')
% 
% ylim([0,9]); xlabel('time(ms)');
% legend('new','new','new','new','old','old', 'old', 'old');
% hold off


%% change phase according to distribution defined here

    function phase_change = phase_distribution_pol(x)
        x_0 = in.params.PHYS_length_dec/2;
        phase_change = a*(x.^5 -(x_0^5)/6) + b*(x.^4 -(x_0^4)/5) + c*(x.^3 -(x_0^3)/4) + d*(x.^2 -(x_0^2)/2) + e*(x-x_0/2); % divivde by r to set scale of fucntion back start r=100 needed to bring phase cvhange to -4,4 range
    end % r lets the optimizer also adjust the strength of the change of the phase


% Defintions of the functions, the function array refers to when th synch.
% molecule is proagated above (integration part above)
    function dxdt = dydtNormVerticalOn(t,x)
        dxdt = [x(2); ax_norm_interpl(x(1),0,0)];
       
    end

    function dxdt = dydtNormHorizontalOn(t,x)
         dxdt = [x(2); -ax_norm_interpl(in.params.PHYS_length_dec -x(1), 0, 0)];
    end

    function dxdt = dydtNegVerticalOn(t,x)
        dxdt = [x(2); ax_neg_interpl(x(1),0,0)];     
    end

    function dxdt = dydtNegHorizontalOn(t,x)
        dxdt = [x(2);  -ax_neg_interpl(in.params.PHYS_length_dec - x(1),0,0)]; 
       
    end

    function dxdt = dydtPosVerticalOn(t,x)
         dxdt = [x(2);ax_neg_interpl(x(1),0,0)]; 

    end
        
    function dxdt = dydtPosHorizontalOn(t,x)
         dxdt = [x(2); -ax_neg_interpl(in.params.PHYS_length_dec -x(1),0,0)]; 
        
    end
                                                             
                        



end