function phase_changer_cos(in,a,b)
arguments
in
a  double
b  double 
end

if nargin < 2
a = 0;
end
% load fields loacally
M_trigger_pattern_old = in.M_trigger_pattern;
M_time_old = in.M_time_vec;
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

in.params.CALC_phase_distance = in.params.CALC_phase_distance;
fprintf('making new switching sequence using polynomial paramters given....');

switch_pos_NM = in.params.PHYS_distance_stages/2 + in.params.CALC_phase_distance: in.params.PHYS_distance_stages : in.params.PHYS_length_dec-in.params.PHYS_distance_stages/2 + in.params.CALC_phase_distance;


switch_pos_Foc_pulse = in.params.PHYS_distance_stages*3/2.0 -in.params.CALC_phase_distance: in.params.PHYS_distance_stages:in.params.PHYS_length_dec - in.params.PHYS_distance_stages/2 - in.params.CALC_phase_distance;
switch_pos_Foc_pulse = switch_pos_Foc_pulse - phase_distribution_cos(switch_pos_Foc_pulse); %here actuially not phase dist but phase dist change fucntion directl easier to deal with

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

        if  any(safer{i,3} <0)
            fprintf('synch. molecule bounced back get new parameters \n')
            in.M_synch_velocity = safer{i,3};
            return
        end

        if  safer{i,1}(end) >= 5e-3
            fprintf('synch moleucle travaled longer than 5 ms get new parameters \n ')
            in.M_synch_velocity = safer{i,3};
            return
        end

        opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-8, 'Events', @(t,x) EventsFcn_FM(t,x,i),InitialStep=2e-9);
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

in.M_synch_time = union(in.M_synch_time,t_end);
in.M_synch_position = union(in.M_synch_position,x_free);
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


%% change phase according to distribution defined here, we start with arctan dist (lower phase beginning, higher towards the end) but can be changed

    function phase_change = phase_distribution_cos(x)
        phase_change = a*cos(b*pi*(x-switch_pos_Foc_pulse(1))/in.params.PHYS_distance_stages);
    end % r lets the optimizer also adjust the strength of the change of the phase

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
