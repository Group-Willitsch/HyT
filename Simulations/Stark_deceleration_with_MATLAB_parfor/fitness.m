%% Fitness function for optimization parametrs used have to be set in file where we run CMAES
%% and phase_dist has to be set in optimizer function (maybe better to make different optimizers)
% 
% here paramter come in formn of matrix and we use parfor yo evaluatew
function fitness_val= fitness_atan(parameter_matrix,index_phase_changer) %output needs to be a small number since searches minimum
% tic;
%set new M_time_vec by using optimizer
[~,M] = size(parameter_matrix); % get Population size (number new paramters)
fitness_val = zeros(1,M);

% DEFINE SIMULATION PARAMETERS
target_vel = 30; % Define target velocity here
Phase_val = 48.83; % Define pahse corresponding to target velocity
num_part = 50e3; % Define number particles for simulation



% Define strings fro all phase_changer functions
pahse_changer_caller =["phase_changer_atan","phase_changer_cos","phase_changer_poly_4","phase_changer_poly_5","phase_changer_poly_6"];
% Safe the one we want to run in order to avoid brodcast variables in parfor loop
pahse_changer_used = pahse_changer_caller(index_phase_changer);

if M==1
    % Needs to be inside parfor loop such that all cores cqan see it
    in=InputParameters(12.5, true, 450, target_vel, 'Phase',Phase_val, "FortranSeqBool", false, 'Verbose', false);
    in.num_particles = num_part;

    % Is called at each core indivuidally with the corresponding paramters
    % makes new in.M_time_vec for each core in order to run them there
    feval(pahse_changer_used ,in,parameter_matrix(:,1))

    % make check if bounced back or if vel to high
    if any(in.M_synch_velocity < 0) || in.M_synch_velocity(end) > target_vel + 10  || in.M_synch_velocity(end) < target_vel -10
        fitness_val(1,1) = NaN;
        return
    end
    % propagate particles using the sequence given by optimizer
    in.propagateParticles_euler();

    fitness_val(1,1) = fitness_value_calculator(in.output{4,3},target_vel);
    [n_entrance,~]= size(in.output{1,3}); %take number of particles at entrance decc
    fitness_val(1,1) = 1e4*fitness_val(1,1)/n_entrance; %rescale fitness to make it more independent on the number of particles
%  
else

parfor i=1:M
    % Needs to be inside parfor loop such that all cores cqan see it
    in=InputParameters(12.5, true, 450, target_vel, 'Phase', Phase_val, "FortranSeqBool", false, 'Verbose', false);
    in.num_particles = num_part;

    % Is called at each core indivuidally with the corresponding paramters
    % makes new in.M_time_vec for each core in order to run them
    feval(pahse_changer_used,in,parameter_matrix(:,i))

    % make check if bounced back or if vel to highs
    if any(in.M_synch_velocity < 0) || in.M_synch_velocity(end) > target_vel + 10  || in.M_synch_velocity(end) < target_vel -10
        fitness_val(1,i) = NaN;
        continue
    end
    % propagate particles using the sequence given by optimizer
    in.propagateParticles_euler();

    fitness_val(1,i) = fitness_value_calculator(in.output{4, 3}, target_vel) % compute fitness, using particles at detection point that made it into laser volume
    [n_entrance, ~] = size( in.output{1,3} ); % take number of particles at entrance decc
    fitness_val(1,i) = 1e4*fitness_val(1,i)/n_entrance; %rescale fitness to make it more independent on the number of particles we run the sim with
%  
end

end
% Display some results to see what CMAES is doing
fitness_val
[min_val,index] = min(fitness_val);
num2str(min_val)
min_till_now = num2str(parameter_matrix(:,index))
end