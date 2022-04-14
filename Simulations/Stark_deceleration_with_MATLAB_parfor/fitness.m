%% Fitness function for optimization parametrs used have to be set in file where we run CMAES
%% and phase_dist has to be set in optimizer function (maybe better to make different optimizers)
% 
% here paramter come in formn of matrix and we use parfor yo evaluatew
function fitness_val= fitness(parameter_matrix) %output needs to be a small number since searches minimum
% tic;
%set new M_time_vec by using optimizer
[~,M] = size(parameter_matrix); % get Population size (number new paramters)
fitness_val = zeros(1,M);
target_vel = 30;
num_part = 10e3;
% suggested action of matlab to avoid errors when calling
% parameter_matrix(1,i) something in parfor loop error is due to brodacast
% variable which is solved by using temporary fariables below

a_param = parameter_matrix(1,:);
b_param = parameter_matrix(2,:);
if M==1
    % Needs to be inside parfor loop such that all cores cqan see it
    in=InputParameters(12.5, true, 450, target_vel, 'Phase', 48.83, "FortranSeqBool", false, 'Verbose', false);
    in.num_particles = num_part;

    % Is called at each core indivuidally with the corresponding paramters
    % makes new in.M_time_vec for each core in order to run them
    optimizer_speed_test(in,a_param(1),b_param(1));

    % make check if bounced back or if vel to highs
    if any(in.M_synch_velocity < 0) || in.M_synch_velocity(end) > 40
        fitness_val(1,1) = NaN;
        return
    end
    % propagate particles using the sequence given by optimizer
    in.propagateParticles_euler();
%     fprintf('time to spropagate ')
%     toc;
    fitness_val(1,1) = Gaussian(in.output{4,3},target_vel);
%  
else

parfor i=1:M
    % Needs to be inside parfor loop such that all cores cqan see it
    in=InputParameters(12.5, true, 450, target_vel, 'Phase', 48.83, "FortranSeqBool", false, 'Verbose', false);
    in.num_particles = num_part;

    % Is called at each core indivuidally with the corresponding paramters
    % makes new in.M_time_vec for each core in order to run them
    optimizer_speed_test(in,a_param(i),b_param(i));

    % make check if bounced back or if vel to highs
    if any(in.M_synch_velocity < 0) || in.M_synch_velocity(end) > 40
        fitness_val(1,i) = NaN;
        continue
    end
    % propagate particles using the sequence given by optimizer
    in.propagateParticles_euler();
%     fprintf('time to spropagate ')
%     toc;
    fitness_val(1,i) = Gaussian(in.output{4,3},target_vel)
%  
end

end

fitness_val
end