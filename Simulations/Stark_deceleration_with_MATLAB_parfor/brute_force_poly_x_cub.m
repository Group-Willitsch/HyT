
%% Brute Force script poly
in=InputParameters(12.5, true, 450, 30, 'Phase', 48.83, "FortranSeqBool", false, 'Verbose', false);
num_particles= 10e4;


a = -100:2:100;
b = -100:2:100;
A = length(a);
B = length(b);
fitness_safer = zeros(A,B);


M = zeros(A,B);
M2 = zeros(A,B);

for i =1:B
    for j =1:A
        M(i,j) = a(j);
        M2(:,j) = b;
    end
end


parfor i=1:A*B
    % Needs to be inside parfor loop such that all cores cqan see it
    in=InputParameters(12.5, true, 450, 30, 'Phase', 48.83, "FortranSeqBool", false, 'Verbose', false);
    in.num_particles = num_particles;

    % Is called at each core indivuidally with the corresponding paramters
    % makes new in.M_time_vec for each core in order to run them
    optimizer_poly_3(in,M(i),M2(i),c_temp);

    % make check if bounced back or if vel to highs
    if any(in.M_synch_velocity < 0) || in.M_synch_velocity(end) > 40 || in.M_synch_velocity(end) < 20
        fitness_safer(i) = 0.5;
        continue
    end
    % propagate particles using the sequence given by optimizer
    in.propagateParticles_euler();
    %     fprintf('time to spropagate ')
    %     toc;
    fitness_safer(i) = Gaussian(in.output{4,3},30)
    %
end

filename = append('results/','brute_force_x_cubic_2');
save(filename,'a','b','c','fitness_safer')

delete(gcp('nocreate'));
