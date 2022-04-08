%% Attention when using the fitness fucntion that gets a mtrix we need to give
%% fitness_Val as a row vector otherwise there will be problems with NaN statem

function fitness_val = parallel_fit_test(parameter_matrix,in)
%paramter vector is NxM where N is number of parameters
[~,M] = size(parameter_matrix); % get Population size (number new paramters)
fitness_val = zeros(1,M);
sigma_vel = [3,2,2];
for i=1:M

    optimizer(in,'r',parameter_matrix(1,i),'Amplitude',parameter_matrix(2,i),'x_tilde',parameter_matrix(3,i));

    if in.M_synch_velocity(end) < 0 || in.M_synch_velocity(end)>40 
        fitness_val(i) = NaN;
        continue
    end

    p = SetOfParticles(in, 10e4);
    p.createParticles();
    p.propagateParticles_euler()

    [size_check,~] = size(p.output{4,3}(2:end,4));
    if size_check == 0
        fitness_val(i) = NaN;
        continue
    end

    fitness_val(i) = Gaussian(p.output{4,3},sigma_vel,in.params.FLY_target_velocity);

end
end