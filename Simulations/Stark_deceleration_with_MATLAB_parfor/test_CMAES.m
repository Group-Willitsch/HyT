
clearvars ; 
paramter_vector_0 = [0;0];
stand_dev = [10;10];
% num_cpus_bash = string2num(getenv('SLURM_CPUS_PER_TASK'))
num_cpus_bash = 4;
if isempty(num_cpus_bash)
num_cpus_bash = 24;
end
opts=cmaes;
% opts.StopFitness=0.01;
% opts.MaxFunEvals = 3;
opts.PopSize = 2*num_cpus_bash;
% opts.StopIter =10;
opts.DispModulo = 0;
opts.SaveVariables = 'off';
% opts.TolX = 0.1; %stop if change sigma smaller than that
opts.Noise.on = 0;
opts.EvalParallel = 'on';
% opts.EvalInitialX = 'off';
tic;
delete(gcp('nocreate')); % shut down parallel pool if one exists
maxNumCompThreads(num_cpus_bash);
parpool(num_cpus_bash) 
[x_min, f_min, counterval, out, bestever] = cmaes('fitness',paramter_vector_0,stand_dev,[opts]);
t_end = toc;
filename = append('./output/',regexprep(datestr(datetime('now')),' ','_'),'_','12.5','kV','_','FM','true','_','30');
fprintf(append('saving x_min, f_min, bestever, out, counterval',filename))
save(filename,'x_min', 'f_min','bestever', 'out','counterval','t_end');
