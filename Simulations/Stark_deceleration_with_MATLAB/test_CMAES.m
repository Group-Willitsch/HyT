%% Apperantely there is constarint handling look at line 945 CMAES.m

clearvars ; 
in = InputParameters(12.5, true, 450, 30, 'Phase', 48.83, "FortranSeqBool", false, 'Verbose', false, 'AlwaysGenerateMSeq', false);
paramter_vector_0 = [3;0;0];
stand_dev = [0.5;8;3];
opts=cmaes;
% opts.StopFitness=0.01;
% opts.MaxFunEvals = 12; 
opts.PopSize = 30;
% opts.StopIter =10;
% opts.DispModulo = 1;
opts.SaveVariables = 'off';
% opts.TolX = 0.1; %stop if change sigma smaller than that
opts.Noise.on = 0;
opts.EvalParallel = 'yes';

[x_min, f_min, counterval, out, bestever] = cmaes('parallel_fit_test',paramter_vector_0,stand_dev,[opts],in);

filename = append('./output/',regexprep(datestr(datetime('now')),' ','_'),'_',num2str(obj.in.params.FLY_voltage_on_electrodes),'kV','_','FM',num2str(obj.in.params.FLY_focusing_mode_bool),'_',num2str(obj.in.params.FLY_target_velocity));
fprintf(append('saving x_min, f_min, bestever, out, counterval',filename))
save(filename,'x_min', 'f_min','bestever', 'out','counterval')