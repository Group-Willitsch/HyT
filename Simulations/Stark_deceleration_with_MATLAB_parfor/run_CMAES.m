% SCRIPT TO RUN CMAES ALGORTIHM TO FIND MINIMUM OF OBJECTIVE FUNCTION
clearvars ; 
% cell variables fro the parameters and a streing array are used to be able
% to different simualtions using only an index i, one can define in the bash
% script

% define your intial guess of the paramters you want to optimize will be
% mean of normal distribuiton (Took standard sequence thus all zeros)
parameters = cell(5,1);
parameters{1,1} = [0;0];
parameters{2,1} = [0;0]; % refers to cos fucntion if you ran this it makes sense to use boundaries for the amplitude to avoid overlaps between pulses
parameters{3,1} = [0;0;0;0];
parameters{4,1} = [0;0;0;0;0];
parameters{5,1} = [0;0;0;0;0;0];

% define standard deviations for the normal dist. if you have a guess in what range the best solution should be set it to 1/3 of the intial search region
% these are rough guesses of 1/3 of search regions where there are
% solutions within that will give non NaN results
standard_deviation = cell(5,1);
standard_deviation{1,1} = [5;3];
standard_deviation{2,1} = [5;3]; % refers to cos fucntion if you ran this it makes sense to use boundaries for the amplitude to avoid overlaps between pulses
standard_deviation{3,1} = [5;5;5;5];
standard_deviation{4,1} = [5;5;5;5;5];
standard_deviation{5,1} = [5;5;5;5;5;5];

% fitness caller array used to define which function to use for
% changing the sequence, by calling the corresponding fitness function


% get # of cpus from bash script in order to create parpool of correct size
num_cpus_bash = str2num(getenv("SLURM_CPUS_PER_TASK"));
index_function = str2num(getenv("index_funct"));

if isempty(num_cpus_bash)
num_cpus_bash = 3; % to try out stuff locally define # cpu here
end
if isempty(index_function)
index_function =1; % to try out stuff locally define here which function to call in fitness to change the pahse along the decc.
% number refers to indexing of
% ["phase_changer_atan","phase_changer_cos","phase_changer_poly_4","phase_changer_poly_5","phase_changer_poly_6"]
% in fitness function, pass index to fitness.m in CMAES call
end

%% Useful optional arguments to adapt CMAES function to see all go beginning CMAES file

% GENERAL optional arguments
opts.PopSize = 2*num_cpus_bash; % number of offsrping created in each generation, shoul scale linearly with duration of 
% optimization higher numbers give better prob. of finding global minima
% (default (4 + floor(3*log(N))), also increases parent number mu
% opts.DispModulo = 100  % [0:Inf], disp messages after every i-th iteration, shows #function_evals,iterations_number,x_mean etc.;
opts.SaveVariables = 'off'; % saves variables to .mat file, I did a save command myself for more control (default: off)
%opts.Noise.on = 0; % For noisy fitness function (uncertainties) OPTS.Noise.on = 1. Interferes presumably with sigma termination criteria... (default: off)
opts.EvalParallel = 'on'; % fitness function is now given a matrix with all paramters of one generation in it, access like a = parameter_matrix(1,:),b = parameter_matrix(2,:)
opts.EvalInitialX = 'on'; % also intial guess is evaluated with first gneration togeter
% opts.LBounds = '-Inf % lower bounds, scalar or Nx1-vector'; for parmaeters
% opts.UBounds = 'Inf  % upper bounds, scalar or Nx1-vector'; for parameters


% STOP CRITERIA
opts.TolFun = 0.5;    % stops CMAES if results of fitness fucntion changes smaller than this value (default 1e-12)
opts.TolHistFun = 0.5;  % stop CMAES if max(f.hist)-min(f.hist) <= TolHistFun, history of fintess values is updated
% opts.TolX         = '1e-11*max(insigma) % stop if x-change smaller TolX'(default: 1e-11) ;
% opts.TolUpX       = '1e3*max(insigma) % stop if x-changes larger TolUpX' (default: 1e3);
% opts.StopFitness=0.01; % stop optimization if f(x_min) is smaller than this value (-INF default)
% opts.MaxFunEvals = 5;  % give max number of fitness function evaluationsuseful to test stuff loacly if e.g. result saving works etc (default INF)
% opts.StopIter =10; % give max number of iterations meaning how many generations are formed before algorithm stops (default: Inf)
% opts.DiffMaxChange = 'Inf  % maximal variable change(s), can be Nx1-vector';
% opts.DiffMinChange = '0    % minimal variable change(s), can be Nx1-vector';


delete(gcp('nocreate')); % shut down parallel pool if one exists
maxNumCompThreads(num_cpus_bash);
parpool(num_cpus_bash,IdleTimeout=600) % give parpool correct number cpus and also give long enough idle timeout such that parpool does not close when still running

% CALL THE CMAES FUNCTION attention here save_hist has ben added into CMAES
% lines 968-971 by Nicolas to save always the best fitness and afterwards be able to see how it evolved over time
[x_min, f_min, counterval,stopflag,out,bestever,save_hist] = cmaes('fitness_atan',parameters{index_function,1},standard_deviation{index_function,1},[opts],index_function);

% save results in a file
filename = append('./output/','atan_',regexprep(datestr(datetime('now')),' ','_'),'_','12.5','kV','_','FM','_','30');
save(filename,'x_min', 'f_min','bestever','stopflag', 'out','counterval',"save_hist");
delete(gcp('nocreate')); % shut down created parpool such that it doesnt run for 600 min!

%% Explenation input/output arguments given at beginning CMAES.m file
% Input arguments: 
%
%  FUN is a string function name like 'frosen'. FUN takes as argument
%     a column vector of size of X0 and returns a scalar. An easy way to
%     implement a hard non-linear constraint is to return NaN. Then,
%     this function evaluation is not counted and a newly sampled
%     point is tried immediately.
%
%   X0 is a column vector, or a matrix, or a string. If X0 is a matrix,
%     mean(X0, 2) is taken as initial point. If X0 is a string like
%     '2*rand(10,1)-1', the string is evaluated first.
%
%   SIGMA is a scalar, or a column vector of size(X0,1), or a string
%     that can be evaluated into one of these. SIGMA determines the
%     initial coordinate wise standard deviations for the search.
%     Setting SIGMA one third of the initial search region is
%     appropriate, e.g., the initial point in [0, 6]^10 and SIGMA=2
%     means cmaes('myfun', 3*rand(10,1), 2).  If SIGMA is missing and
%     size(X0,2) > 1, SIGMA is set to sqrt(var(X0')'). That is, X0 is
%     used as a sample for estimating initial mean and variance of the
%     search distribution.
%
%   OPTS (an optional argument) is a struct holding additional input
%     options. Valid field names and a short documentation can be
%     discovered by looking at the default options (type 'cmaes'
%     without arguments, see above). Empty or missing fields in OPTS
%     invoke the default value, i.e. OPTS needs not to have all valid
%     field names.  Capitalization does not matter and unambiguous
%     abbreviations can be used for the field names. If a string is
%     given where a numerical value is needed, the string is evaluated
%     by eval, where 'N' expands to the problem dimension
%     (==size(X0,1)) and 'popsize' to the population size. 

% Ouput:

% [XMIN, FMIN, COUNTEVAL, STOPFLAG, OUT, BESTEVER] = ...
%    CMAES(FITFUN, X0, SIGMA)
% returns the best (minimal) point XMIN (found in the last
% generation);

% fitness function value FMIN of fitness(XMIN);

% the number of needed function evaluations COUNTEVAL;

% STOPFLAG value as cell array,
% where possible entries are 'fitness', 'tolx', 'tolupx', 'tolfun',
% 'maxfunevals', 'maxiter', 'stoptoresume', 'manual',
% 'warnconditioncov', 'warnnoeffectcoord', 'warnnoeffectaxis',
% 'warnequalfunvals', 'warnequalfunvalhist', 'bug' (use
% e.g. any(strcmp(STOPFLAG, 'tolx')) or findstr(strcat(STOPFLAG,
% 'tolx')) for further processing);
% 
% record struct OUT with some
% more output, where the struct SOLUTIONS.BESTEVER contains the overall
% best evaluated point X with function value F evaluated at evaluation
% count EVALS. 

% The last output argument BESTEVER equals 
% OUT.SOLUTIONS.BESTEVER. Moreover a history of solutions and 
% parameters is written to files according to the Log-options. 

