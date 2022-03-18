clear all; close all;
maxNumCompThreads(6);
set(0, 'DefaultLineLineWidth', 3);
myInput = InputParameters(12.5, true, 450, 30, 'Phase', 48.83, "FortranSeqBool", false);
% myInput = InputParameters(10, false, 450, 92);
% myInput.loadFortranTimeSequence2();
% myInput.loadAccelerationFields();


% myInput.InterpolateAccField();
% myInput.plotAccelerationFields();
% myInput.generateMatlabTimeSequence();
% myInput.params.CALC_phase_degrees = 10;
% myInput.params.CALC_save_sequence_bool = false;
% 
myInput.params.BEAM_avg_velocity_beam = 470; % Average velocity incoming package (m/s)
myInput.params.BEAM_radius_of_nozzle = 0.25e-3; % radius of the valve nozzle (m)
myInput.params.BEAM_long_pos_spread = 10e-3; % longitudinal position spread (m) - along x aixs or beam propagation
myInput.params.BEAM_long_vel_spread = 100; % relative velocity spread along beam axis 0.112 0.12
myInput.params.BEAM_trans_velocity_spread = 5; % velocity spread perp. to beam axis/average velocity 0.10

% t = InputParameters('./inputs/input_matlab.mat');
test = SetOfParticles(myInput);
test.num_particles = 100000;
test.createParticles();

tic; test.propagateParticles_euler(); toc;
% arrival_time_euler=test.arrival_time;
% xyzVxyz_euler=test.xyzVxyz;
% 
% tic; test.propagateParticles_ode45(); toc;
% arrival_time_ode45=test.arrival_time;
% xyzVxyz_ode45=test.xyzVxyz;
% 
% figure;scatter(xyzVxyz_ode45(:,1), xyzVxyz_ode45(:,4), 'filled');hold on;scatter(xyzVxyz_euler(:,1), xyzVxyz_euler(:,4), 'filled');legend("ode45","my euler")
% % figure;histogram(arrival_time_ode45,'BinEdges', 3.28e-3:1e-6:3.36e-3);hold on;histogram(arrival_time_euler, 'BinEdges', 3.28e-3:1e-6:3.36e-3);legend("ode45","my euler");
% figure;histogram(arrival_time_ode45,'BinEdges', 3.7e-3:1e-6:4.3e-3);hold on;histogram(arrival_time_euler, 'BinEdges', 3.7e-3:1e-6:4.3e-3);legend("ode45","my euler");
% test.plotTOF();
% tic; test.propagateParticlesAndSaveTrajectories(); toc;