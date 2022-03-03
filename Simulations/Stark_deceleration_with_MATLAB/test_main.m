clear all; close all;
set(0, 'DefaultLineLineWidth', 3);
myInput = InputParameters(12.5,true, 490, 48, 'Phase', 57.72);
% myInput = InputParameters(12.5,true, 450, 92, 'Phase', 47);
% myInput.loadFortranTimeSequence2();
% myInput.loadAccelerationFields();


% myInput.InterpolateAccField();
% myInput.plotAccelerationFields();
% myInput.generateMatlabTimeSequence();
% myInput.params.CALC_phase_degrees = 10;
% myInput.params.CALC_save_sequence_bool = false;
% 
% myInput.params.BEAM_avg_velocity_beam = 470; % Average velocity incoming package (m/s)
% myInput.params.BEAM_radius_of_nozzle = 0.25e-3; % radius of the valve nozzle (m)
% myInput.params.BEAM_long_pos_spread = 10e-3; % longitudinal position spread (m) - along x aixs or beam propagation
% myInput.params.BEAM_long_vel_spread = 100; % relative velocity spread along beam axis 0.112 0.12
% myInput.params.BEAM_trans_velocity_spread = 5; % velocity spread perp. to beam axis/average velocity 0.10
% 
% % t = InputParameters('./inputs/input_matlab.mat');
% test = SetOfParticles(myInput);
% test.num_particles = 100000;
% test.createParticles();
% tic; test.propagateParticles(); toc;
% % tic; test.propagateParticlesAndSaveTrajectories(); toc;