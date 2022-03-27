clear all; close all;

set(0, 'DefaultLineLineWidth', 3); % for plotting
in = InputParameters(12.5, true, 450, 30, 'Phase', 48.83, "FortranSeqBool", false, 'Verbose', false); 

p = SetOfParticles(in, 1000);
% p.createParticles();
% 
tic; p.propagateParticles_euler(); toc;
% arrival_time_euler=p.arrival_time;
% xyzVxyz_euler=p.xyzVxyz;
% 
% tic; p.propagateParticles_ode45(); toc;
% arrival_time_ode45=p.arrival_time;
% xyzVxyz_ode45=p.xyzVxyz;
% 
% figure;scatter(xyzVxyz_ode45(:,1), xyzVxyz_ode45(:,4), 'filled');hold on;scatter(xyzVxyz_euler(:,1), xyzVxyz_euler(:,4), 'filled');legend("ode45","my euler")
% % figure;histogram(arrival_time_ode45,'BinEdges', 3.28e-3:1e-6:3.36e-3);hold on;histogram(arrival_time_euler, 'BinEdges', 3.28e-3:1e-6:3.36e-3);legend("ode45","my euler");
% figure;histogram(arrival_time_ode45,'BinEdges', 3.7e-3:1e-6:4.3e-3);hold on;histogram(arrival_time_euler, 'BinEdges', 3.7e-3:1e-6:4.3e-3);legend("ode45","my euler");
% p.plotTOF();
% tic; p.propagateParticlesAndSaveTrajectories(); toc;