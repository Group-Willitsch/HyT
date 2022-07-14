%% This script generate the default Matlab sequences for the experiment for
% a set of initial and final velocities

% the whole script is based on the boolean 'AutomFindFinalVelocity', which
% uses the bisecant optimizer coded in InputParamters.m -> generateMatlabTimeSequence

% It's made to be slow but safe: an atan function makes sure that the
% closer the phase to the final one, the smaller the step. It can take half
% a day to generate all the sequences, but it never fails.

%folder tree:
%  ./sequence/ 
% /focusing12p5kV    or /focusing13p5kV or /focusing10kV or whatever                       
clear all; close all; clc;

in = InputParameters(13.5, true, 450, 30, 'Phase', 40, 'Verbose', true, 'AutomFindFinalVelocity', false);
in.autom_find_final_vel = true; % set to true later
in.always_generate_M_seq = true;


%% initial and final velocities
start_vels = [415, 420, 430, 440, 450, 460, 470, 490, 500, 510];  % beam velocities
start_vels = [470, 490, 500, 510];  % beam velocities
start_vels = [425, 435, 445, 455, 465, 475, 480, 485, 495, 505];
start_vels = 415:5:510;
target_vels = [30, 49, 90, 202, 274, 328, 374, 414]; % target velocities
% target_vels = [92]; % target velocities

target_phases_FM_12p5 = [48.836, 48.7, 47.5, 41, 35, 26, 18, 10]; % target phases, FM, 12.5 kV from Nicolas table on OneNote
target_phases_FM_12p5 = [38, 37, 36, 28, 23, 20, 10, 5]; % target phases, FM, 12.5 kV from Nicolas table on OneNote
% target_phases_FM_12p5 = [30]; % target phases, FM, 12.5 kV from Nicolas table on OneNote
target_phases_FM_10 = [61.55, 61.05, 59.14, 49, 38.9, 29.3, 19.4, 10]; % target phases, FM, 10 kV from Nicolas table on OneNote
target_phases_FM_13p5 = []; % to be done

target_phases = target_phases_FM_12p5 ; % CHANGE HERE when you change voltage or mode

for i = 1:length(start_vels) % cycle on beam velocities
   for j = 1:length(target_vels) % cycle on target velocities
       %% set new velocities and generate the new sequence
        in.changeVelocities(start_vels(i), target_vels(j) , target_phases(j) )
        % (saving is called in InputParameters.m)
    end
end