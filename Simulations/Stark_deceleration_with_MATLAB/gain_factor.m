clear all; closmyInput. ll;
% script that can be used to get the area TOF, height peak TOF as well as
% corresponding gain factors of NM and FM at different voltages and final
% velocities
%The only thing one needs to chnage to look at more or different points are
%the final velocities and the corresponding phases below for each case
%(NM/FM ,10/12.5 kV)

%% 12.5 kV NM not relevant since same results as 10 kV NM to check just uncomment the 12.5 kV NM parts

maxNumCompThreads(6);
vel_10 = [30,49,90,202,274,328];
vel_12 = [30,49,92,202,274,328];
phase_10_NM = [62.5,62,60.1,50,39.8,30];
phase_10_FM = [61.59,61.05,59.14,49,38.9,29.3];
% phase_12_NM  = [48.83,47.03,39.4,31.5,23.6,15.5];
phase_12_FM = [48.83,48.48,47.03,39.4,31.4,23.60,15.5];

[v_10_FM,v_10_NM,] = deal(zeros(size(vel_10)));
[v_12_FM,v_12_NM,] = deal(zeros(size(vel_12)));

safer_10_FM = cell(length(vel_10),4); % cell array to save all the neseccary data for plotting intialized such that correct length
safer_10_NM = cell(length(vel_10),4); % we safe a flag which mode at which voltage as well as height, area and area of TOF and a cut TOF
safer_12_FM = cell(length(vel_12),4);
% safer_12_NM = cell(length(vel_12),4);

safe_time = cell(2*(length(vel_10)+length(vel_12)),3); %safe the time we need for each step to see how fst code is when
%intialize everzthing
myInput = InputParameters(10,true, 450, 30, 'Phase',61.59, "FortranSeqBool", false,'Verbose',false,'AlwaysGenerateMSeq',true);

% % arrival_time_euler=test.arrival_time;
% xyzVxyz_euler=test.xyzVxyz;
% [h,a,TOF]=test.gain_TOF;
% TOF_save_eul = test.TOF_save;
% output_eul=test.output;
% electrode_sequences=test.myInput.electrode_sequences;_velocity_beam; % relative velocity spread along beam axis 0.112 0.12
myInput.params.BEAM_trans_velocity_spread = 5; % velocity spread perp. to beam axis/average velocity 0.10

compare = SetOfParticles(myInput);  
compare.num_particles = 10e4;                             % below loops over all different configurations and propagate particles
compare.createParticles();                                 % always save the height TOF, area TOF and cut TOF in corresponding  cell array
safe_output = {};
                                                          % also safe times of each step in cell array labeled and also save end velocity
tic;compare.propagateParticles_euler();tEnd=toc;
[h,a,TOF] = compare.gain_TOF;                               
safer_10_FM(1,:) = {'10 kV FM',h,a,TOF};
safe_time(1,:) = { 'FM 10 kV',tEnd, 30};



%% FM 10 kV
for i = 2 : length(vel_10)                                  
fprintf('load new velocitiy FM \n')
myInput.changeVelocities(450, vel_10(i), phase_10_FM(i))  
compare = SetOfParticles(myInput);
compare.num_particles = 10e4;   
compare.createParticles();
tic;compare.propagateParticles_euler();tEnd=toc;
[h,a,TOF] = compare.gain_TOF;
safer_10_FM(i,:) = {'10 kV FM',h,a,TOF};
safe_time(i,:) = {'FM 10 kV',tEnd,vel_10(i)};

end


 
%% FM 12.5 kV
fprintf('loading 12.5 kV FM mode \n')
myInput.changeFieldConfig(12.5,true) 
for i = 1:length(vel_12)
myInput.changeVelocities(450,vel_12(i),phase_12_FM(i))
compare = SetOfParticles(myInput);
compare.num_particles = 10e4;   
compare.createParticles();
tic;compare.propagateParticles_euler();tEnd=toc;
[h,a,TOF] = compare.gain_TOF;
safer_12_FM(i,:) = {'12.5 kV FM',h,a,TOF};
safe_time(i+length(vel_10),:) = {'FM 12 kV',tEnd,vel_12(i)};

end

%% NM 10 kV
fprintf('loading 10 kV NM mode \n')
myInput.changeFieldConfig(10,false)
for i = 1:length(vel_10)
fprintf('load new velocity NM \n')
myInput.changeVelocities(450, vel_10(i), phase_10_NM(i))
compare = SetOfParticles(myInput);
compare.num_particles = 10e4;   
compare.createParticles();
tic;compare.propagateParticles_euler();tEnd=toc;
[h,a,TOF] = compare.gain_TOF;
safer_10_NM(i,:) = {'10 kV NM',h,a,TOF};
safe_time(i+length(vel_10)+length(vel_12),:) = {'NM 10 kV',tEnd,vel_10(i)};

end

%% NM 12.5 kV
% fprintf('loading 12.5 kV NM mode \n')
% myInput.changeFieldConfig(10,false)
% for i=1:length(vel_12)
% myInput.changeVelocities(450,vel_12(i),phase_12_NM(i))
% compare = SetOfParticles(myInput);
% compare.createParticles();
% tic;compare.propagateParticles_euler();tEnd=toc;
% [h,a,TOF]=compare.gain_TOF;
% safer_12_NM(i,:)={'12.5 kV NM',h,a,TOF};
% safe_time(i+2*length(vel_10)+length(vel_12),:)={'NM 12.5 kV',tEnd, vel_12(i)};
% 
% end

a_10_NM = [safer_10_NM{:,3}];       % save all variables locally needed for plotting
a_10_FM = [safer_10_FM{:,3}];       % could also jsut use right side in plot fucntion
a_12_FM = [safer_12_FM{:,3}];

h_10_NM = [safer_10_NM{:,2}];
h_10_FM =[safer_10_FM{:,2}];
h_12_FM=[safer_12_FM{:,2}];

%% NM 12.5 kV part
% a_12_NM=[safer_12_NM{:,3}];
% h_12_NM=[safer_12_NM{:,2}];


%% Plot cut TOF for all cases
figure()
hold on
for i=1:6
    plot(safer_10_NM{i,4}(:,1)*1e6,safer_10_NM{i,4}(:,2),'Color','#0072BD')
    plot(safer_10_FM{i,4}(:,1)*1e6,safer_10_FM{i,4}(:,2),'r')
    plot(safer_12_FM{i,4}(:,1)*1e6,safer_12_FM{i,4}(:,2),'Color','#A2142F')
end
xlabel('time ( /mu s)'); ylabel('detected molecules'); legend('10kV_{NM}','10kV_{FM}','12.5kV_{FM}');
hold off

%% Plot the area of the TOF as fucntion of final vel.
figure('Name','area TOF')
hold on
plot(vel_10,a_10_NM,'-o')
plot(vel_10,a_10_FM,'-o')
% plot(vel_12,a_12_NM,'-o')
plot(vel_12,a_12_FM,'-o')
xlabel('final velocity (m/s)'); ylabel('area TOF-profile');legend('10kV_{NM}','10kV_{FM}','12.5kV_{FM}');
hold off

%% Plot gain factor area as fucntion of final vel.
 figure('Name','gain factor area')
hold on
plot(vel_10,a_10_FM./a_10_NM,'-o')
plot(vel_12,a_12_FM./a_10_NM,'-o')
xlabel('final velocity (m/s)'); ylabel('gain factor area');legend('10kV_{FM}','12.5kV_{FM}');
hold off

%% Plot the height of the TOF peak as fucntion of final vel.
figure('Name','height peak TOF')
hold on
plot(vel_10,h_10_NM,'-o')
plot(vel_10,h_10_FM,'-o')
% plot(vel_12,h_12_NM,'-o')
plot(vel_12,h_12_FM,'-o')
xlabel('final velocity (m/s)'); ylabel('peak TOF-profile');legend('10kV_{NM}','10kV_{FM}','12.5kV_{FM}');
hold off

%% Plot gain factor of the TOF peak as fucntion of final vel.
figure('Name','gain factor peak TOF')
hold on
plot(vel_10,h_10_FM./h_10_NM,'-o')
plot(vel_12,h_12_FM./h_10_NM,'-o')
xlabel('final velocity (m/s)'); ylabel('gain factor peak');legend('10kV_{FM}','12.5kV_{FM}');
hold off




