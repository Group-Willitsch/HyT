%% Redone gain factor

vel_10 = [30,49,90,202,274,328];
vel_12 = [30,49,90,202,274,328];
phase_10_NM = [62.5,62,60.1,50,39.8,30];
phase_10_FM = [61.59,61.05,59.14,49,38.9,29.3];
% phase_12_NM  = [48.83,47.03,39.4,31.5,23.6,15.5];
phase_12_FM = [48.83,48.48,47.12,39.4,31.4,23.6];

[v_10_FM,v_10_NM,] = deal(zeros(size(vel_10)));
[v_12_FM,v_12_NM,] = deal(zeros(size(vel_12)));
n_part = 1e5;

safer_10_FM = cell(length(vel_10),4); % cell array to save all the neseccary data for plotting intialized such that correct length
safer_10_NM = cell(length(vel_10),4); % we safe a flag which mode at which voltage as well as height, area and area of TOF and a cut TOF
safer_12_FM = cell(length(vel_12),4);
save_output = {};

in = InputParameters(10,true, 450, 30, 'Phase',61.55, "FortranSeqBool", false,'Verbose',false,'AlwaysGenerateMSeq',true);
in.num_particles = n_part;
optimizer_atan(in,0,0)
in.propagateParticles_euler();
[h,a,TOF] = in.gain_TOF;
safer_10_FM(1,:) = {'10 kV FM',h,a,TOF};



for i = 2 : length(vel_10)                                  
fprintf('load new velocitiy FM \n')
in.changeVelocities(450, vel_10(i), phase_10_FM(i))  
optimizer_atan(in,0,0)
in.propagateParticles_euler();
[h,a,TOF] = in.gain_TOF;
safer_10_FM(i,:) = {'10 kV FM',h,a,TOF};
h = [];
a = [];
TOF = [];
end

%% FM 12.5 kV
fprintf('loading 12.5 kV FM mode \n')
in.changeFieldConfig(12.5,true) 
for i = 1:length(vel_12)
in.changeVelocities(450, vel_12(i), phase_12_FM(i)) 
optimizer_atan(in,0,0)
in.propagateParticles_euler();
[h,a,TOF] = in.gain_TOF;
safer_12_FM(i,:) = {'12.5 kV FM',h,a,TOF};
h = [];
a = [];
TOF = [];
end

%% NM 10 kV
fprintf('loading 10 kV NM mode \n')
in.changeFieldConfig(10,false)
for i = 1:length(vel_10)
fprintf('load new velocity NM \n')
in.changeVelocities(450, vel_10(i), phase_10_NM(i))
optimizer_atan(in,0,0)
in.propagateParticles_euler();
[h,a,TOF] = in.gain_TOF;
safer_10_NM(i,:) = {'10 kV NM',h,a,TOF};
h = [];
a = [];
TOF = [];
end



a_10_NM = [safer_10_NM{:,3}];       % save all variables locally needed for plotting
a_10_FM = [safer_10_FM{:,3}];       % could also jsut use right side in plot fucntion
a_12_FM = [safer_12_FM{:,3}];

h_10_NM = [safer_10_NM{:,2}];
h_10_FM =[safer_10_FM{:,2}];
h_12_FM=[safer_12_FM{:,2}];



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


