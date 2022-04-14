%% Compare Phase space of reference and optimization
% Script used to run and compare 

%drop workspace from file in documents/output in here better naming needed

num_p =10e4;


in=InputParameters(12.5, true, 450, 30, 'Phase', 48.83, "FortranSeqBool", false, 'Verbose', false, 'AlwaysGenerateMSeq',true);
in.num_particles = num_p;
in.propagateParticles_euler();
in.output
ref_out = in.output{4,3};
figure()
hold on
in.plot_TOF_laser();

f_ref = Gaussian(ref_out,30)

optimizer_speed_test(in,x_min(1),x_min(2))
in.num_particles = num_p;
in.propagateParticles_euler();
in.output
out_x_min = in.output{4,3};
% in.plot_Phase_dist_atan(x_min(1),x_min(2))
in.plot_TOF_laser();
f_min

optimizer_speed_test(in,bestever.solutions.bestever.x(1),bestever.solutions.bestever.x(2))
in.num_particles = num_p;
in.propagateParticles_euler();
in.output
best_out = in.output{4,3};
% in.plot_Phase_dist_atan(bestever.solutions.bestever.x(1),bestever.solutions.bestever.x(2))
in.plot_TOF_laser();
hold off
f_best = bestever.solutions.bestever.f


ref_out(:,1) = ref_out(:,1) -ref_out(1,1);
ref_out(:,4) = ref_out(:,4) -ref_out(1,4);
label_ref = strings(length(ref_out(:,1)),1)+'Refrence';

out_x_min(:,1) = out_x_min(:,1) - out_x_min(1,1);
out_x_min(:,4) = out_x_min(:,4) - out_x_min(1,4);
label_x_min = strings(length(out_x_min(:,1)),1)+ 'x_{min}';

best_out(:,1) = best_out(:,1) - best_out(1,1);
best_out(:,4) = best_out(:,4) - best_out(1,4);
label_best = strings(length(best_out(:,1)),1) + 'bestever';


comp1 = cat(1,ref_out,out_x_min);
comp2 = cat(1,ref_out,best_out);
label1 = cat(1,label_ref,label_x_min);
label2 = cat(1,label_ref,label_best);

labeling = ["x (m)" , "y (m)", "z (m)"; "v_x (m/s)", "v_y (m/s)", "v_z (m/s)"];
name = [' x vs v_x', 'y vs v_y', 'z vs v_z'];


for i = 1:3
figure()
title(name(i))
scatterhist(comp1(:,i),comp1(:,i+3),'Group',label1,'Color','br','Marker','..',NBins =50)
xlabel(labeling(1,i));ylabel(labeling(2,i));
end


for i = 1:3
figure()
title(name(i))
scatterhist(comp2(:,i),comp2(:,i+3),'Group',label2,'Color','br','Marker','..',NBins =50)
xlabel('y (m)');ylabel('v_y (m/s)');
xlabel(labeling(1,i));ylabel(labeling(2,i));
end

% 
% figure()
% hold on
% 
% xlabel('x (m)');ylabel('v_x (m/s)');
% hold off
% figure()
% hold on
% scatterhist (ref_out(:,2) - ref_out(1,2),ref_out(:,5)-ref_out(1,5))
% xlabel('y (m)');ylabel('v_y (m/s)'); 
% hold off
% figure()
% hold on
% scatterhist (ref_out(:,3) - ref_out(1,3),ref_out(:,6)-ref_out(1,6))
% scatterhist (out_x_min(:,3) - out_x_min(1,3),out_x_min(:,6)-out_x_min(1,6))
% xlabel('z (m)');ylabel('v_z (m/s)'); 
% hold off
% 
% figure()
% hold on
% scatterhist (ref_out(:,1) - ref_out(1,1),ref_out(:,4)-ref_out(1,4))
% scatterhist (best_out(:,1) - best_out(1,1),best_out(:,4)-best_out(1,4))
% xlabel('x (m)');ylabel('v_x (m/s)'); 
% hold off
% figure()
% hold on
% scatterhist (ref_out(:,2) - ref_out(1,2),ref_out(:,5)-ref_out(1,5),'o')
% scatterhist (best_out(:,2) - best_out(1,2),best_out(:,5)-best_out(1,5),'o')
% xlabel('y (m)');ylabel('v_y (m/s)'); 
% hold off
% figure()
% hold on
% scatterhist (ref_out(:,3) - ref_out(1,3),ref_out(:,6)-ref_out(1,6),'o')
% scatterhist (best_out(:,3) - best_out(1,3),best_out(:,6)-best_out(1,6),'o')
% xlabel('z (m)');ylabel('v_z (m/s)');
% hold off

