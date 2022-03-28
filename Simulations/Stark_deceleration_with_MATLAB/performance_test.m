clearvars; close all;

% num particles     time(s)
% normal euler
% 1e4
% 1e5 (2000 at beg.) 15
%                
% faster euler
% 1e4               2.18
% 1e5               5 
% 
% GPU euler     MAYBE ALL WRONG, THE N PARTCILES IS VERY LITTLE
% 1e4               219 (all for last stages...)
% 1e5               218  
% 1e6               219
% 1e7               
    

in = InputParameters(10, false, 450, 30, 'Phase', 62.5, "FortranSeqBool", ...
    false,'Verbose',false, 'AlwaysGenerateMSeq', false);

test = SetOfParticles(in, 1e6);

test.createParticles();

%fprintf('Old Euler runs in ... \t\t')
% tic; test.propagateParticles_euler(); toc;
%fprintf('%i s\n', el_t)


% %% NORMAL MODE 10 kV -- CPU version with local variables and lighter executions
% fprintf('Optimized CPU Euler runs in ... \t\t')
xyzVxyz_0 = test.xyzVxyz_0;
% first propagate to entrance of decelerator
xyzVxyz = xyzVxyz_0;
% 
% ind_particles = transpose(1:size(xyzVxyz,1));
% xyzVxyz(:,1:3) = xyzVxyz(:,1:3) + xyzVxyz(:,4:6) * in.params.FLY_incoupling_time;
% xyzVxyz = xyzVxyz(~(xyzVxyz(:,1) < 0.6765 & (abs(xyzVxyz(:,2)) > 0.001 | abs(xyzVxyz(:,3)) > 0.001)), :); % remove hit particles
% xyzVxyz = xyzVxyz((abs(xyzVxyz(:,1) - xyzVxyz(1,1)) < 11e-3), :); % remove particles not in central fish
% 
% % propagate inside dec
dt = 5e-8;
% % local fields
ax_norm_interpl = in.ax_norm_interpl;
ay_norm_interpl = in.ay_norm_interpl;
az_norm_interpl = in.az_norm_interpl;
ax_norm_H_interpl = in.ax_norm_H_interpl;
ay_norm_H_interpl = in.ay_norm_H_interpl;
az_norm_H_interpl = in.az_norm_H_interpl;
% 
M_time_vec = in.M_time_vec;
% % clearvars in
% dxyzVxyz = {@(y) [y(:, 4:6), ax_norm_interpl(y(:, 1), y(:, 2), y(:, 3)),...
%                                  ay_norm_interpl(y(:, 1), y(:, 2), y(:, 3)),...
%                                  az_norm_interpl(y(:, 1), y(:, 2), y(:, 3)) ]                   
%             @(y) [y(:, 4:6), ax_norm_H_interpl(y(:, 1), y(:, 3), -y(:, 2)),...
%                                  -az_norm_H_interpl(y(:, 1), y(:, 3), -y(:, 2)), ...
%                                  ay_norm_H_interpl(y(:, 1), y(:, 3), -y(:, 2)) ] };
% x = xyzVxyz(:, 1:3);
% V = xyzVxyz(:, 4:6);
% % x = x';
% % V = V';
% 
% tic;
% 
% for i=M_time_vec(1):dt:M_time_vec(end-2)
%     V = V( (abs(x(:,2)) < 0.001 & abs(x(:,3)) < 0.001 & abs(x(:, 1) - x(1, 1)) < 5.5e-3 ), :); % watch out not to put the x before!
%     x = x( (abs(x(:,2)) < 0.001 & abs(x(:,3)) < 0.001 & abs(x(:, 1) - x(1, 1)) < 5.5e-3 ), :);
% 
%     x = x + V*dt;
%         if mod(i, 2) == 0
%         V = V + [ ax_norm_H_interpl(x(:, 1), x(:, 3), -x(:, 2)),...
%                  -az_norm_H_interpl(x(:, 1), x(:, 3), -x(:, 2)), ...
%                  ay_norm_H_interpl(x(:, 1), x(:, 3), -x(:, 2)) ] * dt ;
%         else
%             V = V + [ax_norm_interpl(x(:, 1), x(:, 2), x(:, 3)),...
%                  ay_norm_interpl(x(:, 1), x(:, 2), x(:, 3)),...
%                  az_norm_interpl(x(:, 1), x(:, 2), x(:, 3)) ] * dt ;
%         end
% end
% 
% 
% x = xyzVxyz(:, 1);
% y = xyzVxyz(:, 2);
% z = xyzVxyz(:, 3);
% Vx = xyzVxyz(:, 4);
% Vy = xyzVxyz(:, 5);
% Vz = xyzVxyz(:, 6);
% fprintf('Start with # particles: \t\t\t' + string(length(x)) + '\n')
% 
% for i = 1: 1: (length(M_time_vec) - 2)
% %     xyzVxyz = xyzVxyz( (abs(xyzVxyz(:,2)) < 0.001 & abs(xyzVxyz(:,3)) < 0.001 & abs(xyzVxyz(:, 1) - xyzVxyz(1, 1)) < 5.5e-3 ), :);
% %     V = V( (abs(x(:,2)) < 0.001 & abs(x(:,3)) < 0.001 & abs(x(:, 1) - x(1, 1)) < 5.5e-3 ), :); % watch out not to put the x before!
% %     x = x( (abs(x(:,2)) < 0.001 & abs(x(:,3)) < 0.001 & abs(x(:, 1) - x(1, 1)) < 5.5e-3 ), :);
%     hit = (abs(y) < 0.001 & abs(z) < 0.001 & abs(x - x(1) ) < 5.5e-3 );
%     Vx = Vx( hit ); % watch out not to put the x before!
%     Vy = Vy( hit ); % watch out not to put the x before!
%     Vz = Vz( hit ); % watch out not to put the x before!
%     x = x(hit);
%     y = y(hit);
%     z = z(hit);
% %     length(x                      )
% 
% %     if mod(i, 2) == 0
%     for t = M_time_vec(i):dt: M_time_vec(i+1)
%             x = x + Vx*dt;
%             y = y + Vy*dt;
%             z = z + Vz*dt;
%             if mod(i, 2) == 0
%                 Vx = Vx + dt * ax_norm_H_interpl(x, z, -y);
%                 Vy = Vy + dt * -az_norm_H_interpl(x, z, -y);
%                 Vz = Vz + dt * ay_norm_H_interpl(x, z, -y);
% %             V = V + [ ax_norm_H_interpl(x(:, 1), x(:, 3), -x(:, 2)),...
% %                      -az_norm_H_interpl(x(:, 1), x(:, 3), -x(:, 2)), ...
% %                      ay_norm_H_interpl(x(:, 1), x(:, 3), -x(:, 2)) ] * dt ;
%             else
%                 Vx = Vx + dt * ax_norm_interpl(x, y, z);
%                 Vy = Vy + dt * -ay_norm_interpl(x, y, z);
%                 Vz = Vz + dt * az_norm_interpl(x, y, z);
% %                 V = V + [ax_norm_interpl(x(:, 1), x(:, 2), x(:, 3)),...
% %                      ay_norm_interpl(x(:, 1), x(:, 2), x(:, 3)),...
% %                      az_norm_interpl(x(:, 1), x(:, 2), x(:, 3)) ] * dt ;
%             end
%    end
%             
% %             xyzVxyz = xyzVxyz + [xyzVxyz(:, 4:6), ax_norm_H_interpl(xyzVxyz(:, 1), xyzVxyz(:, 3), -xyzVxyz(:, 2)),...
% %                                  -az_norm_H_interpl(xyzVxyz(:, 1), xyzVxyz(:, 3), -xyzVxyz(:, 2)), ...
% %                                  ay_norm_H_interpl(xyzVxyz(:, 1), xyzVxyz(:, 3), -xyzVxyz(:, 2)) ] .* dt ;
% % 
% %             xyzVxyz = xyzVxyz + [xyzVxyz(:, 4:6), ax_norm_interpl(xyzVxyz(:, 1), xyzVxyz(:, 2), xyzVxyz(:, 3)),...
% %                                  ay_norm_interpl(xyzVxyz(:, 1), xyzVxyz(:, 2), xyzVxyz(:, 3)),...
% %                                  az_norm_interpl(xyzVxyz(:, 1), xyzVxyz(:, 2), xyzVxyz(:, 3)) ] .* dt;
% %             x = x + V*dt;
% %             xyzVxyz = xyzVxyz + [xyzVxyz(:, 4:6), ax_norm_interpl(xyzVxyz(:, 1), xyzVxyz(:, 2), xyzVxyz(:, 3)),...
% %                                  ay_norm_interpl(xyzVxyz(:, 1), xyzVxyz(:, 2), xyzVxyz(:, 3)),...
% %                                  az_norm_interpl(xyzVxyz(:, 1), xyzVxyz(:, 2), xyzVxyz(:, 3)) ] .* dt;
% end
% toc;
% % fprintf('%i s\n', el_t)
% fprintf('Particles survived: \t\t\t' + string(length(x)) + '\n')



%% TRYING A GPU ACCELERATED VERSION
% mayeb usefull: https://ch.mathworks.com/matlabcentral/fileexchange/47437-3d-linear-interpolation-for-gpu
% module load cuda/8.0.44 # add this to the bash, but check the cuda
% version!
gpuDevice(2); % select GPU muber 2, so the screen can run normally...

% (on GPU single precision should be much faster)

gM_time_vec = gpuArray(M_time_vec);

% x, y, z coordinates of big matrix
num_grids_x = 111 + 110 * 122 + 20; % AT LEAST 123, MAYBE 124??
num_grids_y = 41 + 20;
num_grids_z = 41 + 20;
gridded_x = gpuArray( linspace(-10/in.params.SIMION_grid_units_p_meter, in.params.PHYS_length_dec + 10/in.params.SIMION_grid_units_p_meter, num_grids_x)' );
gridded_y = gpuArray( linspace(-10/in.params.SIMION_grid_units_p_meter - in.params.PHYS_seperation_pins/2.0, in.params.PHYS_seperation_pins/2.0 + 10/in.params.SIMION_grid_units_p_meter, num_grids_y)');
gridded_z =	gpuArray( linspace(-10/in.params.SIMION_grid_units_p_meter - in.params.PHYS_seperation_pins/2.0, in.params.PHYS_seperation_pins/2.0 + 10/in.params.SIMION_grid_units_p_meter, num_grids_z)');
% big matrices
gax_norm_extended = gpuArray( in.ax_norm_extended );
gay_norm_extended = gpuArray( in.ay_norm_extended );
gaz_norm_extended = gpuArray( in.az_norm_extended );
gax_norm_H_extended = gpuArray( in.ax_norm_H_extended );
gay_norm_H_extended = gpuArray( in.ay_norm_H_extended );
gaz_norm_H_extended = gpuArray( in.az_norm_H_extended );


gxyzVxyz_0 = gpuArray(xyzVxyz_0); 
gxyzVxyz = gpuArray(xyzVxyz_0);
gxyzVxyz(:,1:3) = gxyzVxyz(:,1:3) + gxyzVxyz(:,4:6) * in.params.FLY_incoupling_time;
gxyzVxyz = gxyzVxyz( (abs(gxyzVxyz(:,2)) < 0.001 & abs(gxyzVxyz(:,3)) < 0.001 & abs(gxyzVxyz(:, 1) - gxyzVxyz(1, 1)) < 5.5e-3 ), :);

times = gpuArray( zeros(length(gM_time_vec)) );
n_par = gpuArray( zeros(length(gM_time_vec)) );
% interpn(X1,X2,...,Xn,V,Xq1,Xq2,...,Xqn) 
size(gxyzVxyz)
for i = 1: 1: (length(gM_time_vec) - 2)
    tic;
    size(gxyzVxyz)
%     gxyzVxyz = gxyzVxyz( (abs(gxyzVxyz(:,2)) < 0.001 & ( abs(gxyzVxyz(:,3)) < 0.001)  & abs(gxyzVxyz(:, 1) - gxyzVxyz(1, 1)) < 5.5e-3 ), :); % remove hit particles
    if mod(i, 2) == 0
        for t = gM_time_vec(i) : dt : gM_time_vec(i+1)
            
            gxyzVxyz = gxyzVxyz + [gxyzVxyz(:, 4:6), ...
                                interpn(gridded_x, gridded_y, gridded_z, gax_norm_H_extended, gxyzVxyz(:, 1), gxyzVxyz(:, 3), gxyzVxyz(:, 2), 'linear' ), ...
                                interpn(gridded_x, gridded_y, gridded_z, gaz_norm_H_extended, gxyzVxyz(:, 1), gxyzVxyz(:, 3), gxyzVxyz(:, 2), 'linear' ), ...
                                interpn(gridded_x, gridded_y, gridded_z, gay_norm_H_extended, gxyzVxyz(:, 1), gxyzVxyz(:, 3), gxyzVxyz(:, 2), 'linear' )] .*dt ;
%                                  ax_norm_H_interpl(gxyzVxyz(:, 1), gxyzVxyz(:, 3), -gxyzVxyz(:, 2)),...
%                                  -az_norm_H_interpl(gxyzVxyz(:, 1), gxyzVxyz(:, 3), -gxyzVxyz(:, 2)), ...
%                                  ay_norm_H_interpl(gxyzVxyz(:, 1), gxyzVxyz(:, 3), -gxyzVxyz(:, 2)) ] .* dt ;
        end
    else
        for t = gM_time_vec(i) : dt : gM_time_vec(i+1)
            gxyzVxyz = gxyzVxyz + [gxyzVxyz(:, 4:6), ...
                interpn(gridded_x, gridded_y, gridded_z, gax_norm_extended, gxyzVxyz(:, 1), gxyzVxyz(:, 2), gxyzVxyz(:, 3), 'linear' ), ...
                interpn(gridded_x, gridded_y, gridded_z, gay_norm_extended, gxyzVxyz(:, 1), gxyzVxyz(:, 2), gxyzVxyz(:, 3), 'linear' ), ...
                interpn(gridded_x, gridded_y, gridded_z, gaz_norm_extended, gxyzVxyz(:, 1), gxyzVxyz(:, 2), gxyzVxyz(:, 3), 'linear' )] .*dt;
%                 ax_norm_interpl(xyzVxyz(:, 1), xyzVxyz(:, 2), xyzVxyz(:, 3)),...
%                 ay_norm_interpl(xyzVxyz(:, 1), xyzVxyz(:, 2), xyzVxyz(:, 3)),...
%                 az_norm_interpl(xyzVxyz(:, 1), xyzVxyz(:, 2), xyzVxyz(:, 3)) ] .* dt;

        end
    end
    % at every decelerator stpe, take time and number of particles for
    % later plotting
    times(i) = toc;
    n_par(i) = length(gxyzVxyz(:, 1));
end
result = gather(gxyzVxyz);
times = gather(times);
n_par = gather(n_par);
figure
plot(times, n_par); xlabel('Comp. time per stage (s)'); ylabel('Number of particle at the end of that deceleration stage')
