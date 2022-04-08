% Input: 
% x will be the xyzVxyz matrice of the moelcuels that were detected so p.output{4,3} so we first
% need to cut off the positions as well as the synch molecule and
% then subtarct the target velocity from v_x for every molecule.
%target vel is zero in x, y direction
% 

% From PHD thesis Haas page 187 (small book): 
% good estimation of sigma accepted in trap is [~8,5,5] (v_x,v_y,v_z)
% actually there will be more molecules with vel > target since the beeam
%is hotter than the synch molecule

% sigma is the standart derivation used to define the gaussian and defines
% which molecules we want to accept as close enough to the target vel

function particles_close_2target_vel =Gaussian(x,sigma,target_vel)
x = x(2:end,4:end); % cut of positions and synch molecule
x(:,1) = x(:,1) - target_vel;
gauss_vel = zeros(length(x),1);  % intialize vector for results in right dim.
k = length(sigma);    % dimension of mutlivariante gaussian
Cov=eye(length(sigma)).*sigma.^2; % we take a non correlated covariance matrix thus diagonal

% loop over all row entries correspondig always to one molecule

for i=1:length(x(:,1))
gauss_vel(i) = exp(-0.5*x(i,:)*(Cov^(-1))*transpose(x(i,:)));
end
% set a boundary of the gaussian like here sigma away from center 0
% boundary = 1/(sqrt((2*pi)^k*det(Cov)))*exp(-0.5*sigma*(Cov^(-1))*transpose(sigma));

% then give back the fitness as the 1/number of particles that were above the
% boundary given by sigma in velocity space
particles_close_2target_vel = gauss_vel;

% this value is the fitness and should be minimizede by having more
% particles in the accpeted range of the gaussian given by sigma.

end