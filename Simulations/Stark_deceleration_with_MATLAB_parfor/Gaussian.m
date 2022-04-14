% Input: 
% x will be the xyzVxyz matrice of the moelcuels that were detected so p.output{4,3} 
% 
% so we first
% need to subtarct the target velocity from v_x for every molecule as well as the pos. of detcetion of the synch molecule
% then we remove the synch molecule and make a multivariante gaussian dist of the pos. and vel.
% Which is not normalized (no prefactor to give a result between 0 and 1
% which then is summed up and 1 is divided by it, to give a result which
% should become smaller if more molecules are closer to the target velocity
% and the pos of the synch molecule
% 

% From PHD thesis Haas page 187 (small book): 
% good estimation of sigma accepted in trap is sigma _vel [~4,2.5,2.5] (v_x,v_y,v_z)
% actually there will be more molecules with vel > target since the beeam
%is hotter than the synch molecule for pos. from same source sigma_pos
%[2,1,1] (pos may need improval)

%


function gauss_vel =Gaussian(x,target_vel)
x(:,4) = x(:,4) - target_vel; % subtarct target velocity (longnitudial)
x(:,1) = x(:,1) - x(1,1); % subtract detection pos of synch molecule
x = x(2:end,:); % cut of synch molecule not of interest
x(:,1:3) = x(:,1:3); %(change units to mm if using normalisation/prefatcor will give small values! to be in same range as vel. (not values far belwo zero))
sigma_vel = [4,2,2];
sigma_pos = [1,0.5,0.5]*1e-3; %change to mm (remove 1e-3 if prefactor applied see above)
sigma = [sigma_pos,sigma_vel];
gauss_vel = zeros(length(x),1);  % intialize vector for results in right dim.
k = length(sigma);    % dimension of mutlivariante gaussian
Cov= eye(length(sigma)).*sigma.^2;
% we take a non correlated covariance matrix thus diagonal

% loop over all row entries correspondig always to one molecule

for i=1:length(x(:,1))
gauss_vel(i) = exp(-0.5*x(i,:)*(Cov^(-1))*transpose(x(i,:)));
end
% set a boundary of the gaussian like here sigma away from center 0
% boundary = 1/(sqrt((2*pi)^k*det(Cov)))*exp(-0.5*sigma*(Cov^(-1))*transpose(sigma));

% then give back the fitness as the 1/number of particles that were above the
% boundary given by sigma in velocity space
gauss_vel = 1/sum(gauss_vel);

% this value is the fitness and should be minimizede by having more
% particles in the accpeted range of the gaussian given by sigma.

end