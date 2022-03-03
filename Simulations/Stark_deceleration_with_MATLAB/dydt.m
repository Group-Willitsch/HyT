function dydt = dydt(t, q, my_params)
% ODE that mimics the decelerator
% it is encoded as dq/dt = f (q)
% where q  = (x, y, z, Vx, Vy, Vz)
% and dq/dt = ( Vx, Vy, Vz, -m*omega^2*x, -m*omega^2*y, 0)
% the potential energy is U = 0.5m*omega^2*(x^2+y^2+z^2)
% the kinetic energy is K = 0.5 * m * (Vx^2 + Vy^2 + Vz^2)
% since these are simply vector equation, no need for a "Mass matrix".
% Otherwise we would have to add it.

% parameters is a struct that contains the parameters m and omega
% extract as parameters.m and parameters.omega

% dydt = [q(4), q(5), q(6), -my_params.kappa * 1.05 * q(1), -my_params.kappa * q(2) + 1, 0]';
dydt = [q(4); q(5); q(6); -my_params.kappax * q(1); -my_params.kappa * q(2) ; 0];
% dydt = my_params.matrix * q ; % seems not to change much wrt to normal
% one

end

