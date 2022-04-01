clear; close all;

tic;
parfor i=1:6
    in=InputParameters(10, false, 450, 30, 'Phase', 62.5, "FortranSeqBool", false, 'Verbose', false);
    in.num_particles = 100000;
    in.propagateParticles_euler();
end
toc;

