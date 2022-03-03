#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 09:38:39 2019

@author: claudiovonplanta
"""
###https://github.com/choderalab/openmm-tutorials
from simtk import openmm, unit
from simtk.openmm import app
import numpy as np
import math
import matplotlib
#matplotlib.use('Agg')

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

import scipy.constants as con

from numpy.linalg import norm
import numpy as np
from cazeeman_import import Bmodel
from cazeeman_import import Bvecs

import time
#import h5py
import subprocess

from scipy.integrate import ode
from scipy.linalg import lu_factor, lu_solve

muB = 9.27400994e-24    # J/T
amtoJ = 4.359744e-18    # conversion au to J
hbar = con.hbar
to_mm = 1e3
kb = con.k

def bug():
    plt.show()
    import sys
    sys.exit()


def velT(m, T):
    return np.sqrt(3 * con.k * T /m)

def sigv(T, m):
    return np.sqrt(kb * T/m)

def angle_3d(v1, v2):
    return np.arccos(np.dot(v1,v2) / (norm(v1) * norm(v2)))

def seady_state_OBE(u, alpha, beta, D1, D2, s1, s2, l1, l2):
    Ln = L_Matrix_oberst(u, alpha, beta, D1, D2, s1, s2, l1, l2)
    pops = steadystate_populations(Ln)
    return pops

def u_from_B(B):
    return B * muB/hbar / (2*np.pi)

def doppler(kz, vz):
	return -np.dot(kz, vz) / (2*np.pi)
    #~ return -kz*vz / (2*np.pi)

def steadystates(us, alphas, betas, d1s, d2s, s1, s2, l1, l2):
    pops = np.zeros((us.shape[0],8))

    for i in range(us.shape[0]):
        popslist = seady_state_OBE(us[i], alphas[i], betas[i], d1s[i], d2s[i], s1, s2, l1, l2)
        pops[i,:] = popslist
    return pops

def L_Matrix_oberst(u, alpha, beta, D1, D2, s1, s2, l1, l2):
    id8 = np.identity(8)
    v1, v2, v3, v4, v5, v6, v7, v8 = id8
    v1, v2, v3, v4, v5, v6, v7, v8 = v1.reshape(8,1), v2.reshape(8,1), v3.reshape(8,1), v4.reshape(8,1), v5.reshape(8,1), v6.reshape(8,1), v7.reshape(8,1), v8.reshape(8,1)

    g1 = 21.58
    g2 = 1.35

    O1 = s1*g1
    O2 = s2*g2

    gS = 2
    gP = 2/3
    gD = 4/5
    ca = np.cos(alpha)
    sa = np.sin(alpha)
    cb = np.cos(beta)
    sb = np.sin(beta)

    w3 = np.sqrt(3)

    H = np.diag(np.array([D1, D1, 0, 0, D2, D2, D2, D2]) + 0.5 * u * np.array([-gS, gS, -gP, gP, -3*gD, -gD, gD, 3*gD]))
    H = np.matrix(H)
    HSP = -O1/w3 * np.array([[-ca, sa],[sa, ca]])
    HDP = -O2 / 2 / w3 * np.array([[w3 * sb, 0],[2*cb, sb],[-sb, 2*cb],[0, -w3*sb]])
    HSP = np.matrix(HSP)
    HDP = np.matrix(HDP)
    H[0:2, 2:4] = HSP
    H[2:4, 0:2] = HSP.H
    H[4:8, 2:4] = HDP
    H[2:4, 4:8] = HDP.H

    C1 = np.matrix(math.sqrt(2/3 * g1) * np.kron(v1, v4.T))
    C2 = np.matrix(math.sqrt(2/3 * g1) * np.kron(v2, v3.T))
    C3 = np.matrix(math.sqrt(1/3 * g1) * (np.kron(v1, v3.T) - np.kron(v2, v4.T)))
    C4 = np.matrix(math.sqrt(1/2 * g2) * np.kron(v5, v3.T) + math.sqrt(1/6 * g2) * np.kron(v6, v4.T))
    C5 = np.matrix(math.sqrt(1/6 * g2) * np.kron(v7, v3.T) + math.sqrt(1/2 * g2) * np.kron(v8, v4.T))
    C6 = np.matrix(math.sqrt(1/3 * g2) * (np.kron(v6, v3.T) + np.kron(v7, v4.T)))
    C7 = np.matrix(math.sqrt(2 * l1) * (np.kron(v1, v1.T) + np.kron(v2, v2.T)))
    C8 = np.matrix(math.sqrt(2 * l1) * (np.kron(v5, v5.T) + np.kron(v6, v6.T) + np.kron(v7, v7.T) + np.kron(v8, v8.T)))
    CC = np.dot(C1.H, C1) + np.dot(C2.H, C2) + np.dot(C3.H, C3) + np.dot(C4.H, C4) + np.dot(C5.H, C5) + np.dot(C6.H, C6) + np.dot(C7.H, C7) + np.dot(C8.H, C8)
    H = H - 1j/2 * CC

    Ckron = np.kron(C1, C1) + np.kron(C2, C2) + np.kron(C3, C3) + np.kron(C4, C4) + np.kron(C5, C5) + np.kron(C6, C6) + np.kron(C7, C7) + np.kron(C8, C8)

    L = -1j * (np.kron(H, id8) - np.kron(id8, H.H)) + Ckron
    return L

def L_Matrix_oberst_sinco(u, ca, sa, cb, sb, D1, D2, s1, s2, l1, l2):
    id8 = np.identity(8)
    v1, v2, v3, v4, v5, v6, v7, v8 = id8
    v1, v2, v3, v4, v5, v6, v7, v8 = v1.reshape(8,1), v2.reshape(8,1), v3.reshape(8,1), v4.reshape(8,1), v5.reshape(8,1), v6.reshape(8,1), v7.reshape(8,1), v8.reshape(8,1)

    g1 = 21.58
    g2 = 1.35

    O1 = s1*g1
    O2 = s2*g2

    gS = 2
    gP = 2/3
    gD = 4/5
    #~ ca = np.cos(alpha)
    #~ sa = np.sin(alpha)
    #~ cb = np.cos(beta)
    #~ sb = np.sin(beta)

    w3 = np.sqrt(3)

    H = np.diag(np.array([D1, D1, 0, 0, D2, D2, D2, D2]) + 0.5 * u * np.array([-gS, gS, -gP, gP, -3*gD, -gD, gD, 3*gD]))
    H = np.matrix(H)
    HSP = -O1/w3 * np.array([[-ca, sa],[sa, ca]])
    HDP = -O2 / 2 / w3 * np.array([[w3 * sb, 0],[2*cb, sb],[-sb, 2*cb],[0, -w3*sb]])
    HSP = np.matrix(HSP)
    HDP = np.matrix(HDP)
    H[0:2, 2:4] = HSP
    H[2:4, 0:2] = HSP.H
    H[4:8, 2:4] = HDP
    H[2:4, 4:8] = HDP.H

    C1 = np.matrix(math.sqrt(2/3 * g1) * np.kron(v1, v4.T))
    C2 = np.matrix(math.sqrt(2/3 * g1) * np.kron(v2, v3.T))
    C3 = np.matrix(math.sqrt(1/3 * g1) * (np.kron(v1, v3.T) - np.kron(v2, v4.T)))
    C4 = np.matrix(math.sqrt(1/2 * g2) * np.kron(v5, v3.T) + math.sqrt(1/6 * g2) * np.kron(v6, v4.T))
    C5 = np.matrix(math.sqrt(1/6 * g2) * np.kron(v7, v3.T) + math.sqrt(1/2 * g2) * np.kron(v8, v4.T))
    C6 = np.matrix(math.sqrt(1/3 * g2) * (np.kron(v6, v3.T) + np.kron(v7, v4.T)))
    C7 = np.matrix(math.sqrt(2 * l1) * (np.kron(v1, v1.T) + np.kron(v2, v2.T)))
    C8 = np.matrix(math.sqrt(2 * l1) * (np.kron(v5, v5.T) + np.kron(v6, v6.T) + np.kron(v7, v7.T) + np.kron(v8, v8.T)))
    CC = np.dot(C1.H, C1) + np.dot(C2.H, C2) + np.dot(C3.H, C3) + np.dot(C4.H, C4) + np.dot(C5.H, C5) + np.dot(C6.H, C6) + np.dot(C7.H, C7) + np.dot(C8.H, C8)
    H = H - 1j/2 * CC

    Ckron = np.kron(C1, C1) + np.kron(C2, C2) + np.kron(C3, C3) + np.kron(C4, C4) + np.kron(C5, C5) + np.kron(C6, C6) + np.kron(C7, C7) + np.kron(C8, C8)

    L = -1j * (np.kron(H, id8) - np.kron(id8, H.H)) + Ckron
    return L

def steadystate_populations(M):
    # exchanging one of the equations (row in matrix) by the normalization condition
    num = 1 # this is the second last equation (we subtract 1 later)
    normvec = np.zeros(64)
    for i in range(8):
        normvec[i*8 + i] = 1
    M[num-1,:] = normvec

    # steady state solutions are 0 everywhere but for the normalization 1
    b = np.zeros(64)
    b[num-1] = 1

    #~ smallvalue = 1e-6
    #~ M = M + smallvalue

#    LU decomposition
    lu, piv = lu_factor(M)
    sol = lu_solve((lu, piv), b, trans = 0, check_finite=True)

    # populations (every 8th entry)
    pops = np.array([sol[i*8 + i].real for i in range(8)])
    return pops

def timestep(M, rho, dt):
    drho = np.dot(M, rho)
    rho = rho - drho * dt
    return rho



class PaulTrap():
    def __init__(self):
        self.localPi = np.pi

        self.rffreq = 8e6
        self.wavelength_resonance = 396.950150e-9
        self.freq_resonance = con.c / self.wavelength_resonance
        self.detuning = - 100e6             ### In Hz
        self.freq_abs = self.freq_resonance + self.detuning
        self.wavelength_abs = con.c / self.freq_abs
        self.kz = 2 * self.localPi / self.wavelength_abs
        self.kzvecn = np.array([1, 0, 1]) * 1/np.sqrt(2)
        self.kzvec = self.kz * self.kzvecn


        self.g1 = 21.58 * 1e6 * 2*np.pi
        self.g2 = 1.35 * 1e6 * 2*np.pi
        self.linewidths = np.array([self.g1,self.g1,self.g1,self.g1,self.g2,self.g2,self.g2,self.g2])
        self.ecdista = 12.5e-3
        self.ecdistb = 12.5e-3
        self.VEC1 = 5
        self.VEC2 = 5
        self.VRF = 150


        self.kappa = 0.89
        self.z0 = 2.75e-3
        self.r0 = 3.5e-3

        self.qx = 0.0824
        self.ax = -16e-4
        self.qy = -3e-4
        self.ay = 4e-4
        self.qz = -0.0806
        self.az = 11e-4

        # initial position parameters
        self.xAxisLength = 20e-6
        self.yAxisLength = 20e-6
        self.zAxisLength = 3000e-6

        self.nCa = 1

        self.camassSI = 40 * 1.6e-27
        self.chargeSI = 1.0 * 1.6e-19

        self.charge=  1.0 * unit.elementary_charge
        self.camass = 40 * unit.amu

        print('mass', self.camass)
        print('charge', self.charge)

        #
        self.timestep = 10.0  * unit.nanoseconds
        self.Temp = 1e-3
        self.dt   = 5

        self.currentTime = 0.0

        self.platf = 'CPU'
        #self.platf = 'OpenCL'
        self.sigma = 0
        self.epsilon = 0

        self.numphotons = 0

        self.freezeconst = 0.5

        self.mode = 'read'
        self.mode = 'init'

        #~ laser and polarizations n stuff

        self.d1 = -20
        self.d2 = 0
        self.l1 = 0.3
        self.l2 = 0.3
        self.s1 = 0.54
        self.s2 = 3.9
        wavelength_397 = 397e-9
        self.k397 = 2*np.pi / wavelength_397
        wavelength_866 = 866e-9
        self.k866 = 2*np.pi / wavelength_866
        self.kz2 = 2 * self.localPi / wavelength_866
        self.kz2vecn = np.array([1, 0, 1]) * 1/np.sqrt(2)
        self.kzvec2 = self.kz2 * self.kz2vecn
        self.CaState = 0



    def initialize(self):
        self.ionLevels = np.zeros(self.nCa)
        #initialize
        self.system = openmm.System()

        # create some particles
        for index in range(self.nCa):
            # Particles are added one at a time
            # Their indices in the System will correspond with their indices in the Force objects we will add later
            self.system.addParticle(self.camass)


        # custom Paul trap force, fitted \sum_{i=0,j=0}^2 a_{i,j} * x^i * y^j to comsol data
        to_nm = 1e-9
        Forcestring = "charge*paulTrapPrefactor*("
        Forcestring += "VEC1 * (ddc * z*z + edc * z*z*z*z) * adc * exp(-bdc * x*x - cdc * y*y)"
        Forcestring += "+vrf*rffactor*"
        Forcestring += "(0.0010654868221601826 + 0.0018288104115380941*x*to_nm + 229.48604865396135*x*x*to_nm*to_nm + 0.0030131554295328485*y*to_nm + 849179.5330164384*x*to_nm*y*to_nm - 2394.692253439868*x*x*y*to_nm*to_nm*to_nm +"
        Forcestring += "229.99761958659172*y*y*to_nm*to_nm - 3043.444228836036*x*y*y*to_nm*to_nm*to_nm - 1.6378374076833544e9*x*x*y*y*to_nm*to_nm*to_nm*to_nm))"
        self.paulTrapForce = openmm.CustomExternalForce(Forcestring)
        self.paulTrapForce.addGlobalParameter("rffactor", math.cos(self.rffreq * self.currentTime * 2 * self.localPi))
        self.paulTrapForce.addGlobalParameter("ecdista", self.ecdista)
        self.paulTrapForce.addGlobalParameter("ecdistb", self.ecdistb)
        self.paulTrapForce.addGlobalParameter("paulTrapPrefactor", 1.0/(2.0*self.localPi) *96.32 )
        self.paulTrapForce.addGlobalParameter("VEC1", self.VEC1)
        self.paulTrapForce.addGlobalParameter("vrf", self.VRF)
        self.paulTrapForce.addGlobalParameter("charge", self.charge)
        self.paulTrapForce.addGlobalParameter("to_nm", to_nm)
        self.paulTrapForce.addGlobalParameter("adc", 2.73876265e+02)
        self.paulTrapForce.addGlobalParameter("bdc", 4.71464738e+04 *to_nm*to_nm)
        self.paulTrapForce.addGlobalParameter("cdc", 1.67635287e+05 *to_nm*to_nm)
        self.paulTrapForce.addGlobalParameter("ddc", 3.66919525e+02 *to_nm*to_nm)
        self.paulTrapForce.addGlobalParameter("edc", 2.38492048e+07 *to_nm*to_nm*to_nm*to_nm)


#
        for index in range(self.nCa):
            self.paulTrapForce.addParticle(index)
        self.system.addForce(self.paulTrapForce)

        # define particle particle force
        self.nonbond = openmm.NonbondedForce()

        for index in range(self.nCa): # all particles must have parameters assigned for the NonbondedForce
            # Particles are assigned properties in the same order as they appear in the System object
            self.nonbond.addParticle(self.charge, self.sigma, self.epsilon)

        self.system.addForce(self.nonbond)

#        some more options for p-p interactions
#        nonbond.setCutoffDistance(3.0 * sigma) # set cutoff (truncation) distance at 3*sigma
        #nonbond.setUseSwitchingFunction(True) # use a smooth switching function to avoid force discontinuities at cutoff
        #nonbond.setSwitchingDistance(2.5 * sigma) # turn on switch at 2.5*sigma
        #nonbond.setUseDispersionCorrection(True) # use long-range isotropic dispersion correction
        #nonbond_index = system.addForce(nonbond) # system takes ownership of the NonbondedForce object


        # positions
        thetaRand = np.arccos(2 * np.random.rand(self.nCa) - 1)
        phiRand = 2 * self.localPi * np.random.rand(self.nCa)
        rRand = 1e9 * np.random.rand(self.nCa) **(1/3)

        x0 = self.xAxisLength * rRand * np.sin(thetaRand) * np.cos(phiRand)
        y0 = self.yAxisLength * rRand * np.sin(thetaRand) * np.sin(phiRand)
        z0 = self.zAxisLength * rRand * np.cos(phiRand)

        #x0 = np.array([1])
        #y0 = np.array([1])
        #z0 = np.array([1])


        if self.mode == 'init':
            self.positions = np.array([x0, y0, z0]).reshape([self.nCa, 3])


        #         printf("defining integrator\n");
        #
        #        omm->integrator = new OpenMM::CustomIntegrator(timeStep);
        #        //adds a global velocity kick size variable for use in the integration
        #        kickTermIndex = omm->integrator->addGlobalVariable("velKickSize", velocityKickSize / 1000.0);
        #        omm->integrator->addPerDofVariable("fs", 0); //adds in a dummy variable to store the force at the start of the timestep
        #        omm->integrator->addComputePerDof("fs", "f"); //sets the dummy variable to the current force
        #        omm->integrator->addComputePerDof("x", "x+dt*v+0.5*dt*dt*f/m"); //updates the positions
        #        omm->integrator->addComputePerDof("v", "v+0.5*dt*(f+fs)/m"); //update the velocity - f is now the force with new positions, so we use the saved force from start of timestep too.
        #        //this also includes a random velocity kick
        #        omm->integrator->addUpdateContextState();


        # Create an integrator
        self.integrator = openmm.VerletIntegrator(self.timestep)
        #~ self.integrator = openmm.VariableVerletIntegrator(self.timestep)

        # Create a Context using the default platform (the fastest abailable one is picked automatically)
        # NOTE: The integrator is bound irrevocably to the context, so we can't reuse it
        self.platform = openmm.Platform.getPlatformByName(self.platf)
        self.context = openmm.Context(self.system, self.integrator, self.platform)
        # Set the positions
        self.context.setPositions(self.positions)

        if self.mode == 'init':
            self.context.setVelocitiesToTemperature(self.Temp) # set the velocities to Boltzmann distribution
        if self.mode == 'read':
            self.context.setVelocities(self.velocities)
        self.context.setParameter("rffactor", 1)

        # Retrieve the energy and forces
        self.state = self.context.getState(getPositions = True, getVelocities = True, getForces = True)
        self.populations = np.zeros([self.nCa,1])

        #~ (u, ca, sa, cb, sb, D1, D2, s1, s2, l1, l2)
        self.L0 = L_Matrix_oberst_sinco(0, 0, 0, 0, 0, 0, 0, self.s1, self.s2, self.l1, self.l2)
        self.L0_u = L_Matrix_oberst_sinco(1, 0, 0, 0, 0, 0, 0, self.s1, self.s2, self.l1, self.l2) - self.L0
        self.L0_ca = L_Matrix_oberst_sinco(0, 1, 0, 0, 0, 0, 0, self.s1, self.s2, self.l1, self.l2) - self.L0
        self.L0_sa = L_Matrix_oberst_sinco(0, 0, 1, 0, 0, 0, 0, self.s1, self.s2, self.l1, self.l2) - self.L0
        self.L0_cb = L_Matrix_oberst_sinco(0, 0, 0, 1, 0, 0, 0, self.s1, self.s2, self.l1, self.l2) - self.L0
        self.L0_sb = L_Matrix_oberst_sinco(0, 0, 0, 0, 1, 0, 0, self.s1, self.s2, self.l1, self.l2) - self.L0
        self.L0_d1 = L_Matrix_oberst_sinco(0, 0, 0, 0, 0, 1, 0, self.s1, self.s2, self.l1, self.l2) - self.L0
        self.L0_d2 = L_Matrix_oberst_sinco(0, 0, 0, 0, 0, 0, 1, self.s1, self.s2, self.l1, self.l2) - self.L0

    def obedynamics(self,state, pops):
        kicks1 = np.array([0,0,0])
        kicks2 = np.array([0,0,0])
        randomkicks1 = 0
        randomkicks2 = 0
        
        dt = 5e-9


        if state==0 or state==1:
            rand = np.random.random()
            rates = pops * self.g1 * dt
            if np.sum(rates)>1:
                print('timestep too large 0,1')

            if rand<rates[2]:
                kicks1[2] += 1
                state = 2
            elif rand<(rates[2]+rates[3]):
                kicks1[2] += 1
                state = 3

        elif state==2 or state==3:
            decayprobs = (1) * self.linewidths * dt # only spontaneous emission so far
            rand = np.random.random()
            if rand<decayprobs[0]:
                randomkicks1 += 1
                state = 0
            elif rand<(decayprobs[0]+decayprobs[1]):
                randomkicks1 += 1
                state = 1
            elif rand<(decayprobs[0]+decayprobs[1]+decayprobs[4]):
                randomkicks2 += 1
                state = 4
            elif rand<(decayprobs[0]+decayprobs[1]+decayprobs[4]+decayprobs[5]):
                randomkicks2 += 1
                state = 5
            elif rand<(decayprobs[0]+decayprobs[1]+decayprobs[4]+decayprobs[5]+decayprobs[6]):
                randomkicks2 += 1
                state = 6
            elif rand<(decayprobs[0]+decayprobs[1]+decayprobs[4]+decayprobs[5]++decayprobs[6]++decayprobs[7]):
                randomkicks2 += 1
                state = 7

        elif state>3:
            rand = np.random.random()
            rates = pops * self.g2 * dt
            if np.sum(rates)>1:
                print('timestep too large >3')
            if rand<rates[2]:
                kicks2[2] += 1
                state = 2
            elif rand<rates[2]+rates[3]:
                kicks2[2] += 1
                state = 3

        return state, kicks1, kicks2, randomkicks1, randomkicks2

    def doStep(self, steps):
        self.state = self.context.getState(getPositions = True, getVelocities = True, getEnergy = True)
        self.positions = self.state.getPositions(asNumpy = True)
        self.ekin = self.state.getKineticEnergy()._value / con.N_A * 1e3
        self.epot = self.state.getPotentialEnergy()._value / con.N_A * 1e3
        self.etot = self.ekin + self.epot

        self.velocities = self.state.getVelocities(asNumpy = True)
        self.currentTime = self.state.getTime()._value * 1e-12
        self.context.setParameter("rffactor", math.cos(self.rffreq * self.currentTime * 2 * self.localPi))
        vz = self.velocities._value[:,2] * 1e3 # nm/ps to m/s

        xs = self.positions[:,0]._value * 1e-9 # nm to m
        ys = self.positions[:,1]._value * 1e-9
        zs = self.positions[:,2]._value * 1e-9
        rsq = xs**2 + ys**2

        Bvectors = Bvecs(xs,ys,zs).T
        Bmagnitudes = Bmodel(xs,ys,zs)
        Bmagnitudes = np.ones_like(Bmagnitudes) * 1e-4
        us = u_from_B(Bmagnitudes) * 1e-6
#        print(us)
        alphas = np.ones_like(us) * np.pi/2#
        #~ alphas = angle_3d(Bvectors, self.laser397_pol)
        betas = np.ones_like(us) * np.pi/2
        #~ betas = angle_3d(Bvectors, self.laser866_pol)

        d1s = self.d1 + doppler(self.kzvec, self.velocities._value.T * 1e3 )*1e-6 * 1
        d2s = self.d2 + doppler(self.kzvec2, self.velocities._value.T * 1e3)*1e-6 * 1
        self.d1s = d1s
        self.d2s = d2s

        #~ print(doppler(self.kzvec, self.velocities._value.T ))
        #~ print(d1s)
        #~ print(doppler(self.k397, vz ))
        #~ self.populations = steadystates(us, alphas, betas, d1s, d2s, self.s1, self.s2, self.l1, self.l2)
        ca = np.cos(alphas[0])
        sa = np.sin(alphas[0])
        cb = np.cos(betas[0])
        sb = np.sin(betas[0])
        Ln = self.L0 + us[0] * self.L0_u + ca * self.L0_ca + cb * self.L0_cb + sa * self.L0_sa + sb * self.L0_sb + d1s[0] * self.L0_d1 + d2s[0] * self.L0_d2
        self.populations = np.array([steadystate_populations(Ln)])

        #~ self.populations = np.array([[1/6,1/6,1/6,1/6,1/12,1/12,1/12,1/12]])
        self.CaState, kicks1, kicks2, randomkicks1, randomkicks2 = self.obedynamics(self.CaState, self.populations[0])
        kicks1 = np.linalg.norm(kicks1) * self.kzvecn
        kicks2 = np.linalg.norm(kicks2) * self.kz2vecn

#		# spontaneous emission into S
        for i in range(randomkicks1):
            phi = np.random.rand() * 2 * np.pi
            costheta = 2 * np.random.rand() - 1
            sintheta = np.sqrt(1-costheta**2)
            kicks1[0] += sintheta * np.cos(phi)
            kicks1[1] += sintheta * np.sin(phi)
            kicks1[2] += costheta

#		# spontaneous emission into D
        for i in range(randomkicks2):
            phi = np.random.rand() * 2 * np.pi
            costheta = 2 * np.random.rand() - 1
            sintheta = np.sqrt(1-costheta**2)
            kicks2[0] += sintheta * np.cos(phi)
            kicks2[1] += sintheta * np.sin(phi)
            kicks2[2] += costheta

        vels = self.velocities._value
        self.velocities =  vels + 1 *  con.hbar * (1 * kicks1 * self.kz + 1 * kicks2 * self.k866) / self.camassSI * 1e-3
        self.numphotons += np.linalg.norm(kicks1)
        self.context.setVelocities(self.velocities)
        self.integrator.step(steps)
        return



if __name__ == "__main__":
    
    Trap = PaulTrap()     ###-- Instantation of the class
    
    Trap.mode = 'init'
    if Trap.mode == 'read':
        Trap.positions = np.load('./positions_last.npy')
        Trap.velocities = np.load('./velocities_last.npy')
    
    
    Trap.mode = 'read'
    offset = 100e-6 * 1e9
    Trap.positions = np.array([offset,offset,offset]).reshape((1,3))
    #~ vel = 1/np.sqrt(3) * 558 * 1e-3
    vel = 1/np.sqrt(3) * 100 * 1e-3 * 0
    Trap.velocities = -np.array([vel,vel,vel]).reshape((1,3))
    
    Trap.d1 = -20
    Trap.d2 = -20
    Trap.l1 = 0.3
    Trap.l2 = 0.3
    Trap.s1 = 1
    Trap.s2 = 6
    wavelength_397 = 397e-9
    Trap.k397 = 2*np.pi / wavelength_397
    wavelength_866 = 866e-9
    Trap.k866 = 2*np.pi / wavelength_866
    Trap.Gamma = 2*np.pi*20.7e6
    Trap.laser397_pol = np.array([1/np.sqrt(2),0,1/np.sqrt(2)])
    Trap.laser866_pol = np.array([1/np.sqrt(2),0,1/np.sqrt(2)])
    
    
    dt = 5 # ns
    Trap.currentTime = 0
    Trap.Temp = 100e-3
    Trap.nCa = 1
    Trap.timestep = dt  * unit.nanoseconds
    Trap.rffreq = 8e6
    Trap.VRF = 150
    Trap.VEC1 = 2
    Trap.VEC2 = 0#200 # not used
    Trap.xAxisLength = 200e-6
    Trap.yAxisLength = 200e-6
    Trap.zAxisLength = 1000e-6
    
    
    totaltime = 100e-9
    totalsteps = int(totaltime / (dt*1e-9))
    #totalsteps = int(1e3)
    
    print('Number of Ca atoms', Trap.nCa)
    print('VRF:',Trap.VRF)
    print('RF freq', Trap.rffreq)
    print('Endcaps', Trap.VEC1, Trap.VEC2)
    print('Timestep', Trap.timestep)
    print('Total time', Trap.timestep * totalsteps)
    #Trap.platf = 'OpenCL'
    Trap.platf = 'CPU'
    
    Trap.initialize()
    print('number particles:', Trap.nCa)
    velocities =Trap.state.getVelocities(asNumpy = True)._value
    print('velocities mean', np.mean(np.sqrt(velocities[:,0]**2 + velocities[:,1]**2 + velocities[:,2]**2)))
    forces = Trap.state.getForces(asNumpy = True)
    
    
    import time
    
    dt = Trap.timestep._value * 1e-9
    
    posxlist = []
    posylist = []
    poszlist = []
    populationslist = []
    velxlist = []
    velylist = []
    velzlist = []
    numphotlist = []
    pstateslist = []
    ekinlist = []
    epotlist = []
    etotlist = []
    detuninglist1 = []
    detuninglist2 = []
    
    
    t0 = time.time()
    
    for i in range(totalsteps):
        Trap.doStep(1)
        positions = Trap.positions
        velxlist.append(Trap.velocities[0,0] * 1e3)
        velylist.append(Trap.velocities[0,1] * 1e3)
        velzlist.append(Trap.velocities[0,2] * 1e3)
        populationslist.append(Trap.populations[0])
        pstateslist.append(Trap.populations[0,2] + Trap.populations[0,3])
    
        detuninglist1.append(Trap.d1s[0])
        detuninglist2.append(Trap.d2s[0])
    
    
        posxlist.append(positions[0,0]._value * 1e-9)
        posylist.append(positions[0,1]._value * 1e-9)
        poszlist.append(positions[0,2]._value * 1e-9)
    
        ekinlist.append(Trap.ekin)
        epotlist.append(Trap.epot)
        etotlist.append(Trap.etot)
    
    #    velocities = Trap.velocities
        numphotlist.append(Trap.numphotons)
    
    
    print('computation time of: %.2f'%(time.time()-t0))
    print('total number of photons', Trap.numphotons)
    
    period = 1/Trap.rffreq
    cyclesteps = int(np.floor(period / dt))
    length = len(etotlist)
    numcycles = int(np.floor(length / cyclesteps))
    cycleenergy = [np.mean(etotlist[i*cyclesteps:(i+1)*cyclesteps]) for i in range(numcycles)]
    
    np.save('./numphot', numphotlist)
    np.save('./xpos', posxlist)
    np.save('./ypos', posylist)
    np.save('./zpos', poszlist)
    np.save('./vx', velxlist)
    np.save('./vy', velylist)
    np.save('./vz', velzlist)
    np.save('./pops', populationslist)
    np.save('./pstates', pstateslist)
    np.save('./ekin', ekinlist)
    np.save('./etot', etotlist)
    np.save('./epot', epotlist)
    np.save('./etot_cycle', cycleenergy)
    np.save('./d1s', detuninglist1)
    np.save('./d2s', detuninglist2)
    np.save('./velocities_last', Trap.velocities)
    np.save('./positions_last', Trap.positions)
