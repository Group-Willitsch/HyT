#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 11:41:53 2019

@author: claudiovonplanta
"""

import scipy.constants as con
import numpy as np
import math
from matplotlib import pyplot as plt
from scipy.integrate import ode
from scipy.linalg import lu_factor, lu_solve


# translated to python by Claudio von Planta
# Refs:
#
# Matrices, equations
#
# Thesis Hilmar Oberst Appendix 8.3
# A. D. Gingell J. Chem. Phys. 133, 194302 (2010), https://doi.org/10.1063/1.3505142
# Ziv Matlab Code three_L_ziv_vec.m
#
# Values for Ca+
# M. S. Safronova1 and U. I. Safronova, Phys Rev A, 83, 012503 (2011)
# 
# LU decomposition
# https://stackoverflow.com/questions/40229147/why-do-mathematica-and-pythons-answers-differ-when-dealing-with-singular-matrix

# some constants
hbar = 1#con.hbar         # m**2 kg / s
hbar2 = 1#con.hbar 

muB = 9.27400994e-24    # J/T
amtoJ = 4.359744e-18    # conversion au to J
s3 = math.sqrt(3)       # sqrt(3)




def L_Matrix_oberst(u, alpha, beta, D1, D2, s1, s2, l1, l2):

    id8 = np.identity(8)
    v1, v2, v3, v4, v5, v6, v7, v8 = id8
    v1, v2, v3, v4, v5, v6, v7, v8 = v1.reshape(8,1), v2.reshape(8,1), v3.reshape(8,1), v4.reshape(8,1), v5.reshape(8,1), v6.reshape(8,1), v7.reshape(8,1), v8.reshape(8,1)
#
#    u = B * muB/hbar
#    g1 = 128/2/np.pi
#    g2 = 7.4/2/np.pi
    
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

def drho_f(t, rho, L):
    return np.dot(L, rho)

if __name__ == "__main__":
    
    plt.close('all')
    
    dynamics = False
    
    # settings for comparison with Oberst thesis
#    
#    alpha = 100/360*2*np.pi#np.pi/2
#    beta = alpha
#    s1 = 0.54
#    s2 = 3.9
#    l1 = 0.3
#    l2 = 0.3
#    u = 2.9
#    d2 = -10
#    d1 = -26
    
    
    alpha = np.pi/2#100/360*2*np.pi#np.pi/2
    beta = alpha
    s1 = 1
    s2 = 7
    l1 = 0.3
    l2 = 0.3
    B = 1e-4
    u =  B * muB/con.hbar / (2*np.pi) * 1e-6
    d2 = 0
    d1 = -20
    
    dt = 1e-13
    steps = 10#1e5
    t1 = steps*dt
    rho = np.zeros([64,1])
    rho[0] = 1
    L0 = L_Matrix_oberst(u, alpha, beta, 0, 0, s1, s2, l1, l2)
    dL1 = L_Matrix_oberst(u, alpha, beta, 1, 0, s1, s2, l1, l2)-L0
    dL2 = L_Matrix_oberst(u, alpha, beta, 0, 1, s1, s2, l1, l2)-L0
    Ln = L0 + d2 * dL2 + d1 * dL1
    
    # time dependence
    if dynamics:
        p0 = [[float(rho[i*8 + i].real) for i in range(8)]]
        
        integrator = ode(drho_f).set_integrator('zvode', method='bdf')
        integrator.set_initial_value(rho, 0).set_f_params(Ln)
        
        #while integrator.successful() and integrator.t < t1:
        for i in range(int(steps)):
        #    integrator.integrate(integrator.t + dt)
        #    rho = integrator.y
            rho = timestep(Ln, rho, dt)
            rho = np.round(rho, 11)
            p = [float(rho[i*8 + i].real) for i in range(8)]
            p0.append(p)
            rho = rho / np.sum(p)
            
        p0 = np.array(p0)    
        labels = [r'$\left| %i \right>$'%i for i in range(1,9)]
        xdat = np.array([i*dt for i in range(len(p0.T[0]))]) * 1e6
        
        fig, ax = plt.subplots()
        #plt.plot(np.sum(p0.T, axis = 0))
        for pop in range(len(p0.T)):
            plt.plot(xdat, p0.T[pop], label = labels[pop])
        
        plt.xlabel('time /us')
        plt.ylabel('population')
        plt.legend()
        fig.savefig('./dynamics.png')
    #
    #
    # steady state
    scanrange = 80#*(2*np.pi)
    freqs = np.linspace(-scanrange, scanrange, 1000)
    
    popslist = []
    
    
    
    L0 = L_Matrix_oberst(u, alpha, beta, 0, 0, s1, s2, l1, l2)
    dL1 = L_Matrix_oberst(u, alpha, beta, 1, 0, s1, s2, l1, l2)-L0
    dL2 = L_Matrix_oberst(u, alpha, beta, 0, 1, s1, s2, l1, l2)-L0
    
    for freq in freqs:
        L = L0+ freq * dL1 + d2 * dL2  
    #    L = L_Matrix(B, alpha, beta, d_b, d_r, O_b, O_r, g_b, g_r)
        pops = steadystate_populations(L)
        popslist.append(pops[2] + pops[3] )
        
    fig, ax = plt.subplots()
    plt.xlabel('frequency detuning /MHz')
    plt.ylabel(r'$4 p_{1/2}$ population')
    plt.ylim([0,1.1*np.max(popslist)])
    plt.plot(freqs/1e6, popslist)
    fig.savefig('./spectrum.png')
    
    L0 = L_Matrix_oberst(u, alpha, beta, 0, 0, s1, s2, l1, l2)
    dL1 = L_Matrix_oberst(u, alpha, beta, 1, 0, s1, s2, l1, l2)-L0
    dL2 = L_Matrix_oberst(u, alpha, beta, 0, 1, s1, s2, l1, l2)-L0
    Ln = L0 + d2 * dL2 + d1 * dL1
    #print(Ln)
    pops = steadystate_populations(Ln)
    print(pops)
