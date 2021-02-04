# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 13:45:34 2021

@author: astap
"""
import numpy as np


zmax    = 100       #dybde i meter
n       = 1000      #antall målepunkter i dybden

H       =   5060 
PCO2_0  =   415e-6 #Brukes ikke, men er start-konsentrasjonen for co2 i luften
#conc   =   H*PCO2 = Stabil CO2 konsentrasjon i vannet

a       =   6.97e-7
u       =   10
kw      =   a*u**2

K0      =   1e-3; 
K1      =   2e-2;  #Kalles Ka i teksten
K2      =   5e-2;  #Kalles Kb i teksten

za = 7; zb = 10;



zs, zstep = np.linspace(0, zmax, n, retstep=True)

Ks = K0 + K1*zs/za*np.exp(-zs/za) + K2*(zmax-zs)/zb*np.exp(-(zmax-zs)/zb)
#K-verdier


atmospheric_co2 = lambda t: 415e-6*(np.exp(t/(3600*24*365)))
concentration = lambda t: H*atmospheric_co2(t)


#Hjelpestørrelser:
dagIS = 3600*24


#størrelser en ikke burde endre på

C = np.zeros_like(zs)
#C += concentration(0) setter start-konsentrasjon som dagens konsentrasjon


DK = np.concatenate((np.array([0]),  #kalles K' i teksten
                     Ks[2:] - Ks[:-2], 
                     np.array([0])))




def makeS(delta_t, t = 0):
    alpha = delta_t/(2*zstep**2)
    gamma = 2*alpha*kw*zstep*(1-(Ks[1]-Ks[0])/(2*Ks[0]))
    S = np.zeros(n)
    S[0] = 2*gamma*(concentration(t));
    return S


import scipy.linalg as la
def develop(C, R, L, S, band_matrix=False):   
    if not band_matrix:
        return la.solve(L, R@C + S)
    else:
        C = la.solve_banded((1,1), L, (R@C) + S)
    return C


def make_RL(delta_t, banded=False):
    #banded gjør matrisen kortere (lagrer kun diagonaler).
    #Med spesialiserte funksjoner
    #(linalg.solve_banded) kan matriselikninger
    #på formen Ax = b for banded A løses mye raskere.
    
    #skaper R og L
    diag = np.identity(n)
    u_diag = np.roll(diag, 1, 1)
    u_diag[n-1,0]=0 #upper diagonal
    l_diag = np.roll(np.roll(u_diag, -1, 1), 1, 0) #lower diagonal
    ul_correction = np.zeros(np.shape(diag)); lr_correction = np.zeros(np.shape(diag)) #bonus-element i posisjon [N,N-1] for å kompensere for at [N,N+1] ikke eksisterer 
    lr_correction[n-1,n-2] = 1; ul_correction[0,1] = 1;# tilsvarende forrige, men kompenserer for at element [-1,0] ikke eksisterer ved bonus i [1,0]
    
    alpha = delta_t/(2*zstep**2)
    gamma = 2*alpha*kw*zstep*(1-(Ks[1]-Ks[0])/(2*Ks[0]))
    
    L = -((-2*diag + u_diag + l_diag + ul_correction.T+lr_correction.T) * Ks).T * alpha
    L += (((u_diag-ul_correction - l_diag+lr_correction) * alpha/4).T*DK).T
    L[0,0] += gamma;
    
    R = diag - L
    L = diag + L
    
    if(banded):
        #Skaper banded (diagonal-representasjon) L matrise
        B = np.zeros((3,n))
        B[0,1:] = L.diagonal(1);
        B[1,:] = L.diagonal(0);
        B[2,:-1] = L.diagonal(-1);
        L = B
    return R,L