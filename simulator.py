# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 16:38:29 2021

@author: astap
"""
import matplotlib.pyplot as plt
import numpy as np

zmax = 100; n=600

H = 5080e-3
PCO2 = 415e-6
conc = H*PCO2 

kw = 6.97e-7*10**2

K0 = 1e-3; K1 = 2e-2; K2 = 5e-2;
za = 7; zb = 10;
Kappa = lambda z: K0 + K1*z/za*np.exp(-z/za) + K2*(zmax-z)/zb*np.exp(-(zmax-z)/zb)
zs, zstep = np.linspace(0, zmax, n, retstep=True)
#zs = np.concatenate((zs, np.array([zmax])));
C = np.zeros_like(zs)
Ks = Kappa(zs);
DK = np.concatenate((np.array([0]), (Ks[2:] - Ks[:-2]), np.array([0])))


t = 0
modnum = 0



activedisplay=True
def handle_close(evt):
    global activedisplay
    activedisplay=False
    print("Activedisplay set to False")

fig = plt.figure(0)
fig.canvas.mpl_connect("close_event", handle_close)

def makeS(delta_t, zstep, t = 0):
    alpha = delta_t/(2*zstep**2)
    gamma = 2*alpha*kw*zstep*(1-(Ks[1]-Ks[0])/(2*Ks[0]))
    S = np.zeros(n)
    S[0] = 2*gamma*(conc*(np.cos(t/(3600*24*180/(2*np.pi)))));
    return S


def make_RL(delta_t, Ks, DK, zstep):

    diag = np.identity(n)
    u_diag = np.roll(diag, 1, 1)
    u_diag[n-1,0]=0
    l_diag = np.roll(np.roll(u_diag, -1, 1), 1, 0)
    ul_correction = np.zeros(np.shape(diag)); lr_correction = np.zeros(np.shape(diag))
    lr_correction[n-1,n-2] = 1; ul_correction[0,1] = 1;
    alpha = delta_t/(2*zstep**2)
    gamma = 2*alpha*kw*zstep*(1-(Ks[1]-Ks[0])/(2*Ks[0]))
    
    L = -((-2*diag + u_diag + l_diag + ul_correction.T+lr_correction.T) * Ks).T * alpha
    L += (((u_diag-ul_correction - l_diag+lr_correction) * alpha/4).T*DK).T
    L[0,0] += gamma;
    
    R = diag - L
    L = diag + L
    
    return R,L



create_band = True

import scipy.linalg as la
import scipy
def develop(C, R, L, S, band_matrix=False):   
    if not band_matrix:
        return la.solve(L, R@C + S)
    else:
        C = la.solve_banded((1,1), B, (R@C) + S)
    return C


Cs = []    
ts = []

t = 0; modnum = 0
    
from numba import jit
import scipy


plt.figure(1)
plt.plot(zs, Ks)

plt.figure(0)



C = np.zeros(n); t_show = 0;
for enum in range(10):
    dt = 1800
    
    t = 0; 
    #C = np.zeros(n)
    R,L = make_RL(dt, Ks, DK, zstep)
    
    if(create_band):
        #Skaper diagonal L matrise
        n = np.shape(L)[0]
        B = np.zeros((3,n))
        #B = [[*L.diagonal(1)], [*L.diagonal(0)], [*L.diagonal(-1)]]
        B[0,1:] = L.diagonal(1);
        B[1,:] = L.diagonal(0);
        B[2,:-1] = L.diagonal(-1);
        L = B
    
    while(t<60*60*24*180*2/110*(enum+1) and activedisplay):
        S = makeS(dt, zstep, t_show)
        
        modnum += 1
        t+= dt; t_show += dt
        
        C = develop(C,R,L,S, create_band)
        if(modnum % 100 == 0):
            fig.clf()
            plt.figure(0)
            plt.title("døgn = {:.2f}".format(t_show/60/60/24))
            for t_choice, C_choice in zip(ts,Cs):
                plt.plot(zs, C_choice, label="dt:{:.1e}s".format(t_choice))
            plt.plot(zs, C, label="{:.1e}".format(t_show))
            #plt.ylim([0,conc])
            plt.legend()
            plt.draw()
            plt.pause(1e-3)
    Cs.append(C)
    ts.append(t_show)
    
fig.clf()
plt.figure(0)
plt.title("døgn = {:.2f}".format(t_show/60/60/24))
for t_choice, C_choice in zip(ts,Cs):
    plt.plot(zs, C_choice, label="dt:{:.1e}s".format(t_choice))
plt.legend()
    #print("*", end="")
#break
    
