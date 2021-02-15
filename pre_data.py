# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 13:45:34 2021

@author: astap
"""
import numpy as np


zmax    = 4000       #dybde i meter
n       = 128000      #antall målepunkter i dybden

H       =   5060 
PCO2_0  =   415e-6 #Brukes ikke, men er start-konsentrasjonen for co2 i luften
#conc   =   H*PCO2 #= Stabil CO2 konsentrasjon i vannet

a       =   6.97e-7
u       =   10
kw      =   a*u**2 

K0      =   1e-4 
K1      =   1e-2

dt      =   3600*4 #sekund mellom hver simulasjons-event  

C_mass  =   12
ocean_A =   360e12

z0 = 100;
a  = 0.5;


zs, zstep = np.linspace(0, zmax, n, retstep=True)

Ks = K1+(K0-K1)/(1+np.exp(-a*(zs-z0)))
#K-verdier

atmospheric_co2 = lambda t: 415e-6 + 2.3e-6*t/(3600*24*365)
equi_concentration = lambda t: H*atmospheric_co2(t) #Stabil overflatekonsentrasjon


t_keep_plot = np.array([0, 2.5, 5, 10])*365 #Lagrer de nevnte tidspunktene (i år)
t_keep_plot *= 3600*24 #konverterer til sekunder

linetype =     [(":", "r"), ("-.", "r"), ("--", "r"), ("-", "r"), 
               (":", "b"), ("-.", "b"), ("--", "b"), ("-", "b"),
               (":", "r"), ("-.", "r"), ("--", "r"), ("-", "r"),
               (":", "b"), ("-.", "b"), ("--", "b"), ("-", "b")]



modmax = 400 #antall kjøringer før et plot blir laget
plt_skiprate = 1 #skip-rate under kjøring i plot




C = np.zeros_like(zs)
C += equi_concentration(0) #setter start-konsentrasjon som dagens konsentrasjon


