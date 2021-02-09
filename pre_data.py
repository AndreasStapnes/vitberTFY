# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 13:45:34 2021

@author: astap
"""
import numpy as np


zmax    = 100       #dybde i meter
n       = 1600      #antall målepunkter i dybden

H       =   5060 
PCO2_0  =   415e-6 #Brukes ikke, men er start-konsentrasjonen for co2 i luften
#conc   =   H*PCO2 #= Stabil CO2 konsentrasjon i vannet

a       =   6.97e-7
u       =   10
kw      =   a*u**2 

K0      =   1e-3 
K1      =   2e-2  #Kalles Ka i teksten
K2      =   5e-2  #Kalles Kb i teksten

dt      =   3600 #sekund mellom hver simulasjons-event  

za = 7; zb = 10;



zs, zstep = np.linspace(0, zmax, n, retstep=True)

Ks = K0 + K1*zs/za*np.exp(-zs/za) + K2*(zmax-zs)/zb*np.exp(-(zmax-zs)/zb)
#K-verdier

atmospheric_co2 = lambda t: 415e-6 + 2e-6*t/(24*365*3600)
equi_concentration = lambda t: H*atmospheric_co2(t) #Stabil overflatekonsentrasjon


t_keep_plot = np.array([0, 2, 7, 15, 24, 44, 100, 180]) #Lagrer de nevnte tidspunktene (i døgn)
t_keep_plot *= 3600*24 #konverterer til sekunder


linetype =     [(":", "r"), ("-.", "r"), ("--", "r"), ("-", "r"), 
               (":", "b"), ("-.", "b"), ("--", "b"), ("-", "b"),
               (":", "r"), ("-.", "r"), ("--", "r"), ("-", "r"),
               (":", "b"), ("-.", "b"), ("--", "b"), ("-", "b")]



modmax = 40 #antall kjøringer før et plot blir laget
plt_skiprate = 1 #skip-rate under kjøring i plot




C = np.zeros_like(zs)
C += equi_concentration(0);
#C += concentration(0) #setter start-konsentrasjon som dagens konsentrasjon


