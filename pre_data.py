# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 13:45:34 2021

@author: astap
"""
import numpy as np


zmax    = 4000       #dybde i meter
n       = 128000      #antall målepunkter i dybden

H       =   5060 #M^-1 atm^-3. Proporsjonalitetskonstant mellom atmosfærisk CO2 og DIC_equilibrium-nivå
PCO2_0  =   415e-6 #Brukes ikke, men er start-konsentrasjonen for co2 i luften

a       =   6.97e-7 #mass-transfer-relationship konstant. Referert til i  oppgaveteori
u       =   10 #Antatt gjennomsnittlig vindhastighet over havene. Referert til i oppgaveteori
kw      =   a*u**2 #Mass-transfer proporsjonalitetskonstant ved overflate. Referert til i oppgaveteori

K0      =   1e-4 
K1      =   1e-2

dt      =   3600*4 #sekund mellom hver simulasjons-event  

C_mass  =   12 #Atomisk masse av karbon i g (per mol CO2, eller formulert annerledes; per mol C)
ocean_A =   360e12 #Havets overflateareal. Referert til i oppgaveteksten tilhørende problem 2.

z0 = 100; #m.     Konstant anvendt i oppgaveteoriens K(z)
a  = 0.5; #m^-1.  Konstant anvendt i oppgaveteoriens K(z)


zs, zstep = np.linspace(0, zmax, n, retstep=True) #henter ut alle steg zs i dybden, 
                                                  #Med tilhørende steglende zstep
                                                 

#K-verdier
Ks = K1+(K0-K1)/(1+np.exp(-a*(zs-z0)))     #Definerer K(z). Definisjon gitt i Problem 2 sin oppgavetekst.

#CO2-nivåer
atmospheric_co2 = lambda t: 415e-6 + 2.3e-6*t/(3600*24*365) #Uttrykk for CO2-konsentrasjon i atm
                                                            # i tid t (målt i sek). Satt til å være
                                                            # 415 ppm (t=0), med en økning på
                                                            # 2.3 ppm per år.

#Uttrykk for stabil marin DIC-konsentrasjon av tid. 
equi_concentration = lambda t: H*atmospheric_co2(t) #Uttrykk hentet fra Oppgaveteori. 
                                                    # er funksjon av t i sek.


t_keep_plot = np.array([0, 2.5, 5, 10]) #tidspunkter en ønsker separate plot (målt i år fra t=0)
t_keep_plot *= 365*3600*24 #konverterer fra år til sekunder

#liste av linje-type og farge for plots i økende rekkefølge. (første kurve plottes med punkter :
    #og er rød. Andre kurve punktes rødt med vekselvis linjer og punkter -. osv)
    #Rent kosmetisk variabel
linetype =     [(":", "r"), ("-.", "r"), ("--", "r"), ("-", "r"), 
               (":", "b"), ("-.", "b"), ("--", "b"), ("-", "b"),
               (":", "r"), ("-.", "r"), ("--", "r"), ("-", "r"),
               (":", "b"), ("-.", "b"), ("--", "b"), ("-", "b")]



modmax = 400 #antall kjøringer(tids-steg) før et plot blir laget/animert. 
plt_skiprate = 1 #En plotter ikke hele datakurver, men heller et utvalg datapunkter fra et datasett
                 # gitt med datapunkter[::plt_skiprate] for å gjøre plotting lettere på pc-en.




C = np.zeros_like(zs) #Skaper konsentrasjons-vektoren
C += equi_concentration(0) #setter start-konsentrasjon som dagens konsentrasjon (t=0, år 2020)


