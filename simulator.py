# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 16:38:29 2021

@author: astap
"""
import matplotlib.pyplot as plt
import numpy as np

from pre_data import dt, C, t_keep_plot, equi_concentration, zs, linetype, plt_skiprate, modmax, kw, zstep, ocean_A, C_mass


from functionals import makeS, develop, make_RL, special_color

from scipy.integrate import simps



plot = True #Bestemmer om koden skal returnere feedback-animasjon under kjøring
#For å få brukbare animasjoner, husk å ha graphic-settings satt til å plotte i 
#separat pop-up vindu.



activedisplay=True
def handle_close(evt): #Handle-funksjon for å tvangs-stoppe kjøring når man trykker kryss på en graf
    global activedisplay
    activedisplay=False
    print("Activedisplay set to False")

fig1, concplot = plt.subplots(1,1,figsize=(10,5))
fig1.canvas.mpl_connect("close_event", handle_close)



minmax_lvls = []; #Lagrer minimums og maximums-konsentrasjoner i havet for skjellige tider t_keep_plot
Cs = [] #Lagrer konsentrasjons-vektorer for spesifikke tidspunkter (angitt i pre_data)   
masses = [] #Lagrer absorbert DIC ved forskjellige tidspunkt

dagIS = 3600*24 #sekunder i en dag
t = 0; #Start-tid
mass = simps(C)*zstep; #Start-masse for havene

R,L = make_RL(dt) #Skaper R og L matrisene

for t_stop in t_keep_plot: #Kjører gjennom tidspunktene angitt i pre_data
    
    modnum = 0; #Brukt for å begrense plotting til å ikke skje for hvert simulasjons-steg
    
    while(t<t_stop and activedisplay):
        S = makeS(dt, t) #Skaper S-vektoren ved dette tidspunktet
        
        masses.append([t,mass, simps(C)*zstep]); #Legger til tid, masse (beregnet fra overflate-kondisjoner)
                                                 # og masse (integrert med simps) i masses
                                                 
        minmax_lvls.append([np.min(C), np.max(C), equi_concentration(t), t]) #Legger til min og max kons
                                                                             #samt equilibrium kons og tid
                                                                             #i minmax_lvls
        
        mass += kw*(1/2*(equi_concentration(t)+equi_concentration(t+dt))-C[0])*dt #Beregne DIC-økning pga.
                                                                                  #absorbsjon ved overflate
    
        t += dt;       modnum += 1 #Gå til neste tidssteg
        C = develop(C,R,L,S) #Føre simulasjonen videre
        
        
        if(modnum % modmax == 0 and plot): #Plotte-funksjoner for animasjonen skjer her
            concplot.cla()
            concplot.set_title("døgn = {:.2f}".format(t/60/60/24))
            Cslen = len(Cs)
            for enum_2, (t_choice, C_choice) in enumerate(zip(t_keep_plot,Cs)):
                concplot.plot(zs[::plt_skiprate], C_choice[::plt_skiprate], 
                         label="{:.0f}D".format(t_choice/dagIS), 
                         color=special_color(enum_2+1, Cslen, linetype[enum_2][1]),
                         linestyle=linetype[enum_2][0])
            concplot.plot(zs[::plt_skiprate], C[::plt_skiprate], label="{:.0f}D".format(t/dagIS), color="k")
            concplot.set_title(r"DIC konsentrasjoner av dybde $z$ ved forskjellig $t$")
            concplot.set_xlabel("dybde i m")
            concplot.set_ylabel(r"DIC i $\frac{mol}{m^3}$")
            concplot.legend()
            plt.pause(1e-4)
    Cs.append(C) #Legg til C i Cs for hver gang t passerer et tidspunkt i t_keep_plot fra pre_data
    
#Plot konsentrasjonene lagret i Cs
concplot.cla() #clear concplot
concplot.set_title("døgn = {:.2f}".format(t/dagIS))
Cslen = len(Cs)
for enum_2, (t_choice, C_choice) in enumerate(zip(t_keep_plot,Cs)):
    concplot.plot(zs, C_choice, 
             label="{:.0f} døgn".format(t_choice/dagIS),
             color=special_color(enum_2+1, Cslen, linetype[enum_2][1]),
             linestyle=linetype[enum_2][0]) #plotter hver konsentrasjon
    concplot.set_title(r"DIC konsentrasjoner av dybde $z$ ved forskjellig $t$")
    concplot.set_xlabel("dybde i m")
    concplot.set_ylabel(r"DIC i $\frac{mol}{m^3}$")
concplot.legend(loc="upper right")

#Plotter minimum og maximums-konsentrasjoner for havene
fig2, minmax_plot = plt.subplots(1,1,figsize=(10,5))
minmax_plot.cla()
mins, maxs, eqs, tmms = np.array(minmax_lvls).T #tmms  er tidspunkter. De andre er selvforklarende
tmms /= (3600*24) #Konverterer tmms til døgn
minmax_plot.plot(tmms, mins, label="min", color="b")
minmax_plot.plot(tmms, maxs, label="max", color="r")
minmax_plot.plot(tmms, eqs, label="equilibrium", color="g")
minmax_plot.legend()
minmax_plot.set_title(r"min- og max- konsentrasjoner av $t$")
minmax_plot.set_xlabel("t i Døgn")
minmax_plot.set_ylabel(r"DIC i $\frac{mol}{m^3}$")


#Plotter havets absorberte masse
fig3, massplot = plt.subplots(1,1,figsize=(10,5))
massplot.cla()
masses = np.array(masses).T #Her vil masses[0] være tider, masses[1] være beregnet fra overflatekond
                            #og masses[2] vil være beregnet fra integrasjon
masses[1:] *= C_mass * ocean_A   #Konverterer kolonne-DIC til globalt DIC i gram basert på parametre fra pre_data
massplot.plot(masses[0]/(3600*24), masses[1], label="Beregnet DIC fra overflatekondisjoner")
massplot.plot(masses[0]/(3600*24), masses[2], label="Integrert DIC")
massplot.set_title(r"globalt absorbert CO$_2$ i verdenshavene")
massplot.set_ylabel(r"DIC i $g$")
massplot.set_xlabel("t i døgn")
massplot.legend()
print("Start-masse: Integrert={}, overflate={}".format(masses[2][0], masses[1][0]))
print("Slutt-masse: Integrert={}, overflate={}".format(masses[2][-1], masses[1][-1]))



fig1.suptitle(r"Tidsutvikling med økende atmosfærisk CO$_2$-nivå", fontsize=16)
fig2.suptitle(r"Tidsutvikling med økende atmosfærisk CO$_2$-nivå", fontsize=16)

fig1.savefig("concplot.pdf")
fig2.savefig("minmaxplot.pdf")
fig3.savefig("massplot.pdf")
