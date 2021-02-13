# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 16:38:29 2021

@author: astap
"""
import matplotlib.pyplot as plt
import numpy as np

from pre_data import dt, C, t_keep_plot, equi_concentration, zs, linetype, plt_skiprate, modmax, kw, zstep


from functionals import makeS, develop, make_RL, special_color

from scipy.integrate import simps



plot = True



activedisplay=True
def handle_close(evt):
    global activedisplay
    activedisplay=False
    print("Activedisplay set to False")
#fig, (concplot, minmax_plot) = plt.subplots(2,1, figsize=(10,10))
fig1, concplot = plt.subplots(1,1,figsize=(10,5))
fig1.canvas.mpl_connect("close_event", handle_close)
#fester her handle_close som kjøre-event når animasjons-figuren tvangs-lukkes
#dette er for å avslutte animasjon når dette skjer.



minmax_lvls = []; #Lagrer minimums og maximums-konsentrasjoner i havet for skjellige tider t_keep_plot
Cs = [] #Lagrer konsentrasjoner for spesifikke tidspunkter    
masses = []

dagIS = 3600*24
t = 0;
mass = simps(C)*zstep;

current_mass_est = 0
R,L = make_RL(dt)
for t_stop in t_keep_plot:
    
    modnum = 0;
    
    while(t<t_stop and activedisplay):
        S = makeS(dt, t)
        
        masses.append([t,mass, simps(C)*zstep]);
        minmax_lvls.append([np.min(C), np.max(C), equi_concentration(t), t])
        
        mass += kw*(1/2*(equi_concentration(t)+equi_concentration(t+dt))-C[0])*dt 
    
        t += dt;       modnum += 1
        C = develop(C,R,L,S)
        
        
        if(modnum % modmax == 0 and plot):
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
    Cs.append(C)
    
concplot.cla()
concplot.set_title("døgn = {:.2f}".format(t/dagIS))
Cslen = len(Cs)
for enum_2, (t_choice, C_choice) in enumerate(zip(t_keep_plot,Cs)):
    concplot.plot(zs, C_choice, 
             label="{:.0f} døgn".format(t_choice/dagIS),
             color=special_color(enum_2+1, Cslen, linetype[enum_2][1]),
             linestyle=linetype[enum_2][0])
    concplot.set_title(r"DIC konsentrasjoner av dybde $z$ ved forskjellig $t$")
    concplot.set_xlabel("dybde i m")
    concplot.set_ylabel(r"DIC i $\frac{mol}{m^3}$")
concplot.legend(loc="upper right")

fig2, minmax_plot = plt.subplots(1,1,figsize=(10,5))
minmax_plot.cla()
mins, maxs, eqs, tmms = np.array(minmax_lvls).T
tmms /= (3600*24)
minmax_plot.plot(tmms, mins, label="min", color="b")
minmax_plot.plot(tmms, maxs, label="max", color="r")
minmax_plot.plot(tmms, eqs, label="equilibrium", color="g")
#minmax_plot.set_title("atmosfærisk co2 av tid")
minmax_plot.legend()
minmax_plot.set_title(r"min- og max- konsentrasjoner av $t$")
minmax_plot.set_xlabel("t i Døgn")
minmax_plot.set_ylabel(r"DIC i $\frac{mol}{m^3}$")


fig3, massplot = plt.subplots(1,1,figsize=(10,5))
massplot.cla()
masses = np.array(masses).T
massplot.plot(masses[0]/(3600*24), masses[1], label="Beregnet DIC fra overflatekondisjoner")
massplot.plot(masses[0]/(3600*24), masses[2], label="Integrert DIC")
#massplot.plot(masses[0]/(3600*24), masses[2]-masses[1])
massplot.set_title(r"totalt absorbert CO$_2$ av $t$")
massplot.set_ylabel(r"DIC i $\frac{mol}{m^2}$")
massplot.set_xlabel("t i døgn")
massplot.legend()


fig1.suptitle(r"Tidsutvikling med økende atmosfærisk CO$_2$-nivå", fontsize=16)
fig2.suptitle(r"Tidsutvikling med økende atmosfærisk CO$_2$-nivå", fontsize=16)

fig1.savefig("concplot.pdf")
fig2.savefig("minmaxplot.pdf")
fig3.savefig("massplot.pdf")
    
