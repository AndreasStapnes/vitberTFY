# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 16:38:29 2021

@author: astap
"""
import matplotlib.pyplot as plt
import numpy as np

import pre_data as pd

zs = pd.zs
C  = pd.C
zstep= pd.zstep
makeS= pd.makeS
develop = pd.develop
make_RL = pd.make_RL

equilib_conc = pd.concentration;


create_band = True





activedisplay=True
def handle_close(evt):
    global activedisplay
    activedisplay=False
    print("Activedisplay set to False")

fig, (concplot, minmax_plot) = plt.subplots(2,1, figsize=(10,10))
fig.canvas.mpl_connect("close_event", handle_close)
#fester her handle_close som kjøre-event når animasjons-figuren tvangs-lukkes
#dette er for å avslutte animasjon når dette skjer.






minmax_lvls = [];

Cs = [] #Lagrer konsentrasjoner for spesifikke tidspunkter    
ts = [] #Lagrer de nevnte tidspunktene
modmax = 100
plt_skiprate = 30


dagIS = 3600*24

t_show = 0;


def hexify_proportion(num, max_num):
    prop = 255*num//max_num
    hex_s = hex(prop)[2:].rjust(2, "0")
    return hex_s

def special_color(num, num_max):
    return "#" + "ff" + hexify_proportion(num_max-num, num_max)*2

for enum in range(10):
    dt = 1800 #sekund mellom hver simulasjons-event    
    t = 0;   modnum = 0;
    R,L = make_RL(dt, banded=create_band)
    
    while(t<dagIS*180*2/110*(enum+1) and activedisplay):
        S = makeS(dt, t_show)
        
        t+= dt;         t_show += dt;       modnum += 1
        
        C = develop(C,R,L,S, band_matrix=create_band)
        if(modnum % modmax == 0):
            minmax_lvls.append([np.min(C), np.max(C), equilib_conc(t_show), t_show])
            
            #plt.figure(0)
            concplot.cla()
            concplot.set_title("døgn = {:.2f}".format(t_show/60/60/24))
            
            Cslen = len(Cs)
            for enum_2, (t_choice, C_choice) in enumerate(zip(ts,Cs)):
                concplot.plot(zs[::plt_skiprate], C_choice[::plt_skiprate], 
                         label="{:.0f}D".format(t_choice/dagIS), 
                         color=special_color(enum_2+1, Cslen))
            concplot.plot(zs[::plt_skiprate], C[::plt_skiprate], label="{:.0f}D".format(t_show/dagIS), color="k")
            concplot.set_xlabel("dybde i m")
            concplot.legend()
            #plt.pause(1e-4)
            
            minmax_plot.cla()
            mins, maxs, eqs, tmms = np.array(minmax_lvls).T
            tmms /= (3600*24)
            minmax_plot.plot(tmms, mins, label="min", color="b")
            minmax_plot.plot(tmms, maxs, label="max", color="r")
            minmax_plot.plot(tmms, eqs, label="equilib", color="g")
            #minmax_plot.set_title("atmosfærisk co2 av tid")
            minmax_plot.legend()
            minmax_plot.set_xlabel("t i D")
            minmax_plot.set_ylabel("DIC")
            plt.pause(1e-4)
            
            #plt.pause(1e-4)
    Cs.append(C)
    ts.append(t_show)
    
concplot.cla()
concplot.set_title("døgn = {:.2f}".format(t_show/dagIS))
Cslen = len(Cs)
for enum_2, (t_choice, C_choice) in enumerate(zip(ts,Cs)):
    concplot.plot(zs, C_choice, 
             label="{:.0f}D".format(t_choice/dagIS),
             color = special_color(enum_2+1, Cslen))
concplot.legend()
    
