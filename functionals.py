# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 09:33:57 2021

@author: astap
"""
from pre_data import Ks, zstep, kw, n, equi_concentration
import numpy as np
from numba import jit
import scipy.sparse as sparse


#kalles K' (K-prime) i teksten
DK = np.concatenate((np.array([0]),  
                     Ks[2:] - Ks[:-2], 
                     np.array([0])))
#Henter fra oppgaveteori. Begynner og slutter med et null-element for å gjøre den til samme
#dimensjon som konsentrasjonsvektoren (matrisemultiplikasjon blir lettere slik)


@jit(nopython = True)
def tdma_solver(a, b, c, d): #Tridiagonal løsningsmetode hentet fra oppgaveteori.
    N = len(d)
    c_ = np.zeros(N-1)
    d_ = np.zeros(N)
    x  = np.zeros(N)
    c_[0] = c[0]/b[0]
    d_[0] = d[0]/b[0]
    for i in range(1, N-1):
        c_[i] = c[i]/(b[i] - a[i-1]*c_[i-1])
    for i in range(1, N):
        d_[i] = (d[i] - a[i-1]*d_[i-1])/(b[i] - a[i-1]*c_[i-1])
    x[-1] = d_[-1]
    for i in range(N-2, -1, -1):
        x[i] = d_[i] - c_[i]*x[i+1]
    return x

def tdma(A, b): #Tridiagonal løser som tar inn en array av form 
                # [Øvre diagonal, Midtre diagonal, Nedre diagonal]
                #som representerer en tridiagonal matrise A for hvilket man vil løse Ax = b for x.
    x = tdma_solver(A[2], A[1], A[0], b)
    return x


def makeS(delta_t, t = 0): #En funksjon hvilket genererer S ved en gitt tid t, slik definert i oppgaveteori
    alpha = delta_t/(2*zstep**2)
    gamma = 2*alpha*kw*zstep*(1-(Ks[1]-Ks[0])/(2*Ks[0]))
    S = np.zeros(n)
    S[0] = gamma*(equi_concentration(t) + equi_concentration(t+delta_t));
    #Legger til 2*gamma*(equi-conc_1, equi-conc_2 sitt gjennomsnitt)
    #For å kombinere det oppgaven refererer til som 1/2(S^i + S^{i+1}) til en "S"
    return S


import scipy.linalg as la
def develop(C, R, L, S):   
    #Kjører ett tids-steg, altså, beregner konsentrasjonsvektor C^{i+1} fra
    # konsentrasjonsvektor C^i = C. Metoden og matrisene er tilsvarende dem nevnt i oppgaveteorien.
    C = tdma(L, R@C + S)
    return C


def make_RL(delta_t):
    #banded gjør matrisen kortere (lagrer kun diagonaler).
    #Med spesialiserte funksjoner
    #(linalg.solve_banded) kan matriselikninger
    #på formen Ax = b for banded A løses mye raskere.
    
    #skaper R og L
    '''
    diag = np.identity(n)
    u_diag = np.roll(diag, 1, 1)
    u_diag[n-1,0]=0 #upper diagonal
    l_diag = np.roll(np.roll(u_diag, -1, 1), 1, 0) #lower diagonal
    ul_correction = np.zeros(np.shape(diag)); lr_correction = np.zeros(np.shape(diag)) #bonus-element i posisjon [N,N-1] for å kompensere for at [N,N+1] ikke eksisterer 
    lr_correction[n-1,n-2] = 1; ul_correction[0,1] = 1;# tilsvarende forrige, men kompenserer for at element [-1,0] ikke eksisterer ved bonus i [1,0]
    
    alpha = delta_t/(2*zstep**2)
    gamma = 2*alpha*kw*zstep*(1-(Ks[1]-Ks[0])/(2*Ks[0]))
    
    L = -((-2*diag + u_diag + l_diag + ul_correction.T+lr_correction.T) * Ks).T * alpha
    L += (((-u_diag+ul_correction + l_diag-lr_correction) * alpha/4).T*DK).T
    L[0,0] += gamma;
    '''
    #Ny diagonal konstruksjon-----------------------------
    tridiag = (1,0,-1) #Viser til de tre diagonalene i diagonal-matrisene vi bryr oss med
    null = sparse.diags((0,0,0), tridiag, shape=(n,n)) #skaper tridiagonal nullmatrise
    diag = sparse.diags((0,1,0), tridiag, shape=(n,n)) #Identitetsmatrisen, tridiagonal
    u_diag = sparse.diags((1,0,0),tridiag, shape=(n,n)) #Upper identitetsmatrise, tridiagonal
    l_diag = sparse.diags((0,0,1), tridiag, shape=(n,n)) #Bedre identitetsmatrise, tridiagonal
    null_diags = [null.diagonal(1), null.diagonal(0), null.diagonal(-1)] #Henter ut diagonalene i nullmatrisen for å kopiere disse
    ul_corr = [i.copy() for i in null_diags]; ul_corr[0][0] = 1; ul_corr = sparse.diags(ul_corr, tridiag)
    lr_corr = [i.copy() for i in null_diags]; lr_corr[2][-1] = 1; lr_corr = sparse.diags(lr_corr, tridiag)
    #ul_corr = upper-left-correction er en tom matrise unntatt element [0,1], som har verdi 1
    #lr_corr = lower-right-correction er en tom matrise unntatt element [n-1,n-2], som har verdi 1
    
    
    alpha = delta_t/(2*zstep**2)
    gamma = 2*alpha*kw*zstep*(1-(Ks[1]-Ks[0])/(2*Ks[0]))
    
    L = -((-2*diag + u_diag + l_diag + ul_corr.T+lr_corr.T).multiply(Ks)).T * alpha
    #Skaper alle ledd i elementene i matrisen L hvor det inngår alpha multiplisert med en K
    
    L += (((-u_diag+ul_corr + l_diag-lr_corr) * alpha/4).T.multiply(DK)).T
    #Skaper alle ledd i elementene i matrisen L hvor det inngår alpha/4*K'
    
    L[0,0] += gamma; #Legger til gamma-leddet i begynnelsen
    #--------------------------------------------------
    
    
    R = diag - L #Legger til 1-er leddene på diagonalene. Skaper L og R.
    L = diag + L
    
    L = [L.diagonal(1), L.diagonal(0), L.diagonal(-1)] #Lagrer L som sine diagonaler, slik etterspurt av TDMA-solver
    R = sparse.diags((R.diagonal(1), R.diagonal(0), R.diagonal(-1)), (1,0,-1)) 
        #Lagrer R som sparse tridiagonal for rask matrisemultiplikasjon
    return R,L




def hexify_proportion(num, max_num): #Uviktig kosmetisk funksjon for å skape heksadesimale tall
    prop = 255*num//max_num
    hex_s = hex(prop)[2:].rjust(2, "0")
    return hex_s

def special_color(num, num_max, color_type="r"): #Uviktig kosmetisk funksjon for å gi farge på kurver
    if(color_type=="r"):
        return "#" + "ff" + hexify_proportion(num_max-num, num_max)*2
    elif(color_type=="b"):
        return "#" + hexify_proportion(num_max-num, num_max)*2 + "ff"
