# -*- coding: utf-8 -*-
"""

Created on Tue Sep 14 02:12:56 2021

@author: Travis J Czechorski

Github: https://github.com/tjczec01

Email:   travis.czechorski@maine.edu
         tjczec01@gmail.com
        
Advisor: thomas.schwartz@maine.edu
       
Github:  https://github.com/tjczec01 
         https://github.com/TravisCzechorskiUMaine

"""

# from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
import numpy as np
import scipy.optimize as optimize
import matplotlib.pylab as plt
import math
import pprint as pp
import thermo as tc
from tqdm import tqdm
from scipy.integrate import solve_ivp

# for i in range(1, 7, 1):
#     print("dh{} = \nds{} = \ndg{} = ".format(i, i, i))

def kst(dh, ds, dg, dht, dst, dgt, T, R):
    kB = 1.380662 * 10E23   #  m2 kg/s2 K1
    h = 6.626176 * 10E34  #  m2 kg/s
    R = 1.987  #   cal/K mol
    Rk = 1.987  / 1000    #   kcal/K mol
    k_constant = (kB*T)/h
    Keq = np.exp(ds/R) * np.exp((-1*dh)/(R * T))
    Kf = k_constant*np.exp(dst/R) * np.exp((-1*dht)/(R * T))
    Kb = Kf/Keq
    return Kf, Kb, Keq

def TOFo(T, δE):
    #  δE = kcal/mol, T = K
    R = 1.987  #   cal/K mol
    Rk = 1.987  / 1000    #   kcal/K mol
    kB = 1.380662 * 10E23   #  m2 kg/s2 K1
    h = 6.626176 * 10E34  #  m2 kg/s
    K_term = (kB*T)/h
    E_term = δE/(Rk * T)
    finalf = K_term*np.exp(E_term)
    return finalf
    

def TOF(T, T0, δE, TOF0):
    #  δE = kcal/mol, T = K
    R = 1.987  #   cal/K mol
    Rk = 1.987  / 1000    #   kcal/K mol
    kB = 1.380662 * 10E23   #  m2 kg/s2 K1
    h = 6.626176 * 10E34  #  m2 kg/s
    # K_term = (kB*T)/h
    T_term0 = T/T0
    T_term = (1/T) - (1/T0)
    E_term = δE/Rk
    EXP_term = np.exp(E_term*T_term)
    finalf2 = TOF0*T_term*np.exp(E_term)
    return finalf2

def sticking(s0, T, T0, Ea, R):
    tt = 1/T - 1/T0
    er = -Ea/R
    expt = np.exp(er*tt)
    fv = expt*s0
    return fv

def ftheata(thetar):
    return [(1 - thetar), (1 - thetar)**2]


def kcols(s0, thetar, T, T0, P, σ, ν, EaA, EaD, Ma,  nn=1):
    #  Θr = Θ/Θsat
    kB = 1.380662 * 10E23   #  m2 kg/s2 K1
    pi = np.pi
    St = sticking(s0, T, T0, EaA, R=8.31446261815324)
    Tr = ftheata(thetar)[0]
    denom = 2.0 * pi * Ma * kB * T
    df = np.sqrt(denom)
    numer = St*Tr
    kf = numer / df
    kb = -1*ν*(σ**nn)*(thetar**nn)
    return kf, kb
    
def RHS(t, y, K_1f, K_1b, K_2f, K_2b, K_3tf, K_3tb, K_3f, K_3b, K_4f, K_4b, K_5tf, K_5tb, K_5f, K_5b, K_6f, K_6b, s0, thetar, σ, ν, δE, TOF0, EaA, EaD, T0):
    ETOH0, ace0, ety0, h220, h200, C_etoh, C_ace, C_ety, C_h2, C_water, θetoh, θeto, θh, θet, θace, θempty, θwatet, θh20, Petoh, Pace, Ph20, T = y
    R1 = K_1f*(Petoh**-1)*(θempty**-1) - K_1b*(Petoh)*(θetoh) 
    R2 = K_2f*(Petoh**-1)*(θetoh**-1) - K_2b*((Petoh)*(θeto) * (Petoh)*(θh))
    R3t = K_3tf*((Petoh**-1)*(θeto**-1) * (Petoh**-1)*(θh**-1)) - K_3tb*(Petoh)*(θet) 
    R3 = K_3f*Petoh*θeto
    R4 = K_4f*(Pace**-1)*(θace**-1) - K_4b*(Pace)*(θempty)
    R5t = K_5tf*(Petoh**-1)*(θetoh**-1) - K_5tb*(Petoh)*(θwatet)
    R5 = K_5f*(Petoh)*(θwatet)
    R6 = K_6f*(Ph20**-1)*(θh20**-1) * (Ph20**-1)*(θh**-1) - K_6b*(Ph20)*(θempty) 
    dETOHDt = ETOH0/t + C_etoh/t - R1
    daceDt = ace0/t + C_ace/t + R4
    detyDt = ety0/t + C_ety/t + R5
    dh2Dt = h220/t + C_h2/t + R3
    dh20Dt = h200/t + C_water/t + R6
    DθetohDt = R1 + R2 - R5t
    DθetoDt = R2 - R3t
    DθhDt = R2 + R5 - R3t - R6
    DθetDt = R3t - R3
    DθaceDt = R3 - R4
    Dθ5tDt = R5t - R5
    DθohDt = R5 - R6
    return [dETOHDt, daceDt, daceDt, detyDt, dh2Dt, dh20Dt, DθetohDt, DθetoDt, DθhDt, DθetDt, DθaceDt, Dθ5tDt, DθohDt]


sf = 0.016
Mg = 24.305
Si = 28.0855
O = 15.999
MgSi = Mg + Si
gcat = Mg + Si + 3*O
mmetal = 0.0034966141202839003

dh1 = -29.9
ds1 = -26.8
dg1 = -10.5
dh2 = 5.1
ds2 =  9.9
dg2 = -3.0
dh3 = 30.8
ds3 = 15.9
dg3 = 19.1
dh3t = 51.4
ds3t = 15.9
dg3t = 39.7
dh4 = 9.3
ds4 = 24.5
dg4 = -8.3
dh5 = 34.8
ds5 = 48.5
dg5 = -1.8
dh5t = 70.1
ds5t = 48.5
dg5t = 33.5
dh6 = 9
ds6 = 8.3
dg6 = 1.7
T = 450 + 273
R = 1.987  / 1000
k3f, k3b, k3eq = kst(dh3, ds3, dg3, dh3t, ds3t, dg3t, T, R)  
k5f, k5b, k5eq = kst(dh5, ds5, dg5, dh5t, ds5t, dg5t, T, R)  
    
