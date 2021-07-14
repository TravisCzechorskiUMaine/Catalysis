# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 09:24:48 2021

@author: travis.czechorski

Email:   travis.czechorski@maine.edu
         tjczec01@gmail.com
        
Advisor: thomas.schwartz@maine.edu
       
Github:  https://github.com/tjczec01       

"""

'''
GHSV = Hourly  volumetric feed gas flow rate/reaction volume
LHSV = Hourly  volumetric feed liquid flow rate/reaction volume
WHSV = Hourly mass feed flow rate/ Catalyst mass
'''

def xin(ni, n):
    return ni/n

def xip(pi, p):
    return pi/p

def p_i(xi, p):
    return xi*p

T = 723  # K
P = 101325  # pa
Pk = 101.325  # kpa
PE = 5.25 * 1000
mgosio2 = 100.3887 # g/mol
sa = 80.4 # m^2/g
cacm = sa*10000.0  # cm^2/g
sacm = 1.0 / cacm  # g/cm^2
saf = sacm/mgosio2  #  mol/cm^2
pi = 3.1415926535
etohmass = 46.069  # g·mol−1
etohden = 0.78945  # g/cm3 20 °C g/ml 20 °C
outerd = 1/8 # inch # https://www.mcmaster.com/tubing/system-of-measurement~inch/od~1-8/smooth-bore-seamless-stainless-steel-tubing/
odm = outerd * 0.0254
idl = [0.085, 0.069, 0.055, 0.027]  # inch
idm = [x*0.0254 for x in idl]  # m
cam = [pi*(x**2) for x in idm]  # m^2
R = 8.31446261815324	# J⋅K−1⋅mol−1
Ratm = 8.20573660809596*(10**-5)  # m3⋅atm⋅K−1⋅mol−1
whsvl = [3.7, 5.2, 21]  # h^-1
gms_m = [(x*50)/60 for x in whsvl]  # g/min
mls_m = [x/etohden for x in whsvl]  # ml/min
gms_s = [x/60 for x in gms_m]  # g/s
ms_s = [x/mgosio2 for x in gms_m]  # mol/s
n = ms_s[-1] / (PE/P)  # mols
V = (n*Ratm*T)/P  # m^3/s
petoh = [3.39, 8.71]
xetoh = [xip(x, P) for x in petoh]
V_ml = V*1E6   # cm^3/s
V_mle = V_ml * (PE/P)  # ml/s

print([round(n, 3) for n in ms_s])
print([round(n, 3) for n in mls_m])
print("Total Moles: {:.3}".format((ms_s[-1] / (PE/P))))
print("{:.3} %".format(PE/P * 100))
print("Total Flow: {:.3} ml".format(V_ml))
print("Total Ethanol Flow: {:.3} ml".format(V_mle))
print("Ethanol Mole Fraction Range: {:.3} - {:.3}".format(xetoh[0], xetoh[1]))
print("Ethanol Mole Fraction Range: {:.3} - {:.3} %".format(xetoh[0] * 100, xetoh[1] * 100))
