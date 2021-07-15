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

T = 273.15  # K
P = 101325  # pa
Pk = 101.325  # kpa
PE = 5.25 * 1000
xi = PE/P
mgosio2 = 100.3887 # g/mol
sa = 80.4 # m^2/g
cacm = sa*(100.0**2)  # cm^2/g
sacm = 1.0 / cacm  # g/cm^2
saf = sacm/mgosio2  #  mol/cm^2
pi = 3.1415926535
HE = 0.1786/1000   # g/mL
HEG = 4.002602  # g/mol  
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
gcat = 0.0067  # g_cat
gms_m = [(x*(1/gcat))/60 for x in whsvl]  # g/min
mls_m = [x/etohden for x in gms_m]  # ml/min
# gms_s = [x/60 for x in gms_m]  # g/s
# print(gms_s[-1]/etohden)
ms_m = [x/etohmass for x in gms_m]  # mol/min
nl = [x  / (PE/P) for x in ms_m]
# print(nl)
n = ms_m[-1] / xi  # mols
V = (n*Ratm*T)/P  # m^3/min
Vcm = V * 1E6
petoh = [3.39, 8.71]
xetoh = [xip(x*1000, P) for x in petoh]
V_ml = V*1E6   # cm^3/min
# V_mlm = Vcm * 60   # ml/min
V_mle = V_ml * xi  # ml/min
# V_mlem = V_mle * 60  #  mL/min
V_he = V_ml - V_mle  #  mL/min
print((V_mle/V_ml)*100,  xi * 100)
print(V_mle, V_he)
# print(50/16.936267904840758)
# print(0.05*2.95)
# print(15.424652008051691*2.95)
# print(V_he)
# print(V_mle/V_ml * 100)
# print(V_mle*1000)
# print(V_ml*60 - V_mle)
# print(V_ml*60)
# print(V*10E6 - V_mle)
# print((0.925/17.9)*100)
# etflow = 0.925/etohmass  # g/min
# etf = etflow/etohden
# print(etf*1000)
he1 = 16.925*22.4  # L
he2 = he1/1000
# print(he2)
# print([round(n, 3) for n in ms_m])
# print([round(n, 3) for n in mls_m])
print("Total Moles: {:.3}".format((ms_m[-1] / xi)))
print("{:.3} %".format(xi * 100))
print("Total Flow: {:.3} ml/min".format(V_ml))
print("Total Ethanol Flow: {:.4} microL/min".format(V_mle*1000))
print("Total Helium Flow: {:.3} mL/min".format(V_he))
print("Ethanol Mole Fraction Range: {:.3} - {:.3}".format(xetoh[0], xetoh[1]))
print("Ethanol Mole Fraction Range: {:.3} - {:.3} %".format(xetoh[0] * 100, xetoh[1] * 100))
