# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 12:47:46 2021

@author: travis.czechorski@maine.edu

Email:   travis.czechorski@maine.edu
         tjczec01@gmail.com
        
Advisor: thomas.schwartz@maine.edu
       
Github:  https://github.com/tjczec01       

"""

def whsv(w, p, pe, gcat, mass, den):
    Ne = (w * gcat) / (60.0 * mass)  #  molE/min
    Ye = pe/p
    Nt = Ne/Ye
    Nh = Nt - Ne  #  molHe/min
    Vhe = Nh * 22.4 * 1000  # ml/min
    Vt = Vhe/(1.0 - Ye)  # ml/min
    VHe = Vt * (1.0 - Ye)
    VE = ((((Vt - Vhe)*Ye/22.4))/den)  # ml/min
    GE = VE * den  #  g/min
    NE = Nt * Ye  #  mol/min
    Ge = NE * mass  #  gE/min
    Ve = Ge / den   #  mLE/min
    Ce = Ne/Vt  #  mol/mL
    return [Ce, Ce * 1000.0]


def ftoc(vf, den, mass, vh, p=101325.0):
    Nh = vh / (22.4 * 1000)  #  mol/min
    Ne = ((vf/1000) * den) / mass  #  mol/min
    Nt = Ne + Nh
    Ye = Ne/Nt
    PE = Ye * p
    Nt = Ne/Ye  #  total mols
    Nh = Nt - Ne   #  helium mols
    Vt = vh / (1.0 - Ye)  #  ml/min
    Ce = Ne/Vt  #  mol/mL
    return [Ce, Ce * 1000.0, PE]

Ceml = []
Ce = []
PES = []

VF = [1.0, 2.0, 2.97, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]

for i in range(0, len(VF), 1):
    # ce, cel, pes = ftoc(VF[i], 0.78945, 46.069, 11.0)
    Ceml.append(ftoc(VF[i], 0.78945, 46.069, 20.9)[0])
    Ce.append(ftoc(VF[i], 0.78945, 46.069, 20.9)[1])
    PES.append(ftoc(VF[i], 0.78945, 46.069, 20.9)[2])
    
for i, j, k in zip(Ceml, Ce, PES):
    print("Ethonal Concentration: {} mol/mL, {} mol/L, Partial Pressure: {} Pa".format(i, j, k))
