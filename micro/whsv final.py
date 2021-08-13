# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 21:38:28 2021

@author: tjczec01@gmail.com


"""

P = 101.325*1000  # pa
PE = 5.25 * 1000
etohmass = 46.069  # g·mol−1
etohden = 0.78945  # g/cm3 20 °C g/ml 20 °C
whsvl = 21  # h^-1
gcat = 0.0067  # g_cat

Ne = (whsvl*gcat)/(60*etohmass)  #  molE/min
Ye = PE/P
Nh = Ne/Ye - Ne  #  molHe/min
Vhe = Nh * 22.4 * 1000  # ml/min
Ge = Ne*etohmass  #  gE/min
Ve = Ge/etohden  #  mLE/min

print("Ethonal Molar flow rate: {:.3} mol/min".format(Ne))
print("Helium Molar flow rate: {:.3} mol/min".format(Nh))
print("Ethonal Molar fraction: {:.3}%".format(Ye * 100))
print("Helium Volumetric flow rate: {:.3} mL/min".format(Vhe))
print("Ethonal Volumetric flow rate: {:.3} mL/min".format(Ve))
print("Ethonal Volumetric flow rate: {:.3} μL/min".format(Ve * 1000))