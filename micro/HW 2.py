# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 13:17:30 2021

@author: travis.czechorski

Email:   travis.czechorski@maine.edu
         tjczec01@gmail.com
        
Advisor: thomas.schwartz@maine.edu
       
Github:  https://github.com/tjczec01       

"""

import math as mt
import numpy as np
import matplotlib.pylab as plt

# print(mt.pi)

def ur(m1, m2):
    M1 = 1/m1
    M2 = 1/m2
    MSUM = M1 + M2
    MF = 1/MSUM
    return MF


def A(ve, De, μ):
    print("Ve: {}, De: {}, μ: {}".format(ve, De, μ))
    a1 = 2*mt.pi*ve
    print(a1)
    a2 = (2*De)/μ
    print(a2)
    a3 = a2**(-1/2)
    print(a3)
    af = a1*a3
    print(af)
    return af

def V(r, re, a, De):
    term1 = r - re
    term2 = term1*a
    term3 = 1/(mt.exp(term2))
    term4 = 1.0 - term3
    term5 = De*(term4**2)
    return term5

def Vd(delrv, a, De):
    term1 = delrv
    term2 = term1*a
    term3 = mt.exp(-term2)
    term4 = 1.0 - term3
    term5 = De*(term4**2)
    return term5

rlist = [0]
delr = [i for i in np.arange(-1.0E-10, 6.2E-10, 1E-11)]
print(delr)
for i in range(-1, 34, 1):
    rlistv = rlist[-1]
    rlistv += 1E-12
    rlist.append(rlistv)
# print(rlist)
VRl = []
for i in delr:
    Vv = Vd(i, 18721292325.403297, 37390.86015) 
    VRl.append(Vv)
print(VRl)
Hm = 1.0078  #  amu
CLm = 34.9688  #  amu
De = 7.427E-19 #  J
JtoCM = 5.03445E22  # cm**-1
DeCm = De*JtoCM  # cm**-1
# print(DeCm)
ve = 9.004E13 #  Hz
re = 127.4E-12  # m
UR = ur(Hm, CLm)  #  amu
URkg = UR*1.6605402E-27  # kg
print("URkg: {}".format(URkg))
AA = A(ve, De, URkg)  #  1/m
AAkg = AA






fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(delr, VRl, color='tab:blue')
# ax.plot(rlist, ydata2, color='tab:orange')

print(ur(Hm, CLm))
print("A: {:.3}".format(AA))