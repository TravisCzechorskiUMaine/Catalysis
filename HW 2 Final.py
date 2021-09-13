# -*- coding: utf-8 -*-
"""
Created on Sun Sep 12 06:22:24 2021

@author: tjczec01@gmail.com


"""

import math as mt
import numpy as np
import matplotlib.pylab as plt
import os
from scipy.constants import hbar, Planck, speed_of_light, Avogadro

csc, hsc, NA, hbarsc = speed_of_light, Planck, Avogadro, hbar  # SI units
clear = os.system('cls')
cwd = os.getcwd()
dir_path = os.path.dirname(os.path.realpath(__file__))
path_fol = r"{}\HomeWork".format(dir_path)

print("\n")
print("Current working directory:\n")
print("{}\n".format(cwd))

try:
    os.mkdir(path_fol)
    print("New Folder was created\n")
    print("Current working directory - Created Folder Path:\n")
    print("{}\n".format(path_fol))
except Exception:
    print("Current working directory - Current Folder Path:\n")
    print("{}\n".format(path_fol))
    pass

font = {'family': 'serif',
        'color': 'black',
        'weight': 'normal',
        'size': 12,
        'ha': 'center',
        'va': 'bottom'}

def indices(lst, element):
    result = []
    offset = -1
    while True:
        try:
            offset = lst.index(element, offset+1)
        except ValueError:
            return result
        result.append(offset)

def k1(a, De):
    t1 = 2.0 * (a**2)
    t2 = De * t1
    return t2

def k2(v, u):
    pi = 3.141592653589793
    t1 = 4.0 * (pi**2)
    t2 = u * (v**2)
    return t1 * t2

def uf(m1, m2):
    M1 = 1.0 / m1
    M2 = 1.0 / m2
    MSUM = M1 + M2
    MF = 1.0 / MSUM
    return MF


def ur(m1, m2):
    M1 = m1/NA
    M2 = m2/NA
    MSUM = M1 + M2
    MMul = M1*M2
    MF = MMul / MSUM
    return MF   # kg/molecule

def Vx(k, r, re):
    vv0 = r - re
    vx1 = vv0**2
    vx = 0.5*k*vx1
    return vx

def A(ve, De, u):
    pi = 3.141592653589793
    a1 = 2.0 * pi * ve 
    a2 = 2.0 * De
    a3 = mt.sqrt(u/a2)
    af = a1 * a3
    # aff = 2*np.pi*ve * np.sqrt(μ/(2*De)) 
    return af  # m^-1

def V(r, re, a, De):
    term1 = r - re
    term2 = term1 * a
    term3 = 1.0 / (np.exp(term2))
    term4 = 1.0 - term3
    term5 = term4**2
    term6 = term5 * De
    return term6

def Vd(delrv, a, De):
    term1 = delrv
    term2 = term1 * a
    term3 = np.exp(-1 * term2)
    term4 = 1.0 - term3
    term5 = term4**2
    term6 = term5 * De
    return term6

def VE(r, re, a, Dev):
    '''The Morse energy for harmonic oscillator at distance r'''
    return Dev * (1 - np.exp(-a*(r-re)))**2

def Vm(r, re, k, D):
    "Morse potential at r given spectroscopic constants"

    alpha = np.sqrt(k/(2*D))
    return D * (1 - np.exp(-alpha*(r - re)))**2

def Evib(n):
    '''Given an (integer) quantum number v, return the vibrational energy in joules.'''
    return hsc*csc * ( (n+0.5)*we - (n+0.5)**2 * wexe )

def En(n, De, h, v):
    
    Xe1 = h * v
    Xe2 = 4.0 * De
    Xe = Xe1/Xe2
    term1 = h * v
    term2 = Xe * term1
    n12 = 1/2
    nterm = n + n12
    Ee = nterm * term1 - (nterm**2) * term2
        
    return Ee

def Eho(ve, n=0):
    """Harmonic oscillator energy (J) of vibrational level v.
    nu is frequency in s^-1"""

    return h*ve*(n+1/2)

def Eaho(ve, n, D0):
    """Anharmonic energy (J) of vibrational level v.
    nu=frequency (s^-1)
    D0=depth (J/molec)"""

    De = D0 + Eho(ve, 0)
    vv = h*ve*(n+1/2) - (h*ve)**2 / (4*De) * (n+1/2)**2 
    return vv  # *(10**9)  # h*nu*(v+1/2) - (h*nu)**2 / (4*De) * (v+1/2)**2

def calc_turning_pts(a, E, re, DE):
    """Calculate the classical turning points at energy E.

    Returns rm and rp, the classical turning points of the Morse
    oscillator at energy E (provided in J). rm < rp.

    """

    b = np.sqrt(E / DE)
    rm = re - np.log(1+b) / a
    rp = re - np.log(1-b) / a
    return rm, rp

# 1 nm = 10 Å
# 1 J = kg⋅m2⋅s−2
# 8065.73 cm-1 = 96.487 kJ/mol

w = np.sqrt(481/1.57*10E27)
amutokg = 1.66053892173E-27
we = 2.990e3 * 1e2  # cm^-1 -> m^-1
wexe = 52.8186 * 1e2 # cm^-1 -> m^-1
csc, hsc, NA = speed_of_light, Planck, Avogadro
FAC = 100 * hsc * csc # Factor for conversion from cm-1 to J
KFAC = FAC/1000
cmtokj = 96.487 / 8065.73 #  1 cm**-1 = 0.011962587391345855 kJ/mol
kjtocm = 8065.73 / 96.487  #  1 kJ/mol =  83.59395566242084 cm**-1
h = 6.62607004 * (10**-34) # E10-34  # m2 kg / s
Na = 6.0221409 * (10**23)  #  molecules/mol
Hm = 1.00782503207E-3  #  amu
CLm = 0.5*34.96885268E-3  #  amu
De = 7.427 * (10**-19) #  J
DekJ = De / 1000 #  kJ
JtoCM = 5.03445 * (10**22)#E22  # cm**-1
DeCm = De*JtoCM  # cm**-1
ve = 9.004E13 # * (10**13)#E13 #  Hz  or # 1/s
# re = 1.27455e-8 * 1E-2  # cm
re = 127.4E-12 #* (10**-12)  # m

zpe = Eho(ve, 0) 
DF = De + zpe
u = uf(Hm, CLm)     #  amu
ukg = ur(Hm, CLm)  # kg
AA = A(ve, De, u/NA)     #  1/m
AAkg = A(ve, De, ukg) #   * NA
xstart = -1.5E-10 #  m
xend = 6.5E-10 #  m
xstep = 1E-13 #  m
delr = [i for i in np.arange(xstart, xend, xstep)] #  m
delrnm = [i * (10**9) for i in np.arange(xstart, xend, xstep)] #  nm
delra = [i * 10 for i in delrnm] #  nm
VRl = [] #  kJ/mol
ens = [i for i in range(0, 16, 1)]
ensb = [i for i in range(0, 16, 1)]
Elistl = [En(i, De, h, ve) for i in ens]
Elistl0 = [Eho(ve, i) for i in ensb]
kk1 = k1(1.88705191439686622620E10, 7.427E-19) 
kk2 = k2(ve, ur(Hm, CLm))
EE = []
EEF = []
EL0 = [i/FAC for i in Elistl0]

for ii in ens:
    Elist = En(ii, De, h, ve)/FAC
    EE.append(Elist)
    EEF.append(Eho(ve, ii)/(FAC))
for i in delr:
    Vv = Vd(i, 1.88705191439686622620E10, De) 
    VRl.append(Vv*NA*1E-3*83.59395566242084)
        
VRV = np.argmin(np.abs(np.array(VRl) - 45000.0)) 
XRV = np.argmin(np.abs(np.array(VRl) - 0.1))
Vxl = [] #  kJ/mol

for i in delrnm:
    Kk = 2.0*(De)*(AAkg**2)*NA
    KK = (528.9457303499935 + 520.6252625322837)/2
    Vxv = Vx(Kk*1E-19, i, re) 
    Vxl.append(Vxv)

VRlcm = [i * 0.011962587391345855 for i in VRl]
Vxlcm = [i * 0.01196258739134585 for i in Vxl]
VRlk = [i/1000 for i in VRl]
VRlkcm = [i/10 for i in VRlcm]
Vxlkcm = [i/(100 * 0.01196258739134585) for i in Vxlcm]

DeGraph = DeCm #* 0.011962587391345855
xgstart = -0.15/10
xgend1 = 0.65/10
xgstep1 = 0.05/10

xgstart2 = -1.5/10
xgend2 = 3.5/10
xgstep2 = 0.5/10

xgstart3 = -0.15/10
xgend3 = 0.1/10
xgstep3 = 0.05/10

figa = plt.figure()
plt.plot(delrnm, VRl, color='tab:blue', label=r"$\bf{V(r) = D_{e} * (1 - e^{-a*(r - r_{r})})^2}$")
# plt.plot(delrnm, Vxl, color='tab:green', label=r"$\bf{V(r) = k*(r - r_{r})^2}$")
plt.xlabel(r'$\bf{\Delta{R}} \ = \ (r - r_{r})$' + r'$\ _{\it{[nm]}}$')
plt.ylabel(r'$\bf{V(r)}$' +  r'$\ [\frac{1}{cm}$]')
plt.legend(loc='best', fontsize='large')
plt.title(label=r"$\bf{Morse} \ \bf{potential} \ (Full)$", fontdict=font)
plt.xlim([delrnm[0], delrnm[-1]])
plt.grid()
figa.savefig("{}\FullGraphHw2.pdf".format(path_fol), bbox_inches='tight')
figa.savefig("{}\FullGraphHw2.svg".format(path_fol), bbox_inches='tight')
plt.show()
plt.close()

figb = plt.figure()
plt.plot(delrnm, VRl, color='tab:blue', label=r"$\bf{V(r) = D_{e} * (1 - e^{-a*(r - r_{r})})^2}$")
plt.plot(delrnm, Vxl, color='tab:green', label=r"$\bf{V(r) = \frac{1}{2}*k*(r - r_{r})^2}$")
plt.axhline(y=float(DeCm), color='k', linestyle='--', label=r"$De$")
plt.xlabel(r'$\bf{\Delta{R}} \ = \ (r - r_{r})$' + r'$\ _{\it{[nm]}}$')
plt.ylabel(r'$\bf{V(r)}$' +  r'$\ [\frac{1}{cm}$]')
plt.legend(loc='best')
plt.title(label=r"$\bf{Morse} \ \bf{potential}$", fontdict=font)
plt.xticks(np.arange(xgstart2, xgend2, step=xgstep2))
plt.xlim([-0.065, xgend2])
plt.ylim(0, 40000) # VRl[VRV])
plt.grid()
figb.savefig("{}\FitGraphHw2.pdf".format(path_fol), bbox_inches='tight')
figb.savefig("{}\FitGraphHw2.svg".format(path_fol), bbox_inches='tight')
plt.show()
plt.close()

figc = plt.figure()
plt.plot(delrnm[VRV:-1], VRl[VRV:-1], color='tab:blue', label=r"$\bf{V(r) = D_{e} * (1 - e^{-a*(r - r_{r})})^2}$")
plt.plot(delrnm[VRV:-1], Vxl[VRV:-1], color='tab:green', label=r"$\bf{V(r) = \frac{1}{2}*k*(r - r_{r})^2}$")
plt.xlabel(r'$\bf{\Delta{R}} \ = \ (r - r_{r})$' + r'$\ _{\it{[nm]}}$')
plt.ylabel(r'$\bf{V(r)}$' +  r'$\ [\frac{1}{cm}$]')
plt.legend(loc='best')
plt.title(label=r"$\bf{Morse} \ \bf{potential}$", fontdict=font)
plt.xticks(np.arange(xgstart3, xgend3, step=0.0005))
plt.xlim([-0.0015, 0.0015])
plt.ylim(0, 20)
plt.grid()
figc.savefig("{}\ZoomedGraphHw2.pdf".format(path_fol), bbox_inches='tight')
figc.savefig("{}\ZoomedGraphHw2.svg".format(path_fol), bbox_inches='tight')
plt.show()
plt.close()

figd = plt.figure()
plt.plot(delrnm, VRl, color='tab:blue', label=r"$\bf{V(r) = D_{e} * (1 - e^{-a*(r - r_{r})})^2}$")
plt.plot(delrnm, Vxl, color='tab:green', label=r"$\bf{V(r) = \frac{1}{2}*k*(r - r_{r})^2}$")

for ii in ens:
    Elist = En(ii, De, h, ve)
    rm, rp = calc_turning_pts(AAkg, Elistl[ii], re/100, De)
    rm2, rp2 = calc_turning_pts(AAkg, Elistl0[ii], re/100, De)
    plt.hlines(Elistl[ii]/FAC, rm*1E9, rp*1E9, colors='b', linestyles='-', lw=1)
    # plt.hlines(Elistl0[ii]/FAC, 1.2*rm*1E9, 0.325*rp*1E9, colors='g', linestyles='-', lw=2)
    
plt.xlabel(r'$\bf{\Delta{R}} \ = \ (r - r_{r})$' + r'$\ _{\it{[nm]}}$')
plt.ylabel(r'$\bf{V(r)}$' +  r'$\ [\frac{1}{cm}$]')
plt.legend(loc='best')
plt.title(label=r"$\bf{Morse} \ \bf{potential}$", fontdict=font)
plt.xticks(np.arange(xgstart2, xgend2, step=xgstep2))
plt.xlim([-0.065, xgend2])
plt.ylim(0, 40000) # VRl[VRV])
plt.grid()
figd.savefig("{}\FitGraphHw2.pdf".format(path_fol), bbox_inches='tight')
figd.savefig("{}\FitGraphHw2.svg".format(path_fol), bbox_inches='tight')
plt.show()
plt.close()
