# -*- coding: utf-8 -*-
"""

Created on Mon Sep 13 19:15:21 2021

@author: Travis J Czechorski

Github: https://github.com/tjczec01

Email:   travis.czechorski@maine.edu
         tjczec01@gmail.com
        
Advisor: thomas.schwartz@maine.edu
       
Github:  https://github.com/tjczec01 
         https://github.com/TravisCzechorskiUMaine

"""

from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
import numpy as np
import scipy.optimize as optimize
import matplotlib.pylab as plt
import math
import pprint as pp

# for i in range(1, 7, 1):
#     print("dh{} = \nds{} = \ndg{} = ".format(i, i, i))

def ks(dh, ds, dg, dht, dst, dgt, T, R):
    kB = 1.380662 * 10E23   #  m2 kg/s2 K1
    h = 6.626176 * 10E34  #  m2 kg/s
    R = 1.987  #   cal/K mol
    Rk = 1.987  / 1000    #   kcal/K mol
    k_constant = (kB*T)/h
    Keq = np.exp(ds/R) * np.exp((-1*dh)/(R * T))
    Kf = k_constant*np.exp(dst/R) * np.exp((-1*dht)/(R * T))
    Kb = Kf/Keq
    

def func(P, T, R, dh1, ds1, dg1, dh2, ds2, dg2, dh3, ds3, dg3, dh3t, ds3t, dg3t, dh4, ds4, dg4, dh5, ds5, dg5, dg5t, dh5t, ds5t, dh6, ds6, dg6):
    kB = 1.380662 * 10E23   #  m2 kg/s2 K1
    h = 6.626176 * 10E34  #  m2 kg/s
    Rcal = 1.987  #   cal/K mol
    Rk = 1.987  / 1000 
    k1 = math.exp(ds1/R) * math.exp(dh1/(R*T))
    k2 = math.exp(ds2/R) * math.exp(dh2/(R*T))
    k3 = math.exp(ds3/R) * math.exp(dh3/(R*T))
    k4 = math.exp(ds4/R) * math.exp(dh4/(R*T))
    k5 = math.exp(ds5/R) * math.exp(dh5/(R*T))
    k6 = math.exp(ds6/R) * math.exp(dh6/(R*T))
    # Retoh = 



def ract(Pa, T, dh1, ds1, dg1, dh2, ds2, dg2, dh3t, ds3t, dg3t, totp):
    # T = 723.15
    # To = 298.15
    P = Pa * 1000
    Pp = Pa/totp
    kb = 1.380662E-23   #  m2 kg/s2 K1
    h = 6.626176E-34  #  m2 kg/s
    R = 1.987  #   cal/K mol
    Rk = 1.987  / 1000    #   kcal/K mol
    rkt = Rk* T
    rt = R * T
    k1 = math.exp(ds1/R) * math.exp(-dh1/rkt)
    k1g = math.exp(-dg1/rt)
    # print(k1, k1g)
    k2 = math.exp(ds2/R) * math.exp(-dh2/rkt)
    k2g = math.exp(-dg2/rt)
    # print("K1: {}, K1G: {} | Difference {:.3} %".format(k1, k1g, abs((1 - k1/k1g) * 100)))
    # print("K2: {}, K2G: {} | Difference {:.3} %".format(k2, k2g, abs((1 - k2/k2g) * 100)))
    term1 = (kb*T)/h
    # print("{:.3}".format(term1))
    k3 = term1*(math.exp(ds3t/R) * math.exp(-dh3t/(Rk*T)))
    
    term2 = (((math.exp(ds3t/R) * math.exp((-1*dh3t)/rkt))) * math.exp(ds1/R) * math.exp((-1*dh1)/rkt) * math.exp(ds2/R) * math.exp((-1*dh2)/rkt)) * P
    term3 = term1*term2
    term4 = (1.0 + ((math.exp(ds1/R) * math.exp(-1*dh1/rkt) * math.exp(ds2/R) * math.exp(-1*dh2/rkt)) * P))
    rate = (k3*k1*k2*P)/(1.0 + (k1*k2*P))
    rateg = (k3*k1g*k2g*P)/(1.0 + (k1g*k2g*P))
    rate2 =  term3/(1.0 + (k1*k2*P))
    # kb*t/ h [1/s]
    # print(rate, rateg, rate2)
    # print("Rate 1: {}, Rate 1G: {}, Rate 2: {} | Difference 1-1G: {:.3}%, 1-2: {:.3}%, 1G-2 {:.3} %".format(rate, rateg, rate2, abs((1 - rate/rateg) * 100), abs((1 - rate/rate2) * 100), abs((1 - rateg/rate2) * 100)))

    return rate*2*5.09E-5  #((((kb*T)/h)*(math.exp(ds3t/R) * math.exp((-1*dh3t)/(R*T)))) * math.exp(ds1/R) * math.exp((-1*dh1)/(R*T)) * math.exp(ds2/R) * math.exp((-1*dh2)/(R*T)) * P) / (1 + math.exp(ds1/R) * math.exp((-1*dh1)/(R*T)) * math.exp(ds2/R) * math.exp((-1*dh2)/(R*T)) * P)
    

def rety(Pe, Tb, dh1, ds1, dg1, dh2, ds2, dg2, dh5t, ds5t, dg5t, totp):
    # T = 723
    # To = 298.15
    Pb = Pe * 1000.0
    Ppb = Pe/totp
    kbb = 1.380662E-23   #  m2 kg/s2 K1
    hb = 6.626176E-34  #  m2 kg/s
    R = 1.987   #  cal/K mol
    Rk = 1.987 / 1000    #   kcal/K mol
    rktb = Rk* Tb
    rtb = R * Tb

    k1b = math.exp(ds1/R) * math.exp(-dh1/rktb)
    k1gb = math.exp(-dg1/rktb)

    k2b = math.exp(ds2/R) * math.exp(-dh2/rktb)
    k2gb = math.exp(-dg2/rktb)
    # k3 = math.exp(ds5t/R) * math.exp((-1*dh5t)/(R*T))
    # print("K1: {}, K1G: {} | Difference {:.3} %".format(k1b, k1gb, abs((1 - k1b/k1gb) * 100)))
    # print("K2: {}, K2G: {} | Difference {:.3} %".format(k2b, k2gb, abs((1 - k2b/k2gb) * 100)))
    term1b = (kbb*Tb)/hb
    k5 = math.exp(ds5t/R) * math.exp(-dh5t/rktb) # term1*(math.exp(ds5t/Rk) * math.exp(-1*dh5t/rkt))
    rateb  = (term1b*k5*k1b*Pb)/(1.0 + k1b*k2b*Pb)
    rategb  = (term1b*k5*k1gb*Pb)/(1.0 + (k1gb*k2gb*Pb))
    term2b = term1b*k1b*Pb*k5
    term3b = 1.0 + k1b*k2b*Pb
    rate2b = term2b/term3b  #((kb*T)/h)*(((math.exp(ds5t/R) * math.exp((-1*dh5t)/rkt)) * math.exp(ds1/R) * math.exp((-1*dh1)/rkt) * P) / (1 + math.exp(ds1/R) * math.exp((-1*dh1)/rkt) * math.exp(ds2/R) * math.exp((-1*dh2)/rkt) * P))
    # print(rateb, rategb,  rate2b)
    # print("Rate 1: {}, Rate 1G: {}, Rate 2: {} | Difference 1-1G: {:.3}%, 1-2: {:.3}%, 1G-2 {:.3} %".format(rateb, rategb, rate2b, abs((1 - rateb/rategb) * 100), abs((1 - rateb/rate2b) * 100), abs((1 - rategb/rate2b) * 100)))

    return rateb*5.09E-5

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
TT = 450 + 273
norm = ract(2.565626381383230, TT, dh1, ds1, dg1, dh2, ds2, dg2, dh3t, ds3t, dg3t, 96.5266)
norme = rety(2.565626381383230, TT, dh1, ds1, dg1, dh2, ds2, dg2, dh5t, ds5t, dg5t, 96.5266)
print(norm, norme)
expir = [1.922121143897420, 1.396531653121820, 0.842796383145219, 1.120467480639130, 0.563504289209456, 0.981833673149154, 0.423245926264219, 2.052612309844440, 2.161080277016760, 2.280107727569750, 2.404224886315110, 2.565626381383230]
expir.sort()
# print(expir)
etohp = [1.396531653, 0.842796383, 1.120467481, 0.563504289, 0.981833673, 0.423245926, 1.922121144, 2.05261231, 2.161080277, 2.280107728, 2.404224886, 2.565626381, 1.922121144, 2.05261231, 2.161080277, 2.280107728, 2.404224886, 2.565626381, 1.922121144, 2.05261231, 2.161080277, 2.280107728, 2.404224886, 2.565626381, 1.922121144, 2.05261231, 2.161080277, 2.280107728, 2.404224886, 2.565626381]
etohp2 = [1.922121144, 2.05261231, 2.161080277, 2.280107728, 2.404224886, 2.565626381]
totp = [96.5266, 96.5266, 96.5266, 96.5266, 96.5266, 96.5266, 96.5266, 96.5266, 96.5266, 96.5266, 96.5266, 96.5266, 93.07922, 93.07922, 93.07922, 93.07922, 93.07922, 93.07922, 89.6318, 89.6318, 89.6318, 89.6318, 89.6318, 89.6318, 89.6318, 89.6318, 89.6318, 89.6318, 89.6318, 89.6318]
totp2 = [96.5266, 96.5266, 96.5266, 96.5266, 96.5266, 96.5266]
# R_ety = [rety(x, TT, dh1, ds1, dg1, dh2, ds2, dg2, dh5t, ds5t, dg5t, y)/norme for x, y in zip(expir, totp)]
# R_act = [ract(x, TT, dh1, ds1, dg1, dh2, ds2, dg2, dh3t, ds3t, dg3t, y)/norme  for x, y in zip(expir, totp)]
R_act = [ract(x, TT, dh1, ds1, dg1, dh2, ds2, dg2, dh3t, ds3t, dg3t, y) for x, y in zip(etohp2, totp2)]
R_ety = [rety(x, TT, dh1, ds1, dg1, dh2, ds2, dg2, dh5t, ds5t, dg5t, y) for x, y in zip(etohp2, totp2)]
tofa = (0.0067/mmetal)*(MgSi/sf)
print(R_act[0]*tofa*1E-5)
actf = [x for x in R_act]
actf.sort()
acte = [x -  1 for x in R_ety]
acte.sort()
# pp.pprint(actf)
# pp.pprint(acte)
R_act.sort()
for x in range(0, len(R_act), 1):
    print(R_act[x], R_ety[x])
# pp.pprint(R_act)
R_ety.sort()
# pp.pprint(R_ety)
# def func(kd, p0, l0):
#     return 0.5 * (-1 - ((p0 + l0)/kd) + np.sqrt(4 * (l0/kd) + (((l0 - p0)/kd) - 1)**2))

# def residuals(kd, p0, l0, PLP):
#     return PLP - func(kd, p0, l0)

# N=1000
# kd_guess=3.5  # <-- You have to supply a guess for kd
# p0 = np.linspace(0, 10, N)
# l0 = np.linspace(0, 10, N)
# PLP = func(kd_guess, p0, l0) + (np.random.random(N) - 0.5) * 0.1

# kd, cov, infodict, mesg, ier = optimize.leastsq(residuals, kd_guess, args=(p0, l0, PLP), full_output=True)

# print(kd)

# PLP_fit=func(kd, p0, l0)

# plt.plot(p0, PLP, '-b', p0, PLP_fit, '-r')
# plt.show()

# # create data to be fitted
# x = np.linspace(0, 15, 301)
# data = (5.0 * np.sin(2 * x - 0.1) * np.exp(-x * x * 0.025) + np.random.normal(size=len(x), scale=0.2))

# # define objective function: returns the array to be minimized
# def fcn2min(params, x, data):
#     """ model decaying sine wave, subtract data"""
#     amp = params['amp']
#     shift = params['shift']
#     omega = params['omega']
#     decay = params['decay']
#     model = amp * np.sin(x * omega + shift) * np.exp(-x*x*decay)
#     return model - data

# # create a set of Parameters
# params = Parameters()
# params.add('amp',   value= 10,  min=0)
# params.add('decay', value= 0.1)
# params.add('shift', value= 0.0, min=-np.pi/2., max=np.pi/2)
# params.add('omega', value= 3.0)


# # do fit, here with leastsq model
# minner = Minimizer(fcn2min, params, fcn_args=(x, data))
# kws  = {'options': {'maxiter':10}}
# result = minner.minimize()


# # calculate final result
# final = data + result.residual

# # write error report
# report_fit(result)

# # try to plot results
# try:
#     import pylab
#     pylab.plot(x, data, 'k+')
#     pylab.plot(x, final, 'r')
#     pylab.show()
# except:
#     pass

# fig, ax = plt.subplots()
# ax.scatter(expir, R_act)
