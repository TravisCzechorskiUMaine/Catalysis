# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 13:28:28 2021

@author: travis.czechorski

Email:   travis.czechorski@maine.edu
         tjczec01@gmail.com
        
Advisor: thomas.schwartz@maine.edu
       
Github:  https://github.com/tjczec01       

"""

from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
import numpy as np
import scipy.optimize as optimize
import matplotlib.pylab as plt

def func(kd, p0, l0):
    return 0.5 * (-1 - ((p0 + l0)/kd) + np.sqrt(4 * (l0/kd) + (((l0 - p0)/kd) - 1)**2))

def residuals(kd, p0, l0, PLP):
    return PLP - func(kd, p0, l0)

N=1000
kd_guess=3.5  # <-- You have to supply a guess for kd
p0 = np.linspace(0, 10, N)
l0 = np.linspace(0, 10, N)
PLP = func(kd_guess, p0, l0) + (np.random.random(N) - 0.5) * 0.1

kd, cov, infodict, mesg, ier = optimize.leastsq(residuals, kd_guess, args=(p0, l0, PLP), full_output=True)

print(kd)

PLP_fit=func(kd, p0, l0)

plt.plot(p0, PLP, '-b', p0, PLP_fit, '-r')
plt.show()

# create data to be fitted
x = np.linspace(0, 15, 301)
data = (5.0 * np.sin(2 * x - 0.1) * np.exp(-x * x * 0.025) + np.random.normal(size=len(x), scale=0.2))

# define objective function: returns the array to be minimized
def fcn2min(params, x, data):
    """ model decaying sine wave, subtract data"""
    amp = params['amp']
    shift = params['shift']
    omega = params['omega']
    decay = params['decay']
    model = amp * np.sin(x * omega + shift) * np.exp(-x*x*decay)
    return model - data

# create a set of Parameters
params = Parameters()
params.add('amp',   value= 10,  min=0)
params.add('decay', value= 0.1)
params.add('shift', value= 0.0, min=-np.pi/2., max=np.pi/2)
params.add('omega', value= 3.0)


# do fit, here with leastsq model
minner = Minimizer(fcn2min, params, fcn_args=(x, data))
kws  = {'options': {'maxiter':10}}
result = minner.minimize()


# calculate final result
final = data + result.residual

# write error report
report_fit(result)

# try to plot results
try:
    import pylab
    pylab.plot(x, data, 'k+')
    pylab.plot(x, final, 'r')
    pylab.show()
except:
    pass