# -*- coding: utf-8 -*-
"""
Created on Mon Sep  6 15:00:02 2021

@author: tjczec01@gmail.com


"""

from sympy import *
import pprint as pp
init_printing(use_unicode=True)

x, y, a, s = symbols('x y a s')

for n in range(0, 12, 1):
    integ = Integral((x**n)*exp(-a*x**2), (x, 0, oo))
    derl = [sqrt(pi)/(2*sqrt(a)), (1/(2*a)), (sqrt(pi)/(4*a**(3/2))), 1/(2*a**2), 3*sqrt(pi)/(8*a**(5/2)), 15*sqrt(pi)/(16*a**(7/2)), 3/a**4, 105*sqrt(pi)/(32*a**(9/2)), 12/a**5, 945*sqrt(pi)/(64*a**(11/2))]
    
    print(n)
    print(integ)
    print(integ.doit())
    print(diff(derl[n], a))