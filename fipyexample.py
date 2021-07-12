# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 16:15:24 2021

@author: travis.czechorski

Email:   travis.czechorski@maine.edu
         tjczec01@gmail.com
        
Advisor: thomas.schwartz@maine.edu
       
Github:  https://github.com/tjczec01       

"""

from fipy import CellVariable, Grid1D, TransientTerm, VanLeerConvectionTerm, DiffusionTerm, ImplicitSourceTerm, ConvectionTerm, CentralDifferenceConvectionTerm, Viewer
from fipy.tools import numerix

molarWeight = 0.118
ee = -0.455971
gasConstant = 8.314
temperature = 650.
vbar = 1.3e-05
liquidDensity = 7354.3402662299995
vaporDensity = 82.855803327810008

def f(rho):
    return ee * rho**2 / molarWeight**2 + gasConstant * temperature * rho / molarWeight * \
           numerix.log(rho / (molarWeight - vbar * rho))
           
def mu(rho):
    return 2 * ee * rho / molarWeight**2 + gasConstant * temperature / molarWeight * \
           (numerix.log(rho / (molarWeight - vbar * rho)) + molarWeight / (molarWeight - vbar * rho))
           
def P(rho):
    return rho * mu(rho) - f(rho)

print(numerix.allclose(mu(liquidDensity), mu(vaporDensity)))

print(numerix.allclose(P(liquidDensity), P(vaporDensity)))

Lx = 1e-6
nx = 100
dx = Lx / nx
mesh = Grid1D(nx=nx, dx=dx)

density = CellVariable(mesh=mesh, hasOld=True, name=r'$\rho$')
velocity = CellVariable(mesh=mesh, hasOld=True, name=r'$u$')
densityPrevious = density.copy()
velocityPrevious = velocity.copy()

potentialNC = CellVariable(mesh=mesh, name=r'$\mu^{NC}$')

epsilon = 1e-16
freeEnergy = (f(density) + epsilon * temperature / 2 * density.grad.mag**2).cellVolumeAverage

matrixDiagonal = CellVariable(mesh=mesh, name=r'$a_f$', value=1e+20, hasOld=True)
correctionCoeff = mesh._faceAreas * mesh._cellDistances / matrixDiagonal.faceValue
massEqn = TransientTerm(var=density) \
          + VanLeerConvectionTerm(coeff=velocity.faceValue + correctionCoeff \
                                        * (density * potentialNC.grad).faceValue, \
                                  var=density) \
          - DiffusionTerm(coeff=correctionCoeff * density.faceValue**2, var=potentialNC)
          
viscosity = 1e-3
ConvectionTerm = CentralDifferenceConvectionTerm
momentumEqn = TransientTerm(coeff=density, var=velocity) \
              + ConvectionTerm(coeff=[[1]] * density.faceValue * velocity.faceValue, var=velocity) \
              == DiffusionTerm(coeff=2 * viscosity, var=velocity) \
              - ConvectionTerm(coeff=density.faceValue * [[1]], var=potentialNC) \
              + ImplicitSourceTerm(coeff=density.grad[0], var=potentialNC)
              
velocity.constrain(0, mesh.exteriorFaces)
potentialDerivative = 2 * ee / molarWeight**2 + gasConstant * temperature * molarWeight / density / (molarWeight - vbar * density)**2
potential = mu(density)

potentialNCEqn = ImplicitSourceTerm(coeff=1, var=potentialNC) \
                 == potential \
                 + ImplicitSourceTerm(coeff=potentialDerivative, var=density) \
                 - potentialDerivative * density \
                 - DiffusionTerm(coeff=epsilon * temperature, var=density)
                 
potentialNC.faceGrad.constrain(value=[0], where=mesh.exteriorFaces)

coupledEqn = massEqn & momentumEqn & potentialNCEqn

numerix.random.seed(2011)
density[:] = (liquidDensity + vaporDensity) / 2 * \
   (1  + 0.01 * (2 * numerix.random.random(mesh.numberOfCells) - 1))
   
from fipy import input
if __name__ == '__main__':
    viewers = Viewer(density), Viewer(velocity), Viewer(potentialNC)
    for viewer in viewers:
        viewer.plot()
    input('Arrange viewers, then press <return> to proceed...')
    for viewer in viewers:
        viewer.plot()
        
cfl = 0.1
tolerance = 1e-1
dt = 1e-14
timestep = 0
relaxation = 0.5
if __name__ == '__main__':
    totalSteps = 1e10
else:
    totalSteps = 10
    
while timestep < totalSteps:

    sweep = 0
    dt *= 1.1
    residual = 1.
    initialResidual = None

    density.updateOld()
    velocity.updateOld()
    matrixDiagonal.updateOld()

    while residual > tolerance:

        densityPrevious[:] = density
        velocityPrevious[:] = velocity
        previousResidual = residual

        dt = min(dt, dx / max(abs(velocity)) * cfl)

        coupledEqn.cacheMatrix()
        residual = coupledEqn.sweep(dt=dt)

        if initialResidual is None:
            initialResidual = residual

        residual = residual / initialResidual

        if residual > previousResidual * 1.1 or sweep > 20:
            density[:] = density.old
            velocity[:] = velocity.old
            matrixDiagonal[:] = matrixDiagonal.old
            dt = dt / 10.
            if __name__ == '__main__':
                print('Recalculate the time step')
            timestep -= 1
            break
        else:
            matrixDiagonal[:] = coupledEqn.matrix.takeDiagonal()[mesh.numberOfCells:2 * mesh.numberOfCells]
            density[:] = relaxation * density + (1 - relaxation) * densityPrevious
            velocity[:] = relaxation * velocity + (1 - relaxation) * velocityPrevious

        sweep += 1

    if __name__ == '__main__' and timestep % 10 == 0:
        print('timestep: %e / %e, dt: %1.5e, free energy: %1.5e' % (timestep, totalSteps, dt, freeEnergy))
        for viewer in viewers:
            viewer.plot()

    timestep += 1
    
from fipy import input
if __name__ == '__main__':
    input('finished')
print(freeEnergy < 1.5e9)