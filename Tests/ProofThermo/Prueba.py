import numpy as np
import matplotlib.pylab as plt
from numpy import linalg as LA
import sys
import random
import cmath
import math
np.random.seed(int(sys.argv[1]))
hbar=1.0
# epsilon=0.1
epsilon=1.0
nu=8.0
tau=3.0
gamma=5*10**(-3)
# g=0.625
g=6.25
dt=0.01
tmax=tau
N=int(tmax/dt)
sigma=math.sqrt(dt)
mu=0.0
trajectories=1
#----------------------------Pauli matrices ------------------------------
def sigmax():
    return np.matrix([[.0,1.],[1.,.0]])
def sigmay():
    return np.matrix([[.0,.0+1.*1j],[.0-1.*1j,.0]])
def sigmaz():
    return np.matrix([[-1.0,0.0],[0.0,1.0]])
def sigmaplus():
    return np.matrix([[.0,.0],[1.,0.]])
def sigmaminous():
    return np.matrix([[.0,1.],[.0,.0]])
#--------------------------Driving-------------------------------
def Lambda(t):
    return g*(1./(np.cosh(nu*(1.-(t/tau)))))
#---------------------- Hamiltonian---------------------------------
def H(t):
    return (epsilon)*sigmaz()+Lambda(t)*sigmax()
#---------------------- Lindblad aperators ----------------------------
def J():
    return math.sqrt(gamma)*sigmaz()
def Jdag():
    return J().conjugate().transpose()
#-----------------------state function-------------------------------
def state(a,b):
    return ket(a,b)*bra(a,b)
#--------------------- State of the system-------------------------------
def ket(a,b):
    return np.matrix([[a],[b]])
def bra(a,b):
    return np.matrix([[a.conjugate(),b.conjugate()]])
#----------------------method-------------------------------------------
def euler(a,b,t,dt,dW):
    A1=(np.identity(2)-Jdag()*J()*0.5-(1.0*1j/hbar)*H(t))*ket(a,b)
    B1=J()*ket(a,b)
    method=A1*dt+B1*dW
    return method
def Normalise(a,b):
    A=ket(a,b)
    norm=cmath.sqrt(A.conjugate().transpose()*A)
    A=A/norm
    return A

rho1 = state(1.,0.)
dXi = -5.4974*10**-2.0
nuevo = euler(1.0,0.0,0,dt,dXi)
nuevo = Normalise(nuevo[0,0],nuevo[1,0])
rho2 = state(nuevo[0,0],nuevo[1,0])
dH=H(dt)-H(0)
drho=rho2-rho1
#print (rho2*dH).trace()[0,0]
pruebas=np.matrix([[1,2+3*1j],[2-3*1j,2]])
print pruebas
W, V = LA.eig(pruebas)
print W
print "----------------------------------------------"
print V[:,0]
print "----------------------------------------------"
print V[:,1]
