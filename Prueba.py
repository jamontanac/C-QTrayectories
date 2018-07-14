import numpy as np
import matplotlib.pylab as plt
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
    A1=((np.identity(2)-Jdag()*J()))
    return A1
    # /2.0-(1.0*1j/hbar)*H(t))*ket(a,b))[0,0]
    # A2=((np.identity(2)-Jdag()*J()/2.0-(1.0*1j/hbar)*H(t))*ket(a,b))[1,0]
    # B1=(J()*ket(a,b))[0,0]
    # B2=(J()*ket(a,b))[1,0]
    # a=a+A1*dt+B1*dW
    # b=b+A2*dt+B2*dW
    # return a,b
   
# std::complex<double> Numero (1.0,2.0);
# std::complex<double> Dos (1.0,3.0);

Numero=1.0+2.0*1j
Dos=1.0+3.0*1j
print euler(Numero,Dos,10,0.1,0.02)
