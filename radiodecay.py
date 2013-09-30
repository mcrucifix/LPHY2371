# -*- coding: utf-8 -*-
import numpy as np


def F(x): return(-x)

def initX(k,t0,tN, dtype='float'):
  roundoff = 1.e-15
  t = np.arange(t0, tN + roundoff, k, dtype=dtype)
  N = np.size(t)
  X = np.zeros(N, dtype=dtype)
  return t,N,X

def leapfrog(k,t0,tN,X0=1,F=F, dtype='float'):
  t,N,X = initX(k,t0,tN, dtype=dtype)

  X[0] = X0 
  X[1] = 1+k*F(1.)

  for i in range(1,(N-1)):
    X[i+1] = X[i-1] + 2 * k * F(X[i])

  return t,X

def eulerForward(k,t0,tN,X0=1,F=F, dtype='float'):
  t,N,X = initX(k,t0,tN, dtype=dtype)
  X[0] = X0
  for i in range(0,N-2):
    X[i+1] = X[i] * (1 - k)
  return t,X


def eulerBackward(k,t0,tN,X0=1,F=F, dtype='float'):
  t,N,X = initX(k,t0,tN, dtype=dtype)
  X[0] = X0
  for i in range(0,N-1):
    X[i+1] = X[i] / (1 + k)
  return t,X


def eulerTrapez(k,t0,tN,X0=1,F=F, dtype='float'):
  t,N,X = initX(k,t0,tN, dtype=dtype)
  X[0] = X0
  for i in range(0,N-1):
    X[i+1] = X[i] * (1 - k/2.) / (1 + k/2.)
  return t,X


def eulerForward_A(k,t0,tN,X0=1, LastPointOnly = False, dtype='float'):
  if LastPointOnly == True:
    t = tN
    N = int((tN-t0)/k)
    X = np.power(1-k,N)*X0
  else:
    t,N,X = initX(k,t0,tN)
    i = range(0,N)
    X = np.power(1-k,i)*X0
  return t,X


def eulerBackward_A(k,t0,tN,X0=1, LastPointOnly = False, dtype='float'):
  if LastPointOnly == True:
    t = tN
    N = int((tN-t0)/k)
    X = np.power(1+k, -N)*X0
  else:
    t,N,X = initX(k,t0,tN, dtype=dtype)
    i = range(0,N)
    X = np.power((1/(1+k)), i)*X0
  return t,X


def eulerTrapez_A(k,t0,tN,X0=1., LastPointOnly = False, dtype='float' ):
  if LastPointOnly == True:
    t = tN
    N = int((tN-t0)/k)
    X = np.power((1-k/2)/(1+k/2), N)*X0
  else:
    t,N,X = initX(k,t0,tN, dtype=dtype)
    i = range(0,N)
    X = np.power((1-k/2)/(1+k/2), i)*X0
  return t,X

def scheme_error (t0, tN, scheme, dtype='float'):
  lS = np.array([4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768])/(tN-t0)
  Xi = np.zeros(np.size(lS))
  I  = range(0,np.size(lS)-1)

  for i in I:
    t,X = scheme(1/lS[i], t0, tN, LastPointOnly=True, dtype=dtype)
    Xi[i] = np.abs(X - np.exp ( - t ))
  return lS, Xi 

def RK4(k,t0,tN,X0=1., LastPointOnly = False, dtype='float'):
  t,N,X = initX(k,t0,tN, dtype=dtype)
  X[0] = X0
  for i in range(0,N-1):
    Xi = X[i]
    k1 = - k*Xi
    k2 = - k*(Xi+k1/2.)
    k3 = - k*(Xi+k2/2.)
    k4 = - k*(Xi+k3)
    X[i+1] = Xi + (k1+2*k2+2*k3+k4)/6.
  if LastPointOnly == True:
    return t[N-1],X[N-1]
  else:
    return t,X



 
if __name__=="__main__":

  import matplotlib
  import matplotlib.pyplot as plt
  import itertools
  # define scheme names

  leapfrog.frenchName = 'Leap-Frog'
  eulerTrapez.frenchName = u'Centré'
  eulerBackward.frenchName = u'Arrière'
  eulerForward.frenchName = 'Avant'
  RK4.frenchName='RK4'
  eulerForward_A.frenchName = eulerForward.frenchName
  eulerBackward_A.frenchName = eulerBackward.frenchName
  eulerTrapez_A.frenchName = eulerTrapez.frenchName

  # plot trajectory in time

#  matplotlib.rc('text', usetex=True)
  matplotlib.rc('font', size=11)
  import pltpref
  k = 0.3
  t0 = 0.
  tN = 5.
  X0 = 1.

  tA = np.linspace(t0,tN,50)
  Anal = np.exp(-tA)
  plt.plot(tA, Anal, '-', label='exact')

  mc=itertools.cycle('o^vp')

  for scheme in (eulerBackward,eulerForward,eulerTrapez):
    L = scheme(k,t0,tN,X0,F)
    plt.plot(L[0], L[1], marker=mc.next(),ls='None', label=scheme.frenchName )
  plt.xlabel("$t$")
  plt.ylabel("$x$")
  plt.title(u"Décroissance radioactive; $k=0.3$")
  plt.legend( loc = 'upper right')
  plt.savefig("radio.pdf")

  L = leapfrog(k,t0,tN,X0,F)
  plt.plot(L[0], L[1], marker=mc.next(), ls='None')
  plt.legend( loc= 'lower left')
  plt.savefig("radio_leapfrog.pdf")


  # 
  plt.clf()
  tN = 0.5
  mc=itertools.cycle('o^vp')
  for method in (eulerBackward_A, eulerForward_A, eulerTrapez_A):
    L,S = scheme_error(t0, tN, method)
    plt.plot(L, S, marker=mc.next(),ls='None', label=method.frenchName )
  plt.title ('$|x(t) - x_{t/k}|,\ t=1$')
  plt.xlabel("$1/k$")
  plt.xscale('log')
  plt.yscale('log')
  plt.legend( loc= 'lower left')
  plt.axhline(ls="--", color="green")
  plt.savefig("error.pdf")

  L,S = scheme_error(t0, tN, RK4, dtype='float')
  plt.plot(L, S, marker=mc.next(),ls='None', label='RK4' )
  plt.legend( loc= 'lower left')
  plt.savefig("error_withRK4.pdf")

## rounding error graphic
  plt.clf()
  plt.xlabel("$1/k$")
  for type in ('float16','float32','float64'):
    L,S = scheme_error(t0, tN, RK4, dtype=type)
    plt.plot(L, S, marker=mc.next(),ls='None', label=type )
  plt.xscale('log')
  plt.yscale('log')
  plt.title ('$|x_{t/k} - \overline{x}_{t/k}|,\ t=1$')
  plt.legend( loc= 'lower left')
  plt.savefig("rounding_error.pdf")
