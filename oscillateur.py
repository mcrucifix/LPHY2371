# -*- coding: utf-8 -*-
import numpy as np
from numpy.linalg import solve
import matplotlib.pyplot as plt
import itertools

class PMatrix (np.matrix):
  def energy(self):
   return (np.dot(np.power(self,2),np.ones(2))/2.).T
  def P(self):
   return (self[:,1])
  def Q(self):
   return (self[:,2])
  def view(self): np.view(self)

M = np.matrix('0 -1; 1 0')

def initX(k,t0,tN):
  t = np.arange(t0, tN , k)
  N = np.size(t)
  P = np.zeros((N,2))
  return t,N,P

def eulerForward(k,t0,tN,P0=(1,0),M=M):
  t,N,P = initX(k,t0,tN)
  P[0,] = P0
  for i in range(0,N-1):
    P[i+1,] =  np.dot((np.eye(2) + M*k),P[i,])
  return t,P.view(PMatrix)

def eulerBackward(k,t0,tN,P0=(1,0),M=M):
  t,N,P = initX(k,t0,tN)
  P[0,] = P0
  for i in range(0,N-1):
    P[i+1,] =  solve((np.eye(2) - M*k),P[i,])
  return t,P.view(PMatrix)

def eulerTrapez(k,t0,tN,P0=(1,0),M=M):
  t,N,P = initX(k,t0,tN)
  P[0,] = P0
  for i in range(0,N-1):
    tmp =  solve((np.eye(2) - M*k/2),P[i,])
    P[i+1,] =  np.dot((np.eye(2) + M*k/2),tmp)
  return t,P.view(PMatrix)

def RK4(k,t0,tN,P0=(1,0),M=M):
  t,N,P = initX(k,t0,tN)
  P[0,] = P0
  for i in range(0,N-1):
    Pi = np.matrix(P[i,]).T
    k1 = k*np.dot(M,Pi)
    k2 = k*np.dot(M,(Pi+k1/2.))
    k3 = k*np.dot(M,(Pi+k2/2.))
    k4 = k*np.dot(M,(Pi+k3))
    P[i+1,] = P[i,] + (k1+2*k2+2*k3+k4).T/6.
  return t,P.view(PMatrix)

def heun(k,t0,tN,P0=(1,0),M=M):
  t,N,P = initX(k,t0,tN)
  P[0,] = P0
  for i in range(0,N-1):
    tmp =  np.dot((np.eye(2) + M*k),P[i,])
    P[i+1,] =  np.dot((np.eye(2) + M*k/2),P[i,]) \
             + np.dot(           + M*k/2 ,tmp.T ).T 
  return t,P.view(PMatrix)

def verlet(k,t0,tN,P0=(1,0),M=M):
  t,N,P = initX(k,t0,tN)
  P[0,] = P0
  for i in range(0,N-1):
   pi1 = P[i,0] + M[0,1]*P[i,1]*k # predict veloc. first
   qi1 = P[i,1] + M[1,0]*k*(pi1+P[i,0])/2.
   pi1 = P[i,0] + M[0,1]*k*(qi1+P[i,1])/2.
   P[i+1,]=(pi1, qi1)
  return t,P.view(PMatrix)


if __name__=="__main__":

  eulerForward.frenchname=r'Avant'
  eulerBackward.frenchname=u'Arrière'
  eulerTrapez.frenchname=u'Centré'
  heun.frenchname=u'Heun'
  verlet.frenchname=u'Verlet'
  RK4.frenchname=u'RK4'
  
  import pltpref

  plt.clf()

  t0 = 0. 
  k  = 0.2
  tN = 20

  for scheme in (eulerForward, eulerBackward, eulerTrapez):
#  for scheme in (eulerTrapez,verlet,RK4):
   L = scheme(k,t0,tN)
   plt.plot(L[0],L[1].P(), label=scheme.frenchname)


  plt.ylabel("p")
  plt.xlabel("time (AU)")
  plt.legend(loc='lower left' )

  plt.savefig('oscillateur_euler.pdf')
  plt.clf()

  t0 = 0. 
  k  = 0.5
  tN = 10

  for scheme in (eulerTrapez,verlet,RK4, heun):
   L = scheme(k,t0,tN)
   plt.plot(L[0],L[1].energy(), label=scheme.frenchname)

  plt.ylabel("Hamiltonien")
  plt.xlabel("time (AU)")
  plt.legend(loc='upper left' )
  plt.savefig('oscillateur_energy1.pdf')

  plt.clf()
  k  = 0.5 
  tN = 1000
  for scheme in (eulerTrapez,verlet,RK4):
   L = scheme(k,t0,tN)
   plt.plot(L[0],L[1].energy(), label=scheme.frenchname)

  plt.ylabel("Hamiltonien")
  plt.xlabel("time (AU)")
  plt.legend(loc='lower left' )

  plt.savefig('oscillateur_energy2.pdf')


