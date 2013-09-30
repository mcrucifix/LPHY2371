
import numpy as np
import itertools
import td

# valid for general problem y'' + p y' + q y = f
def bvp_Mat (h, q, p):
  a = - 2. + h * h * q
  b = 1 - h/2 * p 
  c = 1 + h/2 * p
  return a,b,c

def problem1 (N):
  l = 1. 
  h = 1./np.float(N+1)
  q = - np.ones(N)
  p = np.zeros(N)
  xi = np.array(range(1,N+1))*h
  f = -   np.sin(2*np.pi*xi)
  a,b,c = bvp_Mat(h, q, p)
  y = td.solve(a,b,c,f*h*h)

  y  = np.concatenate(([0],y,[0]), axis=0)
  xi = np.concatenate(([0],xi,[1]), axis=0)
  return xi, y

def problem2 (N):
# valid for general problem y'' + p y' + q y = f

  l = 1. 
  epsilon = 0.01
  h = 1./np.float(N+1)
  xi = np.array(range(1,N+1))*h
  p = -  xi*xi / epsilon
  q = - np.ones(N) / epsilon
  a,b,c = bvp_Mat(h, q, p)

  f = np.zeros(N)
#  boundary conditions
  f[0]        =  - 1 * b[0] 
  f[N-1]      =  - 1 * c[N-1]

  y = td.solve(a,b,c,f)

  y  = np.concatenate(([1],y,[1]), axis=0) 
  xi = np.concatenate(([0],xi,[1]), axis=0)
  return xi, y

def problem2_varchange (N):
# valid for general problem y'' + p y' + q y = f
  l = 1. 
  epsilon = 0.01
  h = 1./np.float(N+1)
  xi = np.array(range(1,N+1))*h
  p = -  xi*xi / epsilon
  q = - np.ones(N) / epsilon

  a,b,c = bvp_Mat(h, q, p)

  f =  np.ones(N) * h * h / epsilon

  y = td.solve(a,b,c,f)  + np.ones(N)

  y  = np.concatenate(([1],y,[1]), axis=0) 
  xi = np.concatenate(([0],xi,[1]), axis=0)
  return xi, y


def bratu (N, tol=1.e-6, itmax=20, initscale=1):
   h = 1./np.float(N+1)
   xi = np.array(range(1,N+1))*h
   # solution 1.
   z = initscale*xi*(1-xi)
   F = np.zeros(N)
   h2 = h*h
   b = c = np.ones(N)
   err = 999999.
   it = 0
   while  ( (err > tol) & (it < itmax)):
     # finite difference
     F[1:N-1] = z[2:N]-2*z[1:N-1] + z[0:N-2]
     F[N-1]     =  -2*z[N-1] + z[N-2]
     F[0]     =  z[1]-2*z[0]
     # add the nonlinear term
     F += np.exp(z) * h2
     a = -2. + np.exp(z)*h2
     z1 = td.solve(a,b,c,F)
     err = np.max(abs(z1))
     z -= z1
     it+=1
   if (it == itmax): print('warning: max. it. reached')
   
   xi = np.concatenate(([0],xi,[1]), axis=0)
   z  = np.concatenate(([0],z,[0]), axis=0)

   return xi, z

# valid for general problem y'' + p y' + q y = f

if __name__ == '__main__':
  # problem 1 : example 1 p. 53
  import matplotlib.pyplot as plt 
  import pltpref

  N=100
  h=1./100

  plt.clf()
  xexact = np.array(range(0,N+2))*h
  exact = np.sin(2*np.pi*xexact)/(1+4*np.pi*np.pi)
 
  plt.plot(xexact, exact, '-', label='exacte', linewidth=2)

  mc=itertools.cycle('o^vp')
  for N in (2,4,8,16):
     xi, y = problem1(N)
     plt.plot(xi,y, marker=mc.next(), ls='dashed', label='N='+N.__str__())
  plt.legend() 
  plt.xlabel("x")
  plt.ylabel("y")
  plt.xlim((0,1))
  plt.title("Solution de $y^{\prime\prime} - y = - \sin(2\pi x)$")

  plt.savefig('bvp_prob1.pdf')

## compute scheme error
 
  plt.clf()
  Ns = (8,16,32,64,128,256,512,1024,2048,4096,16384)
  Err = np.zeros(np.size(Ns))
  for i in range(0,np.size(Ns)):
     xi, y = problem1(Ns[i])
     exact = np.sin(2*np.pi*xi)/(1+4*np.pi*np.pi)
     Err[i] = np.max(np.abs(exact - y))

  plt.plot(Ns, Err, '.')
  plt.title('erreur maximale')
  plt.xscale('log')
  plt.yscale('log')
  plt.xlabel('N')
  plt.title("Erreur de $y^{\prime\prime} - y = - \sin(2\pi x)$")
  plt.savefig('bvp_error.pdf')

## problem2

  plt.clf()
  for N in (10,20):
     xi, y = problem2(N)
     plt.plot(xi,y, '--', marker=mc.next(), label='N='+N.__str__())
  N = 120
  xi, y = problem2(N)
  plt.plot(xi,y, '-', label='N='+N.__str__())

  plt.legend() 
  plt.xlabel("x")
  plt.ylabel("y")
  plt.xlim((0,1))
 
  plt.savefig('bvp_prob2.pdf')

## problem2

  plt.clf()
  for N in (10,20):
     xi, y = problem2_varchange(N)
     plt.plot(xi,y, '--', marker=mc.next(), label='N='+N.__str__())
  N = 120
  xi, y = problem2_varchange(N)
  plt.plot(xi,y, '-', label='N='+N.__str__())

  plt.legend() 
  plt.xlabel("x")
  plt.ylabel("y")
  plt.xlim((0,1))
 
  plt.savefig('bvp_prob2_varchange.pdf')


## Bratu's nonlinear equation
  plt.clf()
  N=12
  for initscale in (1,4,8):
      xi,y = bratu(N, initscale=initscale)
      plt.plot(xi,y, '-', label='$\mu$='+initscale.__str__())
  plt.legend(); plt.xlabel("x"); plt.ylabel("y")
  plt.ylim(0,0.2)
  plt.title('Equation de Bratus, avec $z_0=\mu x(1-x)$')
  plt.savefig('bratu.pdf')

   
