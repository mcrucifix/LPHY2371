from numpy import * 
import td

import pltpref
def uinit(N):
 
  i=int(N/4)
  j=int(3*N/4)
  u0= zeros(N)
  u0[i:j] = 1.
  return(u0)

def uinit_left(N):
 
  i=int(N/10)
  j=int(2*N/10)
  u0= zeros(N)
  u0[i:j] = 1.
  return(u0)




def upwind(lmbda):
  return 0, 1-lmbda, lmbda

def downwind(lmbda):
  return -lmbda, 1+lmbda, 0 

def laxwend(lmbda):
  return -0.5*lmbda*(1-lmbda), 1-lmbda*lmbda, 0.5 * lmbda*(1+lmbda)

def advect(scheme,lmbda,M,N, uinit=uinit):
  a,b,c = scheme(lmbda)

  u = zeros((M,N))
  u[0,:] = uinit(N)

  for i in range(0,M-1):
    u[i+1,range(1,N-2)] = a*u[i,range(2,N-1)] + b*u[i,range(1,N-2)] + c*u[i,range(0,N-3)]
    u[i+1,0]= a*u[i,1]  + b*u[i,0]
    u[i+1,N-1]= c*u[i,N-2]+ b*u[i,N-1]
  return u


M=10
N=40

k=0.15
h=0.1
a=0.5

lmbda=a*k/h
xi = h * array(range(0,N))

import matplotlib.pyplot as plt
for scheme in (upwind,downwind):

  u = advect(scheme,lmbda,M,N)

  plt.clf()
  plt.plot(xi,u[0,:], label='t=0')
  plt.plot(xi,u[M-1,:], label='t='+" %2i " % (k*M))
  plt.xlabel("x")
  plt.ylabel("u")
  plt.ylim((0,1.2))
  plt.title(scheme.func_name + (" $\lambda=$%5.2f" % lmbda))
  plt.legend(loc='best')
  plt.savefig(scheme.func_name+'.pdf')


plt.clf()
plt.plot(xi,u[0,:], '--', label='t=0' )

for M in (5,10,100):

  k=float(1./M)
  lmbda=a*k/h
  scheme = upwind
  u = advect(scheme,lmbda,M,N)


  plt.plot(xi,u[M-1,:], label="$\lambda=$%5.3f" % lmbda)
  plt.xlabel("x")
  plt.ylabel("u")
  plt.ylim((0,1.2))
  plt.title('signal carre, t=1, N=40')
  plt.legend(bbox_to_anchor=(0.25, 1.05))

plt.savefig(scheme.func_name+'_ks.pdf')

plt.clf()

k=0.07
h=0.1
a=1.0
lmbda=a*k/h
M = int(36/k)
N = int(50/h)
xi = h * array(range(0,N))

for scheme in (upwind,laxwend):
  u = advect(scheme,lmbda,M,N, uinit=uinit_left)
  plt.plot(xi,u[M-1,:], label=scheme.func_name)


plt.xlabel("x")
plt.ylabel("u")
#plt.ylim((0,1.2))
plt.title("t= %i , a=%.2f, $\lambda$=%.2f" % ((k*M), a, lmbda))
plt.legend(loc='best')
plt.savefig('compare_up_lax.pdf')



