import numpy as np
import td


def heat_step_exact(xi, t, a,b, N=50):
  sol = np.zeros(np.size(xi))
  for n in range(1,N+1):
    An = 2/(np.pi * n) * (np.cos(a*np.pi*n) - np.cos(b*np.pi*n))
    lamn = np.pi * n 
    sol += An * np.exp(-lamn * lamn * t ) * np.sin( lamn * xi)
  return sol

def heat (M, N, f, g, k, scheme='explicit'):
  h = 1./np.float (N)
  lbda = k/(h*h)
  a = 2. * lbda
  b = c = np.zeros(N-2) + lbda

  id = np.ones(N-1)
  u = np.zeros((M+1, N-1))

  u[0, :]  = g

  xi = np.array(range(1,N)) * h
  ti = np.array(range(0,M+1)) * k

  if scheme == 'explicit':
    for j in range (0,M):
      u[j+1, :] = td.dot (id-a,b,c,u[j, :]) - k*f[j]
  elif scheme == 'implicit':
    for j in range (0,M):
      u[j+1, :] = td.solve (id+a,-b,-c,u[j, :]) - k*f[j+1]
  elif scheme == 'cn':
    for j in range (0,M):
      u[j+1, :] = td.solve (id+a/2,-b/2,-c/2, \
                  td.dot(id-a/2,b/2,c/2, u[j,:])) - \
                  (0.5*k*f[j+1] + 0.5*k*f[j])
  else:
    print ('unsupported scheme')

  # pack up
  u = np.concatenate ((np.zeros(M+1).reshape(M+1,1), 
                       u, 
                       np.zeros(M+1).reshape(M+1,1)), axis=1 )
  xi = np.concatenate(([0], xi, [1]))

  return ti, xi, u

def hatfunc(N):
  g=np.zeros(N-1)
  a1=np.int((N/4)-1. )
  a2=np.int(3*(N/4)-1.)
  # the range function in python does not inclued the last point
  g[a1] = 0.5
  g[a2] = 0.5
  g[range(a1+1,a2)] = 1
  return(g)



if __name__ == '__main__':

  import itertools
  mc=itertools.cycle('o^vp')

  tN = 0.05
  N=12

  # problem1 : hat function as initial condition
  g = hatfunc(N)
  g_full = np.concatenate(([0], g, [0]))

  xexact = np.linspace(0,1,200)
  texact = tN

  uexact = heat_step_exact (xexact, texact, 0.25, 0.75)

  
  import matplotlib.pyplot as plt
  import pltpref

  for nl, lmbda in enumerate((0.5, 1.5, 2.5)):
    plt.clf()
    h = 1./np.float(N)
    k = lmbda*h*h
    ## corrects lmbda such that integer nr of timesteps
    M = np.int(tN/k)
    k = tN / np.float(M)
    ## and corrects lmbda accordingy
    lmbda = k / (h*h)
    f = np.zeros(M+1)

    plt.plot (xexact, uexact, label='exact', lw=3)
    plt.title("lambda = %0.3f; tN= %0.3f" % (lmbda, tN))
    plt.xlabel("x"); plt.ylabel("$\\bf{u}(tN)$")    
    for ns,scheme in enumerate(('implicit', 'cn','explicit')):

      ti, xi,  u = heat (M, N, f, g, k, scheme=scheme)
      utN = u[u.shape[0]-1, ]
      if (ns==0): plt.plot (xi, g_full,marker='p', label='g', ls='-')
      plt.plot (xi, utN, mc.next(), label=scheme, ls='--')
      if (ns==1): 
        plt.legend( prop={"size":9})
        plt.savefig('heat_xl_'+nl.__str__()+'.pdf')
    plt.legend( prop={"size":9})
    plt.savefig('heat_'+nl.__str__()+'.pdf')


  from mpl_toolkits.mplot3d.axes3d import Axes3D
  from matplotlib import cm

  k=0.005
  N=40
  M=36
  g=hatfunc(N)
  f=np.zeros(M + 1)
  ti, xi,  u = heat (M, N, f, g, k, scheme='implicit')
  fig = plt.figure()
  ax = fig.gca(projection='3d')
  x, t = np.meshgrid(xi,ti)

  surf = ax.plot_surface(
      x, t, u, rstride=1, cstride=1,
      cmap=cm.jet,
      linewidth=0, antialiased=False)


  for tick in ax.xaxis.get_major_ticks():
      tick.label1.set_fontsize(8)
  for tick in ax.yaxis.get_major_ticks():
      tick.label1.set_fontsize(8)
  for tick in ax.zaxis.get_major_ticks():
      tick.label1.set_fontsize(8)

  plt.xlabel("time"); plt.ylabel("x")
  fig.savefig('heat_3d.pdf')
