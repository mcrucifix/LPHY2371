import numpy as np

# solve tridiagonal system
# could be replaced by a call to scipy. ??? (need to find out)

#import scipy.linalg as lin

#def solve_s (a,b,c,z):
#  N  = np.size(a)
#  B  = np.concatenate((b[range(1,N)], [0]))
#  C  = np.concatenate(([0], c[range(0,N-1)]))
#  ab = np.matrix([C, a, B])
#  y  = lin.solve_banded ((1,1),  ab, z )
#  return(y)

def solve (a,b,c,z):
  N = np.size(a)
  y = np.zeros(N)
  v = np.zeros(N)
  w = a[0] 
  y[0]=z[0]/w
  for i in range(1, N):
    v[i-1] = c[i-1] / w
    w = a[i] - b[i-1] * v[i-1]
    y[i] = (z[i] - b[i-1]*y[i-1])/w
  for i in range(N-2,-1,-1):
    y[i]  = y[i] - v[i]* y[i+1]
  return y

def dot (a,b,c,z):
  N = np.size(z)
  y = a*z 
  y [0:N-1] += c[0:N-1] * z [1:N]
  y [1:N] += b[0:N-1] * z [0:N-1]
  return y


