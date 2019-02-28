# Implementation of the examle problem for 
# M. Zaccariotto, T. Mudric, D. Tomasi, A. Shojaei, U. Galvanetto, Coupling of fem meshes with peri-
# dynamic grids, Computer Methods in Applied Mechanics and Engineering 330 (2018) 471-497.
#@author patrickdiehl@lsu.edu
#@date 02/09/2019
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rcParams.update({'font.size': 15})

#define properties
E=1
Area=1
F=1.0

#define the material constant for the finte element discretization
h=0.2409090909090909
a=E*Area/h

#define the material constant for the peridyanmic discretization 
delta=2*h
V=h*Area
c=(2*E)/(Area*delta*delta)
b=c/h*V*V

# Coupling approach presented in the paper

MCoupled = np.array(
 [ [1 , 0 , 0 , 0 , 0, 0 , 0 , 0 , 0 , 0 , 0  ],
  [  -a , 2*a , -a , 0 , 0 , 0, 0 , 0 , 0 , 0 , 0 ],
  [   0 , -a , 2*a , -a , 0 , 0 , 0, 0 , 0 , 0 , 0  ],
  [   0 , -b/4 , -b , 2.5*b , -b , -b/4 , 0 , 0 , 0 , 0 , 0 ],
  [   0 , 0 ,  -b/4 , -b , 2.5*b ,  -b , -b/4 , 0 , 0 , 0 , 0 ],
  [   0 , 0 , 0 ,  -b/4 , -b , 2.5*b , -b , -b/4 , 0 , 0 , 0 ],
  [   0 , 0 , 0 , 0 , -b/4 , -b , 2.5*b , -b , -b/4 , 0 , 0 ],
  [   0 , 0 , 0 , 0 , 0 , -b/4 , -b , 2.5*b , -b , -b/4 , 0 ],
  [   0 , 0 , 0 , 0 , 0 , 0 , 0 , -a , 2*a , -a , 0 ],
  [   0 , 0 , 0 , 0 , 0 , 0 , 0 , 0  , -a , 2*a , -a ],
  [   0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , -a , a ]])


# Pure finite element approach

MFem = np.array(
 [ [1 , 0 , 0 , 0 , 0, 0 , 0 , 0 , 0 , 0 , 0  ],
  [  -a , 2*a , -a , 0 , 0 , 0, 0 , 0 , 0 , 0 , 0 ],
  [   0 , -a , 2*a , -a , 0 , 0 , 0, 0 , 0 , 0 , 0  ],
  [   0 , 0 , -a , 2*a , -a , 0 , 0 , 0 , 0 , 0 , 0 ],
  [   0 , 0 ,  0 , -a , 2*a , -a , 0 , 0 , 0 , 0 , 0 ],
  [   0 , 0 , 0 ,  0 , -a , 2*a , -a , 0 , 0 , 0 , 0 ],
  [   0 , 0 , 0 , 0 , 0 , -a , 2*a , -a , 0 , 0 , 0 ],
  [   0 , 0 , 0 , 0 , 0 , 0 , -a , 2*a , -a ,0 , 0 ],
  [   0 , 0 , 0 , 0 , 0 , 0 , 0 , -a , 2*a , -a , 0 ],
  [   0 , 0 , 0 , 0 , 0 , 0 , 0 , 0  , -a , 2*a , -a ],
  [   0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , -a , a ]])

# Pure peridyanmic approach

MPeridynamics = np.array(
 [ [1, 0 , 0 , 0 , 0, 0 , 0 , 0 , 0 , 0 , 0  ],
  [  -b ,  2*b , -b , 0 , 0 , 0, 0 , 0 , 0 , 0 , 0 ],
  [  -b/4 , -b , 2.5*b , -b , -b/4 , 0 , 0, 0 , 0 , 0 , 0  ],
  [   0 , -b/4 , -b , 2.5*b , -b , -b/4 , 0 , 0 , 0 , 0 , 0 ],
  [   0 , 0 ,  -b/4 , -b , 2.5*b ,  -b , -b/4 , 0 , 0 , 0 , 0 ],
  [   0 , 0 , 0 , -b/4 , -b , 2.5*b ,  -b , -b/4 , 0 , 0 , 0 ],
  [   0 , 0 , 0 , 0 , -b/4 , -b , 2.5*b , -b , -b/4 , 0 , 0 ],
  [   0 , 0 , 0 , 0 , 0 , -b/4 , -b , 2.5*b , -b , -b/4 , 0 ],
  [   0 , 0 , 0 , 0 , 0 , 0 , -b/4 , -b , 2.5*b , -b , -b/4 ],
  [   0 , 0 , 0 , 0 , 0 , 0 , 0 , 0  , -b , 2*b , -b ],
  [   0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , -b , b ]])

 
# generate the left-hand side
lenF = np.shape(MFem)[0]
f = np.zeros(lenF)
f[0]=0
f[len(f)-1]=F

# generate the position in the bar
x= np.arange(0,h*11,0.1)

# solve coupled approach
ucoupled = np.linalg.solve(MCoupled,f)

# solve fem approach
ufem  = np.linalg.solve(MFem,f)

# solve pd approach
upd  = np.linalg.solve(MPeridynamics,f)

# plot the results
plt.plot(x,ucoupled,label="Coupled",lw=2)
plt.plot(x,ufem,label="FEM",lw=2)
plt.plot(x,upd,label="PD",lw=2)
plt.legend()
plt.grid()
plt.xlabel("Node position")
plt.ylabel("Displacement")
plt.savefig("plot_all.pdf")
