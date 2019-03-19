# Implementation of the examle problem for 
# M. Zaccariotto, T. Mudric, D. Tomasi, A. Shojaei, U. Galvanetto, Coupling of fem meshes with peri-
# dynamic grids, Computer Methods in Applied Mechanics and Engineering 330 (2018) 471-497.
#@author patrickdiehl@lsu.edu
#@date 02/09/2019
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, cm
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rcParams.update({'font.size': 15})

##############################################################################
#define non-linear loading functions
##############################################################################

def f1(x):
    return x;

def f2(x):
    return x*x;

def f3(x):
    return np.sin(2*np.pi*x)
    
def f4(x):
    return np.sqrt(x)

def diff(u,x):
    return np.diff(u) / np.diff(x)

##############################################################################
#define functions
##############################################################################

def solve(M,f):
    return np.linalg.solve(M,f)

def plotu(x,ucoupled,ufem,upd,name,title):
    plt.plot(x,ucoupled,label="Coupled",lw=2)
    plt.plot(x,ufem,label="FEM",lw=2)
    plt.plot(x,upd,label="PD",lw=2)
    plt.legend()
    plt.grid()
    plt.title(title)
    plt.xlabel("Node position")
    plt.ylabel("Displacement")
    plt.savefig(name)
    plt.cla()
    
    
def plotd(x,diffcoupled,difffem,diffpd,name):
    plt.plot(x[:-1],diffcoupled,label="Coupled",lw=2)
    plt.plot(x[:-1],difffem,label="FEM",lw=2)
    plt.plot(x[:-1],diffpd,label="PD",lw=2)
    plt.legend()
    plt.grid()
    plt.title("$u'(x)$")
    plt.xlabel("Node position")
    plt.ylabel("Displacement")
    plt.savefig(name)
    plt.cla()
    
##############################################################################
#define properties
##############################################################################
E=1
Area=1
F=1.0

##############################################################################
#define the material constant for the finte element discretization
##############################################################################
h=0.2409090909090909
a=E*Area/h

##############################################################################
#define the material constant for the peridyanmic discretization 
##############################################################################
delta=2*h
V=h*Area
c=(2*E)/(Area*delta*delta)
b=c/h*V*V

##############################################################################
# Coupling approach presented in the paper
##############################################################################

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

##############################################################################
# Pure finite element approach
#############################################################################

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

##############################################################################
# Pure peridyanmic approach
##############################################################################

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

##############################################################################
#Loading 1
##############################################################################

# generate the position in the bar
x= np.arange(0,h*11,h)


# generate the left-hand side
lenF = np.shape(MFem)[0]
f = np.zeros(lenF)

for i in range(0,4):
    f[i] = -f1(x[i])
    
for i in range(5,lenF):
    f[i] = f1(x[i])

f[0]=0

# solve coupled approach
ucoupled = solve(MCoupled,f)
ufem  = solve(MFem,f)
upd  = solve(MPeridynamics,f)

#compute first derivative
difffem = diff(ufem,x) 
diffpd = diff(upd,x) 
diffcoupled = diff(ucoupled,x)

# plot the results
plotu(x,ucoupled,ufem,upd,"paper1_f1.pdf","$f(x)=x$")
plotd(x,diffcoupled,difffem,diffpd,"paper1_f1_diff.pdf")

##############################################################################
#Loading 2
##############################################################################

# generate the left-hand side
lenF = np.shape(MFem)[0]
f = np.zeros(lenF)

for i in range(0,4):
    f[i] = -f2(x[i])
    
for i in range(5,lenF):
    f[i] = f2(x[i])

f[0]=0

# solve the stiffness matric for all approaches
ucoupled = solve(MCoupled,f)
ufem  = np.linalg.solve(MFem,f)
upd  = np.linalg.solve(MPeridynamics,f)

#compute first derivative
difffem = diff(ufem,x) 
diffpd = diff(upd,x) 
diffcoupled = diff(ucoupled,x)


# plot the results
plotu(x,ucoupled,ufem,upd,"paper1_f2.pdf","$f(x)=x^2$")
plotd(x,diffcoupled,difffem,diffpd,"paper1_f2_diff.pdf")

##############################################################################
#Loading 3
##############################################################################

# generate the left-hand side
lenF = np.shape(MFem)[0]
f = np.zeros(lenF)

for i in range(0,4):
    f[i] = -1*f3(x[i])
    
for i in range(5,lenF):
    f[i] = f3(x[i])

f[0]=0

# solve the stiffness matric for all approaches
ucoupled = solve(MCoupled,f)
ufem  = solve(MFem,f)
upd  = solve(MPeridynamics,f)

#compute first derivative
difffem = diff(ufem,x) 
diffpd = diff(upd,x) 
diffcoupled = diff(ucoupled,x)

# plot the results
plotu(x,ucoupled,ufem,upd,"paper1_f3.pdf","$f(x)=sin(2\pi x)$")
plotd(x,diffcoupled,difffem,diffpd,"paper1_f3_diff.pdf")

##############################################################################
#Loading 4
##############################################################################

# generate the left-hand side
lenF = np.shape(MFem)[0]
f = np.zeros(lenF)

for i in range(0,4):
    f[i] = -1*f4(x[i])
    
for i in range(5,lenF):
    f[i] = f4(x[i])

f[0]=0

# solve coupled approach
ucoupled = solve(MCoupled,f)
ufem  = solve(MFem,f)
upd  = solve(MPeridynamics,f)

#compute first derivative
difffem = diff(ufem,x) 
diffpd = diff(upd,x) 
diffcoupled = diff(ucoupled,x)

# plot the results
plotu(x,ucoupled,ufem,upd,"paper1_f4.pdf","$f(x)=\sqrt(x)$")
plotd(x,diffcoupled,difffem,diffpd,"paper1_f4_diff.pdf")
