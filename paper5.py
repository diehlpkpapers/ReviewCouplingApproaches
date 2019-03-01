# Implementation of the examle problem for 
#   W. Sun, J. Fish, Superposition-based coupling of peridynamics and finite 
#   element method, Computational Mechanics (2019) 1--18
#@author patrickdiehl@lsu.edu
#@date 02/26/2019
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rcParams.update({'font.size': 15})
plt.rcParams.update({'figure.autolayout': True})

#define properties
E=200e9
Area=1e-4
F=1.0

#define the material constant for the finte element discretization
h=0.36
a=E*Area/h

# Pure finite element approach

MFem = np.array(
 [ [a , -a , 0 , 0 , 0, 0 , 0 , 0 , 0 , 0 , 0  ],
  [  -a , 2*a , -a , 0 , 0 , 0, 0 , 0 , 0 , 0 , 0 ],
  [   0 , -a , 2*a , -a , 0 , 0 , 0, 0 , 0 , 0 , 0  ],
  [   0 , 0 , -a , 2*a , -a , 0 , 0 , 0 , 0 , 0 , 0 ],
  [   0 , 0 ,  0 , -a , 2*a , -a , 0 , 0 , 0 , 0 , 0 ],
  [   0 , 0 , 0 ,  0 , -a , 0 , -a , 0 , 0 , 0 , 0 ],
  [   0 , 0 , 0 , 0 , 0 , -a , 2*a , -a , 0 , 0 , 0 ],
  [   0 , 0 , 0 , 0 , 0 , 0 , -a , 2*a , -a ,0 , 0 ],
  [   0 , 0 , 0 , 0 , 0 , 0 , 0 , -a , 2*a , -a , 0 ],
  [   0 , 0 , 0 , 0 , 0 , 0 , 0 , 0  , -a , 2*a , -a ],
  [   0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , -a , a ]])

# generate the left-hand side
lenF = np.shape(MFem)[0]
f = np.zeros(lenF)
f[0]=-F
f[len(f)-1]=F

# generate the position in the bar
xFEM  = np.arange(-0.6,h*9,h)

# solve fem approach
uFEM  = np.linalg.solve(MFem,f)

print uFEM

xCoupled = []
uCoupled = []
with open('paper5.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
       xCoupled.append(float(row[0]))
       uCoupled.append(float(row[1]))
       
       
plt.plot(xCoupled,uCoupled,label="Coupling approach",lw=2)
plt.plot(xFEM,uFEM,label="FEM",lw=2)
plt.grid()
plt.xlabel("Node position")
plt.ylabel("Displacement")
plt.legend()
plt.savefig("paper5.pdf")

