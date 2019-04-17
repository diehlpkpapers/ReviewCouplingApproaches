# Implementation of the examle problem for 
#   W. Sun, J. Fish, Superposition-based coupling of peridynamics and finite 
#   element method, Computational Mechanics (2019) 1--18
#@author patrickdiehl@lsu.edu
#@date 02/26/2019
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc , cm
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rcParams.update({'font.size': 15})
plt.rcParams.update({'figure.autolayout': True})
from cycler import cycler
monochrome = (cycler('color', ['k']) * cycler('linestyle', ['-', '--', ':', '=.']))

##############################################################################
#define properties
##############################################################################
E=200e9
Area=1e-4
F=1.0

##############################################################################
#define the material constant for the finte element discretization
##############################################################################
h=0.36
a=E*Area/h

##############################################################################
#define the material constant for the peridyanmic discretization 
##############################################################################
hPD=0.36
delta=2.*hPD
V=hPD*Area
c=(2.*E)/(Area*delta*delta)
b=c/h*V*V

##############################################################################
# Pure finite element approach
##############################################################################

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
 
# Pure PD approach
pdNodes = 11

MPD = np.zeros([pdNodes,pdNodes])

for i in range(0,pdNodes):
    if i == 0:
        MPD[i][0] = b
        MPD[i][1] = -b
    elif i == 1:
        MPD[i][i-1] = -b
        MPD[i][i] = 2.*b
        MPD[i][i+1] = -b
    elif i == pdNodes - 2:
        MPD[i][i-1] = -b
        MPD[i][i] = 2.*b
        MPD[i][i+1] = -b
    elif i == pdNodes - 1:
        MPD[i][i] = b
        MPD[i][i-1] = -b
    else:
        MPD[i][i-2] = -b/4.
        MPD[i][i-1] = -b
        MPD[i][i] = 2.5*b
        MPD[i][i+1] = -b
        MPD[i][i+2] = -b/4.
        
MPD[5][5] = 0
        
# generate the left-hand side
lenF = np.shape(MFem)[0]
f = np.zeros(lenF)
f[0]=-F
f[len(f)-1]=F

# generate the position in the bar
xFEM  = np.arange(-0.6,h*9,h)

# solve fem approach
uFEM  = np.linalg.solve(MFem,f)

# generate the left-hand side for PD
lenFPD = np.shape(MPD)[0]
fPD = np.zeros(lenFPD)
fPD[0]=-F 
fPD[len(fPD)-1]=F


# solve fem approach
uPD  = np.linalg.solve(MPD,fPD)


# generate the position in the bar
xPD  = np.arange(-0.6,3.1,hPD)

xCoupled = []
uCoupled = []
with open('paper5.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
       xCoupled.append(float(row[0]))
       uCoupled.append(float(row[1]))

# plot
fig, ax = plt.subplots(1,1)
ax.set_prop_cycle(monochrome)

ax.plot(xCoupled,uCoupled,label="Coupled",lw=2)
ax.plot(xFEM,uFEM,label="FEM",lw=2)
ax.plot(xPD,uPD,label="PD",lw=2)
plt.grid()
plt.xlabel("Node position")
plt.ylabel("Displacement")
plt.legend()
plt.savefig("paper5.pdf")


mat = np.zeros([52,52])

with open('paper5_stiffness.txt') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=' ')
    i = 0
    for row in csv_reader:
        for j in range(0,len(row)):
            mat[i][j] = float(row[j])
        i = i + 1
        print i
        
plt.cla()

plt.imshow(mat, cmap=cm.binary)
plt.axis('off')
plt.colorbar()

plt.savefig("paper5_matrix.pdf")
