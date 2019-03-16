# Implementation of the examle problem for 
#G. Fang, S. Liu, M. Fu, B. Wang, Z. Wu, J. Liang, A method to couple state-based peridynamics and
#finite element method for crack propagation problem, Mechanics Research Communications 95 (2019)
#89 -- 95.
#@author patrickdiehl@lsu.edu
#@date 02/09/2019
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, cm
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rcParams.update({'font.size': 15})# Implementation of the examle problem for 

#material properties
K = 2
E=1
mu = (9*K*E)/(9*K-E)

print K,E,mu

h = 0.1
delta=2*h
Area=1
V=V=h*Area

F=1

x = 0.5*(K-mu)
y = 6.*mu/(np.pi*h*np.power(delta,4))
z = 2./(np.pi*h*np.power(delta,3))

print x,y,z


a = E*Area/h
b = 4*delta*y*V*V/h
c = 2*delta*delta*z*z*x*V*V*V/h/h

print a,b,c


#Define the neighborhoods
m = np.array([
    [], #2
    [], #3
    [2,3,5,6], #4
    [4,6,7], #5
    [4,5,7,8], #6
    [5,6,8,9], #7
    [6,7,9,10], #8
    [7,8,10,11], #9
    [8,9,11,12], #10
    [9,10,12], #11
    [10,11,13,14], #12
    [], #13
    [], #14
    ]
    )


#Assemble the stiffness matrix

MCoupled = np.zeros([17,17])

#add the entries for the FEM elements on the left-hand side
MCoupled[0][0] = a
MCoupled[0][1] = -a
MCoupled[1][0] = -a
MCoupled[1][1] = 2*a
MCoupled[1][2] = -a
MCoupled[2][1] = -a
MCoupled[2][2] = 2*a
MCoupled[2][3] = -a
#add the fem entries in the left-hand side coupling zone
MCoupled[3][2] = -a
MCoupled[3][3] = 2*a
MCoupled[3][4] = -a
#add the fem entries in the right-hand side coupling zone
MCoupled[13][12] = -a
MCoupled[13][13] = 2*a
MCoupled[13][14] = - a
#add the entries for the FEM elements on the right-hand side
MCoupled[14][13] = -a
MCoupled[14][14] = 2*a
MCoupled[14][15] = -a
MCoupled[15][14] = -a
MCoupled[15][15] = 2*a
MCoupled[15][16] = -a
MCoupled[16][15] = -a
MCoupled[16][16] = a


#add the pd interactions
for i in range(0,len(m)):
    for j in range(0,len(m[i])):
        idI = 2+i
        idJ = m[i][j]
        MCoupled[idI][idI] += b
        MCoupled[idI][idJ] += -b
        for k in range(0,len(m[i])):
            MCoupled[idI][idI] += c
            MCoupled[idI][m[i][k]] += -c
        #print idJ-4
        for l in range(0,len(m[idJ-2])):
            MCoupled[idI][idJ] += -c
            MCoupled[idI][m[idJ-2][l]] += c

#plot the stiffness matrix
plt.imshow(MCoupled, cmap=cm.binary)
plt.axis('off')
plt.colorbar()
plt.savefig("paper4_matrix_coupled.pdf")
plt.clf()


#set boundary conditions
MCoupled[8][8] = 0 

# generate the left-hand side for the coupling approach
lenM = np.shape(MCoupled)[0]
f = np.zeros(lenM)
f[0]=-F
f[len(f)-1]=F

# solve coupled approach
ucoupled = np.linalg.solve(MCoupled,f)

# generate the position in the bar
xCoupled= np.arange(0,h*17,h)


# Prepare the stiffness matrix for the pure fem approach

hFem = 0.16
a = E*Area/hFem

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
lenMFem = np.shape(MFem)[0]
f = np.zeros(lenMFem)
f[0]=-F
f[len(f)-1]=F
 
# solve fem approach
ufem  = np.linalg.solve(MFem,f) 
 
# generate the position in the bar
xFem= np.arange(0,11*hFem,hFem)

# plot the results
plt.plot(xCoupled,ucoupled,label="Coupled",lw=2)
plt.plot(xFem,ufem,label="FEM",lw=2)
#lt.plot(x,upd,label="PD",lw=2)
plt.legend()
plt.grid()
plt.xlabel("Node position")
plt.ylabel("Displacement")
plt.savefig("paper4.pdf")
