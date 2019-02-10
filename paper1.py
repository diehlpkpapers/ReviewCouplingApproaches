import numpy as np
import matplotlib.pyplot as plt

#define properties
E=1
Area=1
F=1.0

#define the material constant for the finte element discretization
h=0.1
a=E*Area/h

#define the material constant for the peridyanmic discretization 
delta=2*h
V=h*Area
c=(2*E)/(Area*delta*delta)
b=c/h*V*V

# coupled approach

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


# finite element

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

# peridyanmic

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

#generate the position in the bar
x= np.arange(0,h*11,0.1)

#solve coupled approach
ucoupled = np.linalg.solve(MCoupled,f)

#solve fem approach
ufem  = np.linalg.solve(MFem,f)

#solve pd approach
upd  = np.linalg.solve(MPeridynamics,f)


plt.plot(x,ucoupled,label="Coupled")
plt.plot(x,ufem,label="FEM")
plt.legend()
plt.grid()
plt.xlabel("Node position")
plt.ylabel("Displacement")
plt.savefig("plot_two.pdf")
plt.plot(x,upd,label="PD")
plt.savefig("plot_all.pdf")
