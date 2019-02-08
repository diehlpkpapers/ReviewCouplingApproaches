import numpy as np
import matplotlib.pyplot as plt

#define properties
E=1
h=0.1
delta=2*0.1
V=0.1*0.1
A=1
F=0.1
#define the material constant for the finte element discretization
a=E/h

#define the material constant for the peridyanmic discretization 
c=(2*E)/(A*delta*delta)
b=c/h*V*V




# coupled approach

MCoupled = np.array(
 [ [a , -a , 0 , 0 , 0, 0 , 0 , 0 , 0 , 0 , 0  ],
  [  -a , 2*a , -a , 0 , 0 , 0, 0 , 0 , 0 , 0 , 0 ],
  [  0  , -a , 2*a , -a , 0 , 0 , 0, 0 , 0 , 0 , 0  ],
  [      0 , b/4 , -b , 2.5*b , -b , b/4 , 0 , 0 , 0 , 0 , 0 ],
  [       0 , 0 ,  b/4 , -b , 2.5*b ,  -b , b/4 , 0 , 0 , 0 , 0 ],
  [     0 , 0 , 0 ,  b/4 , -b , 2.5*b , -b , b/4 , 0 , 0 , 0 ],
  [      0 , 0 , 0 , 0 , b/4 , -b , 2.5*b , -b , b/4 , 0 , 0 ],
  [      0 , 0 , 0 , 0 , 0 , b/4 , -b , 2.5*b , -b , b/4 , 0 ],
  [      0 , 0 , 0 , 0 , 0 , 0 , 0 , -a , 2*a , -a , 0 ],
  [     0 , 0 , 0 , 0 , 0 , 0 , 0 , 0  , -a , 2*a , -a ],
  [     0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , -a , a ]])


# finite element

MFem = np.array(
 [ [a , -a , 0 , 0 , 0, 0 , 0 , 0 , 0 , 0 , 0  ],
  [  -a , 2*a , -a , 0 , 0 , 0, 0 , 0 , 0 , 0 , 0 ],
  [  0  , -a , 2*a , -a , 0 , 0 , 0, 0 , 0 , 0 , 0  ],
  [      0 , 0 , -a , 2*a , -a , 0 , 0 , 0 , 0 , 0 , 0 ],
  [       0 , 0 ,  0 , -a , 2*a ,  -a , 0 , 0 , 0 , 0 , 0 ],
  [     0 , 0 , 0 ,  0 , -a , 2*a , -a , 0 , 0 , 0 , 0 ],
  [      0 , 0 , 0 , 0 , 0 , -a , 2*a , -a , 0 , 0 , 0 ],
  [      0 , 0 , 0 , 0 , 0 , 0 , -a , 2*a , -a ,0 , 0 ],
  [      0 , 0 , 0 , 0 , 0 , 0 , 0 , -a , 2*a , -a , 0 ],
  [     0 , 0 , 0 , 0 , 0 , 0 , 0 , 0  , -a , 2*a , -a ],
  [     0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , -a , a ]])

# peridyanmic

 

# generate the left-hand side
lenF = np.shape(MFem)[0]
f = np.zeros(lenF)
f[0]=-F
f[len(f)-1]=F

#generate the position in the bar
x= np.arange(0,h*11,0.1)

#solve coupled approach
ucoupled = np.linalg.solve(MCoupled,f)

#solve fem approach
ufem = ucoupled = np.linalg.solve(MFem,f)

print ucoupled

plt.plot(x,ucoupled,label="Coupled")
plt.plot(x,ufem,label="FEM")
plt.legend()
plt.grid()
plt.xlabel("Node position")
plt.ylabel("Displacement")
plt.show()
