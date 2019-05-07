# Implementation of the examle problem for 
# T. Ni, M. Zaccariotto, Q.Z. Zhu, U. Galvanetto, Coupling of fem and ordinary
#state-based peridynamics for brittle failure analysis in 3d, Mechanics of 
#Advanced Materials and Structures 0(0), 1 (2019). 
#@author patrickdiehl@lsu.edu
#@date 05/07/2019
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, cm
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rcParams.update({'font.size': 15})
from cycler import cycler
monochrome = (cycler('color', ['k']) * cycler('linestyle', ['-', '--', ':', '=.']))


pos = []
fem = []
bpd = []
spd = []

with open('paper6_coarse.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
       pos.append(float(row[0]))
       fem.append(float(row[1]))
       bpd.append(float(row[2]))
       spd.append(float(row[3]))




# plot the results
fig, ax = plt.subplots(1,1)
ax.set_prop_cycle(monochrome)

ax.plot(pos,fem,label="FEM",lw=2)
ax.plot(pos,bpd,label="Coupled BB",lw=2)
ax.plot(pos,spd,label="Coupled SB",lw=2)
ax.legend()
ax.grid()
plt.xlabel("Node position")
plt.ylabel("Displacement")
plt.tight_layout()
plt.savefig("paper6_coarse.pdf")
