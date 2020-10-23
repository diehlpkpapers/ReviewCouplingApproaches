# Implementation of the examle problem for 
# E. Madenci, A. Barut, M. Dorduncu, N.D. Phan, in 2018AIAA/ASCE/AHS/ASC Structures, 
# Structural Dynamics, and Materials Conference (2018)
#@author patrickdiehl@lsu.edu
#@date 10/24/2020
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, cm
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rcParams.update({'font.size': 15})
from cycler import cycler
monochrome = (cycler('color', ['k']) * cycler('linestyle', ['-', '--', ':', '=.']))

def solution(x):
    return 2*x - 1

pos = []
fem = []

with open('../paper8.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
       pos.append(float(row[0]))
       fem.append(float(row[1]))

pos = np.array(pos)


# plot the results
fig, ax = plt.subplots(1,1)
ax.set_prop_cycle(monochrome)



ax.plot(pos,fem,label="Coupled",lw=2)
ax.plot(pos,solution(pos),label="FEM",lw=2)
ax.legend()
ax.grid()
plt.xlabel("Node position")
plt.ylabel("Displacement")
plt.tight_layout()
plt.savefig("paper8.pdf")



