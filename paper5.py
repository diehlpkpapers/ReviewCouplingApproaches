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

x = []
u = []
with open('paper5.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
       x.append(float(row[0]))
       u.append(float(row[1]))

plt.plot(x,u,label="Coupling approach",lw=2)
plt.grid()
plt.xlabel("Node position")
plt.ylabel("Displacement")
plt.legend()
plt.savefig("paper5.pdf")

