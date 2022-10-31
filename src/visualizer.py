#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

natoms = 5

start = 0
stop = 1000
step = 1*natoms  # Number of atoms or its multiple

data = np.loadtxt("trajectory.dat")

for i in range(natoms):
    plt.plot(data[start+i:stop+i:step, 2], data[start +
             i:stop+i:step, 3],  label="atom " + str(i+1))
plt.legend()
plt.show()
