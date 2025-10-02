import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Set the working directory to the location of the script
file_path = os.path.abspath(__file__)

path_to_csv = os.path.join(os.path.dirname(file_path), 'Trajectories_Front.csv')

# Read the CSV file
I = pd.read_csv(path_to_csv, header=None).values

# Organize the tracks

A = []
B = []
i = -1
Track = []

for k in range(len(I)):
    TrackFrame = I[k, 0]
    if TrackFrame == 0:
        Track = []
        j = 0
        i += 1
    j = len(Track)
    Track.append([I[k, 1], I[k, 2]])
    if len(A) <= i:
        A.append([])
        B.append(0)
    A[i] = np.array(Track)
    B[i] = len(Track)

N = i + 1
C = []

for i in range(N):
    Disp = [0]
    for j in range(1, B[i]):
        d = np.sqrt((A[i][j, 0] - A[i][j-1, 0])**2 + (A[i][j, 1] - A[i][j-1, 1])**2)
        Disp.append(d)
    C.append(Disp)

# Plotting

sns.set(style="ticks")
fig, ax = plt.subplots()
ax.plot(0, 0, 'k')
ax.set_xlabel('X [um]')
ax.set_ylabel('Y [um]')
ax.set_aspect('equal')
ax.set_xlim([-10, 80])
ax.set_ylim([-10, 34])
ax.set_facecolor('k')

for i in range(N):
    j = i  
    Q = j // 5
    R = j % 5
    X = A[i][:, 0] + 7 * Q
    Y = A[i][:, 1] + 7 * R
    for k in range(B[i] - 1):
        if C[i][k+1] < 0.17:
            color = 'magenta'     
        elif C[i][k+1] < 0.23:
            color = 'cyan'      
        elif C[i][k+1] < 0.31:
            color = '#77AC30'       
        elif C[i][k+1] < 0.44:
            color = 'yellow'        
        else:
            color = 'red'        
        ax.plot(X[k:k+2], Y[k:k+2], color=color)

plt.show()
