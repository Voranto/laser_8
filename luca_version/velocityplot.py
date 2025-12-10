import numpy as np
import matplotlib.pyplot as plt

# Replace 'data.dat' with your filename
filename = "velocity.dat"

# Load data: skip the header line
dataVel = np.loadtxt(filename, skiprows=1)

dataDep1 = np.loadtxt("depth1.dat", skiprows=1)

dataDep2 = np.loadtxt("depth2.dat", skiprows=1)



fig, axes = plt.subplots(2, 2, figsize=(8, 6))

# Split into columns
time = dataVel[:, 0]
velocity = dataVel[:, 1]

axes[0,0].plot(time,velocity, marker='o', linestyle='-')
axes[0,0].set_xlabel('Time (s)')
axes[0,0].set_ylabel('Velocity')
axes[0,0].set_title('Velocity vs Time')
axes[0,0].grid(True)
axes[0,0].set_ylim(-8000,8000)




axes[1,1].plot(dataDep2[:,0],dataDep2[:,1],color='blue',label='Depth 1', marker='o', linestyle='-')
axes[1,1].plot(dataDep1[:,0],dataDep1[:,1],color='red',label='Depth 2', marker='o', linestyle='-')
axes[1,1].set_xlabel('Time (s)')
axes[1,1].set_ylabel('Depth2 and Depth1')
axes[1,1].set_title('Depth2 vs Time')
plt.legend()
axes[1,1].grid(True)

plt.tight_layout()

plt.show()

