import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

time_steps = []
temp_rows = []
time_list = []
depth_list = []
#Access DX from state.dat
with open("state.dat") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue

        parts = line.split()
        #GET DX
        if (parts[0] == "TIME"):
            DX = float(parts[3][1:])
depth_list.append(0)

with open("temp_diag.dat") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue

        parts = line.split()

        # skip headers
        if parts[0] in ["STATE", "TIME"]:
            continue

        # time row
        try:
            t = float(parts[0])
            temps = list(map(float, parts[1:]))

            time_list.append(t)
            temp_rows.append(temps)
            depth_list.append(depth_list[-1] + DX)

        except:
            continue

# matrix: (time, depth)
data = np.array(temp_rows)

# transpose → (depth, time)
data = data.T

# min/max for legend
t_min = np.min(data)
t_max = np.max(data)

# continuous red → blue colormap
cmap = LinearSegmentedColormap.from_list(
    "blue_red",
    ["blue", "red"]
)

plt.figure(figsize=(10, 6))

times = np.array(time_list, dtype=float)
depth = np.array(depth_list, dtype=float)
im = plt.imshow(
    data,
    aspect='auto',
    origin='upper',
    cmap=cmap,
    extent=[
        times[0], times[-1],   # x-axis (time)
        depth[-1], depth[0]     # y-axis (depth)
]
)

plt.xlabel("Time (seconds)")
plt.ylabel("Depth (cm)")
plt.title("Temperature Evolution Heatmap")

# colorbar with meaning
cbar = plt.colorbar(im,pad=0.1)
cbar.set_label("Temperature")

# legend showing extremes
plt.text(
    1.02, 0.95,
    f"Hot (red): {t_max:.2f}\nCold (blue): {t_min:.2f}",
    transform=plt.gca().transAxes,
    fontsize=10,
    va='top',
    bbox=dict(facecolor='white', alpha=0.7)
)

plt.show()