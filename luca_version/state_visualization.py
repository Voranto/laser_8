state_info = {
    1: ("CRYSTAL", "blue"),
    2: ("LARGE POLY", "orange"),
    3: ("FINE POLY", "yellow"),
    4: ("AMORPHOUS", "purple"),
    5: ("MUSHY", "brown"),
    8: ("LIQUID", "red"),
    9: ("SUPERCOOLED", "cyan"),
}

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

letter_to_state = {
    'C': 1,
    'P': 2,
    'F': 3,
    'A': 4,
    'M': 5,
    'L': 8,
    'S': 9,
}



# --- parse data (your code unchanged) ---
time_steps = []
time_list = []
depth_list = []

with open("state.dat") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue

        parts = line.split()
        #GET DX
        if (parts[0] == "TIME"):
            print(parts)
            DX = float(parts[3][1:])
            print(DX)

        if all(len(p) == 1 and p.isalpha() for p in parts):
            depth_list.append(DX)
            time_steps.append(parts)

        elif len(parts) > 1 and parts[1].isalpha():
            time_steps.append(parts[1:])
            time_list.append(float(parts[0]))
            depth_list.append(depth_list[-1] + DX)

data = np.array(time_steps)

numeric = np.vectorize(letter_to_state.get)(data).astype(float)
numeric = numeric.T


values = sorted(state_info.keys()) 

# map states into contiguous indices 0..6
value_to_index = {v: i for i, v in enumerate(values)}

indexed = np.vectorize(lambda x: value_to_index.get(x, -1))(numeric)

colors = [state_info[v][1] for v in values]
cmap = ListedColormap(colors)

# boundaries so colors align correctly
norm = BoundaryNorm(np.arange(-0.5, len(values) + 0.5, 1), cmap.N)

# --- plot ---
times = np.array(time_list, dtype=float)
depth = np.array(depth_list, dtype=float)

plt.imshow(indexed, cmap=cmap, norm=norm, aspect='auto', origin='upper', extent=[
        times[0], times[-1],   # x-axis (time)
        depth[-1], depth[0]     # y-axis (depth)
    ])
plt.xlabel("Time(seconds)")
plt.ylabel("Depth(cm)")
plt.title("State evolution")

import matplotlib.patches as mpatches

patches = [
    mpatches.Patch(color=state_info[v][1], label=state_info[v][0])
    for v in values
]

plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc='upper left')

plt.show()