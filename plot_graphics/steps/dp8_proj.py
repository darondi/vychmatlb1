import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

file_paths = [
    r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\dp8_001.txt",
    r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\dp8_002.txt",
    r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\dp8_003.txt",
    r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\dp8_004.txt",
    r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\dp8_005.txt"
]
data_paths = [
    r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\dp8_006.txt",
    r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\dp8_007.txt",
    r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\dp8_008.txt",
    r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\dp8_009.txt",
    r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\dp8_01.txt"
]
method_names = [
    "0.001", "0.002", "0.003", "0.004", "0.005"
]
met_names = [
    "0.006", "0.007", "0.008", "0.009", "0.01"
]
colors = [
    'blue', 'green', 'red', 'purple', 'orange'
]
projections = [("x-y", 2, 3), ("x-z", 2, 4), ("y-z", 3, 4)]

fig1, axes1 = plt.subplots(5, 3, figsize=(12, 16))

for i, file_path in enumerate(file_paths):
    data = np.loadtxt(file_path)
    x = data[:, 2]
    y = data[:, 3]
    z = data[:, 4]

    for j, (title, col1, col2) in enumerate(projections):
        ax = axes1[i, j]
        ax.plot(data[:, col1], data[:, col2], color=colors[i], linewidth=1.5)
        ax.set_title(f"{method_names[i]}: {title}", fontsize=8)
        ax.set_xlabel(title.split("-")[0], fontsize=6)
        ax.set_ylabel(title.split("-")[1], fontsize=6)
        ax.tick_params(axis='both', which='major', labelsize=6)
        ax.grid(True)
        x_min, x_max = np.min(data[:, col1]), np.max(data[:, col1])
        y_min, y_max = np.min(data[:, col2]), np.max(data[:, col2])
        x_step = 0.2
        y_step = 0.3
        ax.xaxis.set_major_locator(MultipleLocator(x_step))
        ax.yaxis.set_major_locator(MultipleLocator(y_step))
        ax.xaxis.set_minor_locator(MultipleLocator(x_step / 2))
        ax.yaxis.set_minor_locator(MultipleLocator(y_step / 2))

plt.tight_layout(pad=2)
plt.subplots_adjust(wspace=0.6, hspace=0.7)

fig2, axes2 = plt.subplots(5, 3, figsize=(12, 16))

for i, data_path in enumerate(data_paths):
    data = np.loadtxt(data_path)
    x = data[:, 2]
    y = data[:, 3]
    z = data[:, 4]

    for j, (title, col1, col2) in enumerate(projections):
        ax = axes2[i, j]
        ax.plot(data[:, col1], data[:, col2], color=colors[i], linewidth=1.5)
        ax.set_title(f"{met_names[i]}: {title}", fontsize=8)
        ax.set_xlabel(title.split("-")[0], fontsize=6)
        ax.set_ylabel(title.split("-")[1], fontsize=6)
        ax.tick_params(axis='both', which='major', labelsize=6)
        ax.grid(True)
        x_min, x_max = np.min(data[:, col1]), np.max(data[:, col1])
        y_min, y_max = np.min(data[:, col2]), np.max(data[:, col2])
        x_step = 0.2
        y_step = 0.3
        ax.xaxis.set_major_locator(MultipleLocator(x_step))
        ax.yaxis.set_major_locator(MultipleLocator(y_step))
        ax.xaxis.set_minor_locator(MultipleLocator(x_step / 2))
        ax.yaxis.set_minor_locator(MultipleLocator(y_step / 2))

plt.tight_layout(pad=2)
plt.subplots_adjust(wspace=0.6, hspace=0.7)
plt.show()
