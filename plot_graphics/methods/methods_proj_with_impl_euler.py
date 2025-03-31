import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

file_paths = [
    r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\rk4.txt",
    r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\dp8.txt",
    r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\euler.txt",
    r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\adapt_runge.txt",
    r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\new_adam_mult.txt",
    r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\pr_cor.txt",
    r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\im_euler.txt",
]
methods_names = [
    "Runge-Kutta 4",
    "Dormand-Prince 8",
    "Euler",
    "Adaptive Runge-Kutta",
    "Implicit method",
    "Predictor-Corrector",
    "Implicit Euler method"
]
colors = ['blue', 'green', 'red', 'purple', 'orange', 'brown', 'red']

fig, axes = plt.subplots(7, 3, figsize=(12, 16))
projections = [("x-y", 2, 3), ("x-z", 2, 4), ("y-z", 3, 4)]

for i, file_path in enumerate(file_paths):
    data = np.loadtxt(file_path)
    x = data[:, 2]
    y = data[:, 3]
    z = data[:, 4]

    for j, (title, col1, col2) in enumerate(projections):
        ax = axes[i, j]
        ax.plot(data[:, col1], data[:, col2], color=colors[i], linewidth=1.5)
        ax.set_title(f"{methods_names[i]}: {title}", fontsize=8)
        ax.set_xlabel(title.split("-")[0], fontsize=6)
        ax.set_ylabel(title.split("-")[1], fontsize=6)
        ax.tick_params(axis='both', which='major', labelsize=5)
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
