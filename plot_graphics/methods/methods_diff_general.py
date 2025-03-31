import numpy as np
import matplotlib.pyplot as plt
from matplotlib.scale import FuncScale

data_euler = np.loadtxt(r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\euler.txt")
data_rk4 = np.loadtxt(r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\rk4.txt")
dp8 = np.loadtxt(r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\dp8.txt")
adapt_runge = np.loadtxt(r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\adapt_runge.txt")
impl_euler = np.loadtxt(r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\new_adam_mult.txt")
pr_cor = np.loadtxt(r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\pr_cor.txt")
time = dp8[:, 1]


def compute_difference(dp8_met, other_data):
    diff_of_x = np.abs(dp8_met[:, 2] - other_data[:, 2])
    diff_of_y = np.abs(dp8_met[:, 3] - other_data[:, 3])
    diff_of_z = np.abs(dp8_met[:, 4] - other_data[:, 4])
    return diff_of_x, diff_of_y, diff_of_z


diff_x_rk4_dp8, diff_y_rk4_dp8, diff_z_rk4_dp8 = compute_difference(dp8, data_rk4)
diff_x_euler_dp8, diff_y_euler_dp8, diff_z_euler_dp8 = compute_difference(dp8, data_euler)
diff_x_adapt_dp8, diff_y_adapt_dp8, diff_z_adapt_dp8 = compute_difference(dp8, adapt_runge)
diff_x_impl_dp8, diff_y_impl_dp8, diff_z_impl_dp8 = compute_difference(dp8, impl_euler)
diff_x_pr_cor_dp8, diff_y_pr_cor_dp8, diff_z_pr_cor_dp8 = compute_difference(dp8, pr_cor)


def filter_small_values(diff, threshold=1e-15):
    return np.where(diff < threshold, np.nan, diff)


diff_x_rk4_dp8, diff_y_rk4_dp8, diff_z_rk4_dp8 = map(lambda d: filter_small_values(d),
                                                     [diff_x_rk4_dp8, diff_y_rk4_dp8, diff_z_rk4_dp8])
diff_x_euler_dp8, diff_y_euler_dp8, diff_z_euler_dp8 = map(lambda d: filter_small_values(d),
                                                           [diff_x_euler_dp8, diff_y_euler_dp8, diff_z_euler_dp8])
diff_x_adapt_dp8, diff_y_adapt_dp8, diff_z_adapt_dp8 = map(lambda d: filter_small_values(d),
                                                           [diff_x_adapt_dp8, diff_y_adapt_dp8, diff_z_adapt_dp8])
diff_x_impl_dp8, diff_y_impl_dp8, diff_z_impl_dp8 = map(lambda d: filter_small_values(d),
                                                        [diff_x_impl_dp8, diff_y_impl_dp8, diff_z_impl_dp8])
diff_x_pr_cor_dp8, diff_y_pr_cor_dp8, diff_z_pr_cor_dp8 = map(lambda d: filter_small_values(d),
                                                              [diff_x_pr_cor_dp8, diff_y_pr_cor_dp8, diff_z_pr_cor_dp8])

fig, axes = plt.subplots(5, 1, figsize=(12, 16), sharex=True)
fig.suptitle("Сравнение разниц между методами и Дормана-Принса", fontsize=16)

methods = [
    ("Метод Рунге-Кутты 4", diff_x_rk4_dp8, diff_y_rk4_dp8, diff_z_rk4_dp8),
    ("Метод Эйлера", diff_x_euler_dp8, diff_y_euler_dp8, diff_z_euler_dp8),
    ("Адаптивный метод Рунге-Кутты", diff_x_adapt_dp8, diff_y_adapt_dp8, diff_z_adapt_dp8),
    ("Неявный метод", diff_x_impl_dp8, diff_y_impl_dp8, diff_z_impl_dp8),
    ("Метод предиктора-корректора", diff_x_pr_cor_dp8, diff_y_pr_cor_dp8, diff_z_pr_cor_dp8),
]

colors = ["blue", "green", "red"]


def exp_forward(x):
    return np.exp(x)


def exp_inverse(x):
    return np.log(x)


for i, (title, diff_x, diff_y, diff_z) in enumerate(methods):
    ax = axes[i]
    ax.plot(time, diff_x, label="Разница по x", color=colors[0])
    ax.plot(time, diff_y, label="Разница по y", color=colors[1])
    ax.plot(time, diff_z, label="Разница по z", color=colors[2])
    ax.set_title(title)
    ax.set_ylabel("Абсолютная\nразница")
    ax.grid()
    ax.legend()
    min_val = np.nanmin([np.nanmin(diff_x), np.nanmin(diff_y), np.nanmin(diff_z)])
    max_val = np.nanmax([np.nanmax(diff_x), np.nanmax(diff_y), np.nanmax(diff_z)])
    scaled_diff_x = np.log(diff_x / min_val)
    scaled_diff_y = np.log(diff_y / min_val)
    scaled_diff_z = np.log(diff_z / min_val)
    ax.set_yscale(FuncScale(ax, (exp_forward, exp_inverse)))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{np.exp(x) * min_val:.1e}"))
    ax.set_ylim([min_val * 0.5, max_val * 2])

axes[-1].set_xlabel("Время")
plt.tight_layout()
plt.subplots_adjust(top=0.90, bottom=0.05)
plt.show()
