import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data_files = {
    "rk4": r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\rk4.txt",
    "dp8": r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\dp8.txt",
    "euler": r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\euler.txt",
    "adapt": r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\adapt_runge.txt",
    "impl": r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\new_adam_mult.txt",
    "impl_euler": r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\im_euler.txt",
    "predictor_cor": r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\pr_cor.txt"
}

data = {}
for method, file_path in data_files.items():
    data[method] = pd.read_csv(file_path, sep=" ", header=None,
                               names=["Шаг", "dt", f"x_{method}", f"y_{method}", f"z_{method}"])


def check_data(data_dict):
    for method, df in data_dict.items():
        if df.isnull().values.any():
            print(f"Ошибка: Данные {method} содержат NaN")
            exit()
        if np.isinf(df.values).any():
            print(f"Ошибка: Данные {method} содержат Inf")
            exit()


check_data(data)


def process_differences(base_method, compare_methods, output_prefix):
    base_data = data[base_method]
    for method in compare_methods:
        compare_data = data[method]

        diff_x = base_data[f"x_{base_method}"] - compare_data[f"x_{method}"]
        diff_y = base_data[f"y_{base_method}"] - compare_data[f"y_{method}"]
        diff_z = base_data[f"z_{base_method}"] - compare_data[f"z_{method}"]

        diff_data = pd.DataFrame({
            "Шаг": base_data["Шаг"],
            "dt": base_data["dt"],
            "Diff_x": diff_x,
            "Diff_y": diff_y,
            "Diff_z": diff_z
        })

        if diff_data.isnull().values.any():
            print(f"Ошибка: DataFrame для {method} содержит NaN")
            continue
        if diff_data.empty:
            print(f"Ошибка: DataFrame для {method} пустой")
            continue

        try:
            diff_data.to_csv(f"{output_prefix}_{method}.txt", sep=" ", index=False)
            print(f"Файл {output_prefix}_{method}.txt успешно создан")
        except Exception as e:
            print(f"Ошибка при записи файла {output_prefix}_{method}.txt: {e}")

        plt.figure(figsize=(15, 10), constrained_layout=True)
        titles = [f"Разница x ({base_method.upper()} - {method.upper()})",
                  f"Разница y ({base_method.upper()} - {method.upper()})",
                  f"Разница z ({base_method.upper()} - {method.upper()})"]
        colors = ['red', 'blue', 'green']

        for i, (diff, title, color) in enumerate(zip([diff_x, diff_y, diff_z], titles, colors)):
            plt.subplot(3, 1, i + 1)
            plt.plot(base_data['dt'], diff, color=color, linewidth=1)
            plt.title(title, fontsize=12)
            plt.xlabel('Время', fontsize=10)
            plt.ylabel('Разница', fontsize=10)
            plt.grid(True, linestyle='--', alpha=0.5)
            plt.ticklabel_format(axis='y', style='scientific', scilimits=(0, 0))

        plt.savefig(f'difference_{method}_{base_method}.png', dpi=300)
        plt.show()
        plt.close()

        # Вывод средних значений разниц
        print(f"\nСредние значения разниц для {method}:")
        print(f"Среднее значение разницы x: {diff_x.mean()}")
        print(f"Среднее значение разницы y: {diff_y.mean()}")
        print(f"Среднее значение разницы z: {diff_z.mean()}")


process_differences("dp8", ["rk4", "euler", "adapt", "impl", "predictor_cor"], "differences")
# process_differences("rk4", ['dp8'], "diff_rk4")
# process_differences("dp8", ['impl_euler'], "diff_dp8")
