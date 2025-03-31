import matplotlib.pyplot as plt
import pandas as pd

file_paths = [
    r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\dp8_001.txt",
    r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\dp8_002.txt",
    r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\dp8_003.txt",
    r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\dp8_004.txt",
    r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\dp8_005.txt",
]
data_paths = [
    r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\dp8_006.txt",
    r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\dp8_007.txt",
    r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\dp8_008.txt",
    r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\dp8_009.txt",
    r"C:\Users\Huawei\source\repos\ConsoleCpp\ConsoleCpp\dp8_01.txt"
]


def plot_all_files(file):
    fig, axes = plt.subplots(len(file), 3, figsize=(18, 10 * len(file)))

    for i, filename in enumerate(file):
        column_names = ['счетчик', 'время', 'x', 'y', 'z']
        try:
            data = pd.read_csv(filename, header=None, names=column_names, sep="\s+")
            print(f"Данные из файла {filename}:")
            print(data.head())
        except Exception as e:
            print(f"Ошибка при чтении файла {filename}: {e}")
            continue

        s = data['счетчик']
        x = data['x']
        y = data['y']
        z = data['z']

        axes[i][0].plot(s, x)
        axes[i][0].set_title(f"x(s), шаг=0.{filename.split('_')[-1].split('.')[0]}")
        axes[i][0].set_xlabel("Количество итераций (s)")
        axes[i][0].set_ylabel("x")
        axes[i][0].grid()

        axes[i][1].plot(s, y)
        axes[i][1].set_title(f"y(s), шаг=0.{filename.split('_')[-1].split('.')[0]}")
        axes[i][1].set_xlabel("Количество итераций (s)")
        axes[i][1].set_ylabel("y")
        axes[i][1].grid()

        axes[i][2].plot(s, z)
        axes[i][2].set_title(f"z(s), шаг=0.{filename.split('_')[-1].split('.')[0]}")
        axes[i][2].set_xlabel("Количество итераций (s)")
        axes[i][2].set_ylabel("z")
        axes[i][2].grid()

    plt.tight_layout(h_pad=3.0)
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, hspace=0.8, wspace=0.4)
    plt.show()


plot_all_files(file_paths)
plot_all_files(data_paths)
