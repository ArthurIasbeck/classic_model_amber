import csv
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


class Dataset:
    def __init__(self, file_path):
        self.file_path = file_path
        self.file_name = Path(file_path).stem
        self.t = None
        self.i = None
        self.d = None
        self.y = None

    def load(self):
        suffix = Path(self.file_path).suffix.lower()

        if suffix == ".txt":
            return self.load_txt()
        if suffix == ".csv":
            return self.load_csv()

        raise ValueError(
            f"Extensão de arquivo não suportada: '{Path(self.file_path).suffix}'"
        )

    def load_txt(self):
        requested_path = Path(self.file_path)
        base_path = (
            requested_path.with_suffix("") if requested_path.suffix else requested_path
        )
        npz_path = base_path.with_suffix(".npz")
        txt_path = base_path.with_suffix(".txt")

        if npz_path.exists():
            with np.load(npz_path) as data:
                self.t = data["t"]
                self.i = data["i"]
                self.d = data["d"]
                self.y = data["y"]
            return self.t, self.i, self.d, self.y

        if not txt_path.exists():
            raise FileNotFoundError(
                f"Nenhum arquivo encontrado: '{npz_path}' ou '{txt_path}'"
            )

        data = np.loadtxt(txt_path, delimiter=",")
        self.t = data[:, 0]
        self.i = data[:, 1:5]
        self.d = data[:, 5:9]
        self.y = data[:, 9:]

        np.savez(npz_path, t=self.t, i=self.i, d=self.d, y=self.y)
        return self.t, self.i, self.d, self.y

    def load_csv(self):
        requested_path = Path(self.file_path)
        base_path = (
            requested_path.with_suffix("") if requested_path.suffix else requested_path
        )
        npz_path = base_path.with_suffix(".npz")
        csv_path = base_path.with_suffix(".csv")

        if npz_path.exists():
            with np.load(npz_path) as data:
                self.t = data["t"]
                self.i = data["i"]
                self.d = data["d"]
                self.y = data["y"]
            return self.t, self.i, self.d, self.y

        if not csv_path.exists():
            raise FileNotFoundError(
                f"Nenhum arquivo encontrado: '{npz_path}' ou '{csv_path}'"
            )

        with csv_path.open("r", encoding="cp1252", newline="") as file:
            reader = csv.reader(file)
            trace_values_start = None

            for row in reader:
                if row and row[0] == "trace_values":
                    trace_values_start = row
                    break

            if trace_values_start is None:
                raise ValueError("Linha 'trace_values' não encontrada no arquivo.")

            data = [[float(value) for value in trace_values_start[1:]]]

            for row in reader:
                if row:
                    data.append([float(value) for value in row[1:]])

        data = np.asarray(data, dtype=np.float64)
        if data.ndim != 2 or data.shape[1] != 13:
            raise ValueError(
                "O arquivo CSV deve conter uma coluna de tempo e 12 colunas de sinais."
            )

        self.t = data[:, 0]
        self.i = data[:, 1:5]
        self.d = data[:, [5, 7, 6, 8]]
        self.y = data[:, [9, 11, 10, 12]]

        np.savez(npz_path, t=self.t, i=self.i, d=self.d, y=self.y)
        return self.t, self.i, self.d, self.y

    def plot(self):
        plots_dir = Path(__file__).resolve().parent.parent / "plots"
        plots_dir.mkdir(exist_ok=True)

        n_plot_points = 1_000_000
        t_plot = np.linspace(self.t[0], self.t[-1], n_plot_points)

        def interpolate(data):
            return np.column_stack(
                [
                    np.interp(t_plot, self.t, data[:, column])
                    for column in range(data.shape[1])
                ]
            )

        i_plot = interpolate(self.i)
        d_plot = interpolate(self.d)
        y_plot = interpolate(self.y)

        plt.figure(figsize=(8, 4), dpi=250)
        plt.plot(t_plot, i_plot)
        plt.title("Current")
        plt.xlabel("Time")
        plt.ylabel("Amplitude")
        plt.legend(["$i_1$", "$i_2$", "$i_3$", "$i_4$"], loc="upper left")
        plt.tight_layout()
        plt.grid()
        plt.savefig(plots_dir / (self.file_name + ".svg"), format="svg")

        plt.figure(figsize=(8, 4), dpi=250)
        plt.plot(t_plot, d_plot)
        plt.title("Disturbance")
        plt.xlabel("Time")
        plt.ylabel("Amplitude")
        plt.legend(["$d_1$", "$d_2$", "$d_3$", "$d_4$"], loc="upper left")
        plt.tight_layout()
        plt.grid()
        plt.savefig(plots_dir / (self.file_name + "_disturbance.svg"), format="svg")

        plt.figure(figsize=(8, 4), dpi=250)
        plt.plot(t_plot, y_plot)
        plt.title("Output")
        plt.xlabel("Time")
        plt.ylabel("Amplitude")
        plt.legend(["$y_1$", "$y_2$", "$y_3$", "$y_4$"], loc="upper left")
        plt.tight_layout()
        plt.grid()
        plt.savefig(plots_dir / (self.file_name + "_output.svg"), format="svg")


if __name__ == "__main__":
    load_data = Dataset(file_path="../data/random_v13.csv")
    load_data.load()
    load_data.plot()
    plt.show()
