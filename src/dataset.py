"""Carregamento, processamento e visualização de conjuntos de dados experimentais."""

import csv
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from utils import interpolate_signals


def compute_abs_phase(X):
    """Calcula magnitude e fase dos elementos de um vetor complexo.

    Zera a fase dos elementos nulos para evitar uma fase indefinida.
    """
    X = np.round(X, 10)
    abs_X = np.absolute(X)
    phase_X = np.array([0 if np.absolute(z) == 0 else np.angle(z) for z in X])
    return abs_X, phase_X


class Dataset:
    """Representa dados temporais experimentais e suas operações de análise."""

    def __init__(self, file_path):
        """Inicializa um conjunto de dados associado a um arquivo.

        Args:
            file_path: Caminho do arquivo TXT, CSV ou NPZ de entrada.
        """
        # Identificação e localização do conjunto de dados.
        self.file_path = file_path  # Caminho do arquivo de dados.
        self.file_name = Path(file_path).stem  # Nome-base usado nos arquivos de saída.

        # Sinais carregados: tempo, correntes, perturbações e saídas.
        self.t = None  # Vetor de tempo.
        self.i = None  # Matriz de correntes de entrada.
        self.d = None  # Matriz de sinais de perturbação.
        self.y = None  # Matriz de sinais de saída.

    def load(self):
        """Carrega os dados usando o leitor correspondente à extensão."""
        suffix = Path(self.file_path).suffix.lower()

        if suffix == ".npz":
            with np.load(self.file_path) as data:
                self.t = data["t"]
                self.i = data["i"]
                self.d = data["d"]
                self.y = data["y"]
            return self.t, self.i, self.d, self.y

        if suffix == ".txt":
            return self.load_txt()
        if suffix == ".csv":
            return self.load_csv()

        raise ValueError(
            f"Extensão de arquivo não suportada: '{Path(self.file_path).suffix}'"
        )

    def load_txt(self):
        """Carrega dados TXT ou seu cache NPZ associado e retorna os sinais."""
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
        self.i = data[:, 1:5].T
        self.d = data[:, 5:9].T
        self.y = data[:, 9:].T

        np.savez(npz_path, t=self.t, i=self.i, d=self.d, y=self.y)
        return self.t, self.i, self.d, self.y

    def load_csv(self):
        """Carrega dados CSV ou seu cache NPZ associado e retorna os sinais."""
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
        self.i = data[:, 1:5].T
        self.d = data[:, [5, 7, 6, 8]].T
        self.y = data[:, [9, 11, 10, 12]].T

        np.savez(npz_path, t=self.t, i=self.i, d=self.d, y=self.y)
        return self.t, self.i, self.d, self.y

    @staticmethod
    def _decimate_signal(signal, decimation_factor):
        """Reduz a taxa de amostragem de um sinal por subamostragem regular."""
        if not isinstance(decimation_factor, (int, np.integer)):
            raise TypeError("O fator de decimação deve ser um número inteiro.")
        if decimation_factor < 1:
            raise ValueError("O fator de decimação deve ser maior ou igual a 1.")

        signal = np.asarray(signal, dtype=float)
        return signal[::decimation_factor]

    def compute_crosscorrelation_function(self, u, y, decimation_factor=1_000):
        """Calcula a correlação cruzada normalizada entre dois sinais.

        Retorna os atrasos em amostras, a correlação e o intervalo de confiança
        aproximado de 95%.
        """
        u = self._decimate_signal(u, decimation_factor)
        y = self._decimate_signal(y, decimation_factor)
        if len(u) != len(y):
            raise ValueError("Os sinais devem ter o mesmo número de amostras.")

        u = u.copy()
        y = y.copy()
        mean_u = np.mean(u)
        sigma_u = np.std(u)
        mean_y = np.mean(y)
        sigma_y = np.std(y)
        n_samples = len(u)

        u -= mean_u
        u = u / sigma_u if sigma_u != 0 else u
        y -= mean_y
        y = y / sigma_y if sigma_y != 0 else y

        r_uy = np.correlate(y, u, mode="full") / n_samples
        confidence_interval = 1.96 / np.sqrt(n_samples)
        k_values = np.arange(-n_samples + 1, n_samples)

        return k_values, r_uy, confidence_interval

    def compute_autocorrelation_function(self, u, decimation_factor=1_000):
        """Calcula a autocorrelação normalizada para atrasos não negativos."""
        k_values, r_uu, confidence_interval = self.compute_crosscorrelation_function(
            u, u, decimation_factor=decimation_factor
        )
        nonnegative_lags = k_values >= 0
        return (
            k_values[nonnegative_lags],
            r_uu[nonnegative_lags],
            confidence_interval,
        )

    def plot_crosscorrelation(self, decimation_factor=500):
        """Gera e salva as funções de correlação entre entradas e saídas."""
        n_inputs = self.i.shape[0]
        n_outputs = self.y.shape[0]
        sampling_period = np.mean(np.diff(self.t))
        figure, axes = plt.subplots(
            n_inputs,
            n_outputs,
            figsize=(4 * n_outputs, 2.5 * n_inputs),
            squeeze=False,
            sharex=True,
        )

        for i_u in range(n_inputs):
            for j_y in range(n_outputs):
                k_values, r_uy, confidence_interval = (
                    self.compute_crosscorrelation_function(
                        self.i[i_u, :],
                        self.y[j_y, :],
                        decimation_factor=decimation_factor,
                    )
                )
                axis = axes[i_u, j_y]
                lag_times = k_values * sampling_period * decimation_factor
                axis.plot(lag_times, r_uy)
                axis.axhline(confidence_interval, color="black", linestyle="--")
                axis.axhline(-confidence_interval, color="black", linestyle="--")
                axis.set_title(f"$u_{i_u + 1}$ x $y_{j_y + 1}$")
                axis.set_ylabel("Correlação")
                axis.grid()

        for axis in axes[-1, :]:
            axis.set_xlabel("Lag (s)")

        figure.tight_layout()
        plots_dir = Path(__file__).resolve().parent.parent / "plots"
        plots_dir.mkdir(exist_ok=True)
        figure.savefig(
            plots_dir / (self.file_name + "_crosscorrelation.svg"),
            format="svg",
        )
        return figure, axes

    def plot_autocorrelation(self, decimation_factor=500):
        """Gera e salva as funções de autocorrelação das entradas."""
        n_inputs = self.i.shape[0]
        sampling_period = np.mean(np.diff(self.t))
        figure, axes = plt.subplots(
            n_inputs,
            1,
            figsize=(10, 2.5 * n_inputs),
            squeeze=False,
            sharex=True,
        )

        for i_u in range(n_inputs):
            k_values, r_uu, confidence_interval = self.compute_autocorrelation_function(
                self.i[i_u, :], decimation_factor=decimation_factor
            )
            axis = axes[i_u, 0]
            lag_times = k_values * sampling_period * decimation_factor
            axis.plot(lag_times, r_uu)
            axis.axhline(confidence_interval, color="black", linestyle="--")
            axis.axhline(-confidence_interval, color="black", linestyle="--")
            axis.set_title(f"Autocorrelação de $u_{i_u + 1}$")
            axis.set_ylabel("Autocorrelação")
            axis.grid()

        axes[-1, 0].set_xlabel("Lag (s)")
        figure.tight_layout()
        plots_dir = Path(__file__).resolve().parent.parent / "plots"
        plots_dir.mkdir(exist_ok=True)
        figure.savefig(
            plots_dir / (self.file_name + "_autocorrelation.svg"), format="svg"
        )
        return figure, axes

    def compute_fft(self):
        """Calcula e salva os espectros de magnitude das perturbações e saídas."""
        if self.t is None or self.d is None or self.y is None:
            raise ValueError("Os dados devem ser carregados antes de computar a FFT.")
        if len(self.t) < 2:
            raise ValueError("O vetor de tempo deve conter pelo menos duas amostras.")
        if self.d.ndim != 2 or self.d.shape[1] != len(self.t):
            raise ValueError(
                "A matriz de perturbações deve ter uma amostra por coluna "
                "correspondente ao vetor de tempo."
            )
        if self.y.ndim != 2 or self.y.shape[1] != len(self.t):
            raise ValueError(
                "A matriz de saídas deve ter uma amostra por coluna "
                "correspondente ao vetor de tempo."
            )

        sampling_period = self.t[1] - self.t[0]
        if sampling_period == 0:
            raise ValueError("O período de amostragem deve ser diferente de zero.")

        n_samples = self.d.shape[1]
        frequencies = np.fft.fftfreq(n_samples, sampling_period)
        positive_frequencies = frequencies >= 0

        plots_dir = Path(__file__).resolve().parent.parent / "plots"
        plots_dir.mkdir(exist_ok=True)

        def plot_fft(signals, signal_name, title, file_suffix):
            """Plota e salva as magnitudes FFT dos canais fornecidos."""
            n_signals = signals.shape[0]
            figure, axes = plt.subplots(
                n_signals,
                1,
                figsize=(10, 2.5 * n_signals),
                sharex=True,
                squeeze=False,
            )

            for channel_index, signal in enumerate(signals):
                spectrum = np.fft.fft(signal) * 2 / n_samples
                magnitude, _ = compute_abs_phase(spectrum)
                axis = axes[channel_index, 0]
                axis.plot(
                    frequencies[positive_frequencies], magnitude[positive_frequencies]
                )
                axis.set_ylabel(rf"$|{signal_name}_{channel_index + 1}|$")
                axis.grid()

            axes[-1, 0].set_xlabel("Frequência (Hz)")
            figure.suptitle(title)
            figure.tight_layout()
            figure.savefig(plots_dir / (self.file_name + file_suffix), format="svg")

        plot_fft(
            self.d,
            "D",
            "Magnitude das FFTs das perturbações",
            "_disturbance_fft.svg",
        )
        plot_fft(
            self.y,
            "Y",
            "Magnitude das FFTs das saídas",
            "_output_fft.svg",
        )

    def plot(self):
        """Interpola e salva gráficos temporais de correntes, perturbações e saídas."""
        plots_dir = Path(__file__).resolve().parent.parent / "plots"
        plots_dir.mkdir(exist_ok=True)

        t_plot, i_plot = interpolate_signals(self.t, self.i)
        _, d_plot = interpolate_signals(self.t, self.d)
        _, y_plot = interpolate_signals(self.t, self.y)

        legend_loc = "upper left"
        plt.figure(figsize=(8, 4), dpi=250)
        plt.plot(t_plot, np.transpose(i_plot))
        plt.title("Current")
        plt.xlabel("Time")
        plt.ylabel("Amplitude")
        plt.legend(["$i_1$", "$i_2$", "$i_3$", "$i_4$"], loc=legend_loc)
        plt.tight_layout()
        plt.grid()
        plt.savefig(plots_dir / (self.file_name + ".svg"), format="svg")

        plt.figure(figsize=(8, 4), dpi=250)
        plt.plot(t_plot, np.transpose(d_plot))
        plt.title("Disturbance")
        plt.xlabel("Time")
        plt.ylabel("Amplitude")
        plt.legend(["$d_1$", "$d_2$", "$d_3$", "$d_4$"], loc=legend_loc)
        plt.tight_layout()
        plt.grid()
        plt.savefig(plots_dir / (self.file_name + "_disturbance.svg"), format="svg")

        plt.figure(figsize=(8, 4), dpi=250)
        plt.plot(t_plot, np.transpose(y_plot))
        plt.title("Output")
        plt.xlabel("Time")
        plt.ylabel("Amplitude")
        plt.legend(["$y_1$", "$y_2$", "$y_3$", "$y_4$"], loc=legend_loc)
        plt.tight_layout()
        plt.grid()
        plt.savefig(plots_dir / (self.file_name + "_output.svg"), format="svg")


if __name__ == "__main__":
    load_data = Dataset(file_path="../data/chirp_v13_0.csv")
    load_data.load()
    load_data.plot()
    load_data.plot_crosscorrelation()
    load_data.plot_autocorrelation()
    load_data.compute_fft()
    plt.show()
