from pathlib import Path

import control
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter
from scipy.signal import butter, sosfiltfilt

from utils import interpolate_signals

plt.rcParams.update(
    {
        "font.family": "TeX Gyre Termes",
        "font.size": 14,
        "axes.labelsize": 16,
        "axes.titlesize": 16,
        "xtick.labelsize": 13,
        "ytick.labelsize": 13,
        "legend.fontsize": 13,
        "figure.titlesize": 18,
        "svg.fonttype": "none",
    }
)


def _format_decimal_comma(value, _):
    return f"{value:g}".replace(".", ",")


class ExperimentalFrequencyResponseSiso:
    def __init__(
        self,
        data_name,
        t_experiments,
        u_experiments,
        y_experiments,
        data_dir="../data",
        plots_dir="../plots",
    ):

        self.data_name = data_name
        self.t = list(t_experiments)
        self.u = list(u_experiments)
        self.y = list(y_experiments)
        self.y_original = None

        self.n_experiments = len(self.t)
        self.angular_frequencies = None
        self.response = None
        self.plots_dir = plots_dir
        self.data_dir = data_dir

        self.clip_threshold = 0.01
        self.exp_duration = 601

    def clip_data(self, plot=True):
        if plot:
            plt.figure(figsize=(6, 5), dpi=180)

        diffs = []

        for i in range(self.n_experiments):
            diff = np.diff(self.t[i])
            if plot:
                plt.plot(diff)
            diffs.append(diff)

        dt_mean = np.mean(np.concatenate(diffs))
        if plot:
            plt.axhline(y=dt_mean, color="k", linestyle="--")

        if plot:
            plt.xlabel("Amostra")
            plt.ylabel("Diferença de tempo")
            plt.title("Diferença entre passos de tempo")
            plt.gca().xaxis.set_major_formatter(FuncFormatter(_format_decimal_comma))
            plt.gca().yaxis.set_major_formatter(FuncFormatter(_format_decimal_comma))
            plt.grid(True)
            plt.tight_layout()
            plt.savefig(
                Path(self.plots_dir) / "mean_dt.svg",
                format="svg",
            )

        for i in range(self.n_experiments):
            t = self.t[i]
            u = self.u[i]
            y = self.y[i]

            t_interp = np.arange(t[0], t[-1], dt_mean)
            u_interp = np.interp(t_interp, t, u)
            y_interp = np.interp(t_interp, t, y)

            t = t_interp
            u = u_interp
            y = y_interp

            dt = t[1] - t[0]
            i_0 = np.argmax(u > self.clip_threshold) - np.round(1 / dt).astype(int)
            i_f = i_0 + np.round(self.exp_duration / dt).astype(int)

            self.t[i] = t[i_0:i_f] - t[i_0]
            self.u[i] = u[i_0:i_f]
            self.y[i] = y[i_0:i_f]

            print("\n===========================================")
            print(f"\nExperimento {i + 1}")
            print(f"dt = {dt} s")
            print(f"i_0 = {i_0} | i_f = {i_f}")
            print(f"t = {self.t[i].shape}")

        print("\n===========================================\n")

    def filter_output(self, cutoff_frequency, plot=False):
        filtered_y_experiments = []

        for t, y in zip(self.t, self.y):
            sampling_frequency = 1 / (t[1] - t[0])
            nyquist_frequency = sampling_frequency / 2
            if cutoff_frequency >= nyquist_frequency:
                raise ValueError(
                    "the cutoff frequency must be below the Nyquist frequency"
                )

            second_order_filter = butter(
                4,
                cutoff_frequency,
                btype="lowpass",
                fs=sampling_frequency,
                output="sos",
            )

            filtered_y = sosfiltfilt(
                second_order_filter, np.asarray(y, dtype=float).copy()
            )

            filtered_y_experiments.append(filtered_y)

        self.y_original = self.y
        self.y = filtered_y_experiments

        if plot:
            self.plot_filter_output()

        return self.y

    def plot_filter_output(self):
        if self.y_original is None:
            raise ValueError(
                "filter_output must be called before plotting filtered data"
            )

        plots_directory = Path(self.plots_dir)
        plots_directory.mkdir(parents=True, exist_ok=True)

        for experiment_index, (t, original_y, filtered_y) in enumerate(
            zip(self.t, self.y_original, self.y), start=1
        ):
            t_plot, original_y_plot = interpolate_signals(t, np.asarray(original_y))
            _, filtered_y_plot = interpolate_signals(t, np.asarray(filtered_y))
            figure, axis = plt.subplots(figsize=(8, 6), dpi=200)
            axis.plot(
                t_plot,
                original_y_plot,
                color="C0",
                alpha=0.8,
                linestyle="-",
                label="Original",
            )

            axis.plot(
                t_plot,
                filtered_y_plot,
                color="C3",
                linestyle="-",
                label="Filtrado",
            )
            axis.set_ylabel(r"$y$ [$\mu$m]")
            axis.set_xlabel("Tempo (s)")
            axis.xaxis.set_major_formatter(FuncFormatter(_format_decimal_comma))
            axis.yaxis.set_major_formatter(FuncFormatter(_format_decimal_comma))
            axis.grid(True)
            axis.legend(loc="upper right")
            figure.tight_layout()
            figure.savefig(
                plots_directory / f"experimental_filtered_data_{experiment_index}.svg",
                format="svg",
            )

    def compute(self):
        self.response = None
        self.angular_frequencies = None

        s_uu_values = []
        s_yu_values = []
        frequency = None

        for experiment_index, (t, u, y) in enumerate(
            zip(self.t, self.u, self.y), start=1
        ):
            dt = t[1] - t[0]
            f_s = 1 / dt
            N = len(u)
            U = np.fft.fft(u)
            Y = np.fft.fft(y)
            S_uu = 1 / N * U * np.conj(U)
            S_yu = 1 / N * Y * np.conj(U)
            frequency = np.fft.fftfreq(N, d=1 / f_s)

            s_uu_values.append(S_uu)
            s_yu_values.append(S_yu)

        S_uu_mean = np.mean(s_uu_values, axis=0)
        S_yu_mean = np.mean(s_yu_values, axis=0)

        self.response = np.full_like(S_yu_mean, np.nan + 1j * np.nan)
        valid_bins = np.abs(S_uu_mean) > np.finfo(float).eps
        self.response[valid_bins] = S_yu_mean[valid_bins] / S_uu_mean[valid_bins]
        self.angular_frequencies = 2 * np.pi * frequency

        np.savez(
            Path(self.data_dir) / f"{self.data_name}_frequency_response.npz",
            angular_frequencies=self.angular_frequencies,
            response=self.response,
        )

        return self.angular_frequencies, self.response

    def plot_exp_data(self):
        plots_directory = Path(self.plots_dir)
        plots_directory.mkdir(parents=True, exist_ok=True)

        for experiment_index, (t, u, y) in enumerate(
            zip(self.t, self.u, self.y), start=1
        ):
            t_plot, u_plot = interpolate_signals(t, np.asarray(u))
            _, y_plot = interpolate_signals(t, np.asarray(y))
            figure, axis = plt.subplots(figsize=(8, 6))
            output_axis = axis.twinx()

            input_line = axis.plot(t_plot, u_plot, label=r"Perturbação ($\mu$m)")
            output_line = output_axis.plot(
                t_plot, y_plot, color="C1", label=r"Deslocamento ($\mu$m)"
            )

            axis.set_xlabel("Tempo (s)")
            axis.set_ylabel(r"Amplitude da perturbação ($\mu$m)")
            output_axis.set_ylabel(r"Amplitude do deslocamento ($\mu$m)")
            axis.xaxis.set_major_formatter(FuncFormatter(_format_decimal_comma))
            axis.yaxis.set_major_formatter(FuncFormatter(_format_decimal_comma))
            output_axis.yaxis.set_major_formatter(FuncFormatter(_format_decimal_comma))
            axis.grid(True)
            axis.legend(
                input_line + output_line,
                [line.get_label() for line in input_line + output_line],
                loc="upper right",
            )

            figure.suptitle(f"Experimento {experiment_index}")
            figure.tight_layout()
            figure.savefig(
                plots_directory / f"experimental_data_{experiment_index}.svg",
                format="svg",
            )

    def plot_freq_resp(self, min_freq=None, max_freq=None):
        min_freq = 0 if min_freq is None else min_freq
        max_freq = np.inf if max_freq is None else max_freq

        if self.response is None or self.angular_frequencies is None:
            raise ValueError("compute must be called before plotting the response")

        fig, axes = plt.subplots(
            2,
            1,
            figsize=(10, 6),
            sharex=True,
            squeeze=False,
            dpi=200,
        )

        angular_frequency = self.angular_frequencies
        positive_frequency = (
            (angular_frequency > min_freq)
            & (angular_frequency < max_freq)
            & np.isfinite(self.response)
        )
        angular_frequency = angular_frequency[positive_frequency]
        response = self.response[positive_frequency]

        magnitude = 20 * np.log10(np.abs(response))
        phase = np.angle(response)
        angular_frequency_plot, magnitude_plot = interpolate_signals(
            angular_frequency, magnitude
        )
        _, phase_plot = interpolate_signals(angular_frequency, phase)

        axes[0, 0].semilogx(angular_frequency_plot, magnitude_plot)
        axes[0, 0].set_ylabel("Magnitude (dB)")
        axes[0, 0].xaxis.set_major_formatter(FuncFormatter(_format_decimal_comma))
        axes[0, 0].yaxis.set_major_formatter(FuncFormatter(_format_decimal_comma))
        axes[0, 0].grid(True, which="both")

        axes[1, 0].semilogx(angular_frequency_plot, phase_plot * 180 / np.pi)
        axes[1, 0].set_ylabel("Fase (graus)")
        axes[1, 0].xaxis.set_major_formatter(FuncFormatter(_format_decimal_comma))
        axes[1, 0].yaxis.set_major_formatter(FuncFormatter(_format_decimal_comma))
        axes[1, 0].grid(True, which="both")

        fig.suptitle("Magnitude da resposta em frequência experimental")
        fig.tight_layout()

        plots_directory = Path(self.plots_dir)
        plots_directory.mkdir(parents=True, exist_ok=True)
        fig.savefig(
            plots_directory / "experimental_frequency_response.svg",
            format="svg",
        )


def main():
    # Produção de sinais sintéticos
    sample_rate = 1000
    duration = 20
    n_samples = int(sample_rate * duration)
    experiment_time = np.arange(n_samples) / sample_rate
    t = np.arange(n_samples) / sample_rate

    s = control.tf("s")
    natural_frequency = 10
    damping_ratio = 0.1
    transfer_function = natural_frequency**2 / (
        s**2 + 2 * damping_ratio * natural_frequency * s + natural_frequency**2
    )
    system = control.ss(transfer_function)

    initial_frequency = 0.5
    final_frequency = 600
    sweep_rate = (final_frequency - initial_frequency) / duration
    phase = (
        2
        * np.pi
        * (initial_frequency * experiment_time + sweep_rate * experiment_time**2 / 2)
    )
    chirp_signal = np.sin(phase)

    rng = np.random.default_rng(0)
    t_experiments = []
    u_experiments = []
    y_experiments = []
    for _ in range(4):
        u = chirp_signal
        y = np.asarray(control.forced_response(system, T=t, U=u).outputs)
        y += rng.normal(0.0, 0.05, size=y.shape)
        t_experiments.append(t)
        u_experiments.append(u)
        y_experiments.append(y)

    # Emprego da classe ExperimentalFrequencyResponseSiso
    experimental_frequency_response = ExperimentalFrequencyResponseSiso(
        "synthetic", t_experiments, u_experiments, y_experiments
    )
    experimental_frequency_response.plot_exp_data()
    experimental_frequency_response.compute()
    experimental_frequency_response.plot_freq_resp(max_freq=600)


if __name__ == "__main__":
    main()
    plt.show()
