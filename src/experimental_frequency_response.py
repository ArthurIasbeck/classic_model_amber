from pathlib import Path

import control
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, sosfiltfilt

from utils import interpolate_signals


class ExperimentalFrequencyResponse:
    def __init__(self, t_experiments, u_experiments, y_experiments):

        if not (len(t_experiments) == len(u_experiments) == len(y_experiments)):
            raise ValueError("t, u, and y must contain the same number of experiments")
        if not t_experiments:
            raise ValueError("at least one experiment is required")

        self.t = list(t_experiments)
        self.u = list(u_experiments)
        self.y = list(y_experiments)
        self.y_original = None

        self.n_experiments = len(self.t)
        self.angular_frequencies = None
        self.response = None
        self.plots_dir = "../plots"

        self.legend_loc = "upper right"

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
            filtered_y = np.asarray(y, dtype=float).copy()
            for output_index in range(filtered_y.shape[0]):
                filtered_y[output_index, :] = sosfiltfilt(
                    second_order_filter, filtered_y[output_index, :]
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
            t_plot, original_y_plot = interpolate_signals(t, original_y)
            _, filtered_y_plot = interpolate_signals(t, filtered_y)
            n_y = original_y.shape[0]
            figure, axes = plt.subplots(
                n_y,
                1,
                figsize=(12, 3 * n_y),
                sharex=True,
                squeeze=False,
                dpi=100,
            )

            for output_index in range(n_y):
                axis = axes[output_index, 0]
                axis.plot(
                    t_plot,
                    original_y_plot[output_index, :],
                    color="C0",
                    alpha=0.8,
                    linestyle="-",
                    label="Original",
                )
                axis.plot(
                    t_plot,
                    filtered_y_plot[output_index, :],
                    color="C3",
                    linestyle="-",
                    label="Filtered",
                )
                axis.set_ylabel(rf"$y_{output_index + 1}$")
                axis.grid(True)
                axis.legend(loc=self.legend_loc)

            axes[-1, 0].set_xlabel("Time (s)")
            figure.suptitle(f"Experiment {experiment_index}: filtered outputs")
            figure.tight_layout()
            figure.savefig(
                plots_directory / f"experimental_filtered_data_{experiment_index}.svg",
                format="svg",
            )

    def compute(self):
        self.response = {}
        self.angular_frequencies = {}

        n_u = self.u[0].shape[0]
        n_y = self.y[0].shape[0]

        for i_u in range(n_u):
            u = self.u[i_u][i_u, :]

            t = self.t[i_u]
            dt = t[1] - t[0]
            f_s = 1 / dt
            N = len(u)

            for j_y in range(n_y):
                y = self.y[i_u][j_y, :]

                U = np.fft.fft(u)
                Y = np.fft.fft(y)
                S_uu = 1 / N * U * np.conj(U)
                S_yu = 1 / N * Y * np.conj(U)
                G = S_yu / S_uu
                frequency = np.fft.fftfreq(N, d=1 / f_s)
                self.response[f"u_{i_u + 1} -> y_{j_y + 1}"] = G
                self.angular_frequencies[f"u_{i_u + 1} -> y_{j_y + 1}"] = (
                    2 * np.pi * frequency
                )

        return self.angular_frequencies, self.response

    def plot_exp_data(self):
        plots_directory = Path(self.plots_dir)
        plots_directory.mkdir(parents=True, exist_ok=True)

        for experiment_index, (t, u, y) in enumerate(
            zip(self.t, self.u, self.y), start=1
        ):
            t_plot, u_plot = interpolate_signals(t, u)
            _, y_plot = interpolate_signals(t, y)
            figure, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True, dpi=200)

            for input_index, input_signal in enumerate(u_plot, start=1):
                axes[0].plot(t_plot, input_signal, label=rf"$u_{input_index}$")

            axes[0].set_title("Input signals")
            axes[0].set_ylabel("Amplitude")
            axes[0].grid(True)
            axes[0].legend(loc=self.legend_loc)

            for output_index, output_signal in enumerate(y_plot, start=1):
                axes[1].plot(t_plot, output_signal, label=rf"$y_{output_index}$")

            axes[1].set_title("Output signals")
            axes[1].set_xlabel("Time (s)")
            axes[1].set_ylabel("Amplitude")
            axes[1].grid(True)
            axes[1].legend(loc=self.legend_loc)

            figure.suptitle(f"Experiment {experiment_index}")
            figure.tight_layout()
            figure.savefig(
                plots_directory / f"experimental_data_{experiment_index}.svg",
                format="svg",
            )

    def plot_freq_resp(self, max_freq):
        if self.response is None or self.angular_frequencies is None:
            raise ValueError("compute must be called before plotting the response")
        if self.response.keys() != self.angular_frequencies.keys():
            raise ValueError("response and frequency relations must match")

        n_u = self.u[0].shape[0]
        n_y = self.y[0].shape[0]
        magnitude_figure, magnitude_axes = plt.subplots(
            n_u,
            n_y,
            figsize=(4 * n_y, 3 * n_u),
            sharex=True,
            squeeze=False,
            dpi=100,
        )
        phase_figure, phase_axes = plt.subplots(
            n_u,
            n_y,
            figsize=(4 * n_y, 3 * n_u),
            sharex=True,
            squeeze=False,
            dpi=100,
        )

        for relation, response in self.response.items():
            input_name, output_name = relation.split(" -> ")
            input_index = int(input_name.split("_")[1]) - 1
            output_index = int(output_name.split("_")[1]) - 1

            angular_frequency = self.angular_frequencies[relation]
            positive_frequency = (angular_frequency > 0) & (
                angular_frequency < max_freq
            )
            angular_frequency = angular_frequency[positive_frequency]
            response = response[positive_frequency]

            magnitude = 20 * np.log10(np.abs(response))
            phase = np.angle(response)
            plot_angular_frequency, magnitude_plot = interpolate_signals(
                angular_frequency, magnitude[np.newaxis, :]
            )
            _, phase_plot = interpolate_signals(angular_frequency, phase[np.newaxis, :])

            magnitude_axis = magnitude_axes[input_index, output_index]
            magnitude_axis.semilogx(plot_angular_frequency, magnitude_plot[0])
            magnitude_axis.set_title(relation)
            magnitude_axis.set_ylabel("Magnitude (dB)")
            magnitude_axis.grid(True, which="both")

            phase_axis = phase_axes[input_index, output_index]
            phase_axis.semilogx(plot_angular_frequency, phase_plot[0] * 180 / np.pi)
            phase_axis.set_title(relation)
            phase_axis.set_ylabel("Phase (degrees)")
            phase_axis.grid(True, which="both")

        for output_index in range(n_y):
            magnitude_axes[-1, output_index].set_xlabel("Angular frequency (rad/s)")
            phase_axes[-1, output_index].set_xlabel("Angular frequency (rad/s)")

        magnitude_figure.suptitle("Experimental frequency-response magnitude")
        magnitude_figure.tight_layout()
        phase_figure.suptitle("Experimental frequency-response phase")
        phase_figure.tight_layout()

        plots_directory = Path(self.plots_dir)
        plots_directory.mkdir(parents=True, exist_ok=True)
        magnitude_figure.savefig(
            plots_directory / "experimental_frequency_response_magnitude.svg",
            format="svg",
        )
        phase_figure.savefig(
            plots_directory / "experimental_frequency_response_phase.svg",
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
    transfer_function_matrix = control.tf(
        [[transfer_function, 0], [0, transfer_function]]
    )
    system = control.ss(transfer_function_matrix)

    initial_frequency = 0.5
    final_frequency = 300
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
    for input_index in range(2):
        u = np.zeros((2, n_samples))
        u[input_index] = chirp_signal
        y = np.asarray(control.forced_response(system, T=t, U=u).outputs)
        y += rng.normal(0.0, 0.005, size=y.shape)
        t_experiments.append(t)
        u_experiments.append(u)
        y_experiments.append(y)

    # Emprego da classe ExperimentalFrequencyResponse
    experimental_frequency_response = ExperimentalFrequencyResponse(
        t_experiments, u_experiments, y_experiments
    )
    experimental_frequency_response.plot_exp_data()
    experimental_frequency_response.compute()
    experimental_frequency_response.plot_freq_resp(max_freq=85)


if __name__ == "__main__":
    main()
    plt.show()
