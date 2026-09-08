from pathlib import Path
import control
import matplotlib.pyplot as plt
import numpy as np


class ExperimentalFrequencyResponse:
    def __init__(self, t_experiments, u_experiments, y_experiments):
        if not (len(t_experiments) == len(u_experiments) == len(y_experiments)):
            raise ValueError("t, u, and y must contain the same number of experiments")
        if not t_experiments:
            raise ValueError("at least one experiment is required")

        self.t = list(t_experiments)
        self.u = list(u_experiments)
        self.y = list(y_experiments)
        self.n_experiments = len(self.t)
        self.angular_frequencies = None
        self.response = None
        self.plots_dir = "../plots"

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
                G = S_uu / S_yu
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
            figure, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True, dpi=200)

            for input_index, input_signal in enumerate(u, start=1):
                axes[0].plot(t, input_signal, label=rf"$u_{input_index}$")
            axes[0].set_title("Input signals")
            axes[0].set_ylabel("Amplitude")
            axes[0].grid(True)
            axes[0].legend()

            for output_index, output_signal in enumerate(y, start=1):
                axes[1].plot(t, output_signal, label=rf"$y_{output_index}$")
            axes[1].set_title("Output signals")
            axes[1].set_xlabel("Time (s)")
            axes[1].set_ylabel("Amplitude")
            axes[1].grid(True)
            axes[1].legend()

            figure.suptitle(f"Experiment {experiment_index}")
            figure.tight_layout()
            figure.savefig(
                plots_directory / f"experimental_data_{experiment_index}.svg",
                format="svg",
            )

    def plot_freq_resp(self):
        pass


def main():
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

    experimental_frequency_response = ExperimentalFrequencyResponse(
        t_experiments, u_experiments, y_experiments
    )
    experimental_frequency_response.plot_exp_data()
    experimental_frequency_response.compute()


if __name__ == "__main__":
    main()
    plt.show()
