from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import control


def chirp(t, initial_frequency, final_frequency, duration):
    """Generate a linear chirp with frequencies in hertz."""
    sweep_rate = (final_frequency - initial_frequency) / duration
    phase = 2 * np.pi * (initial_frequency * t + sweep_rate * t**2 / 2)
    return np.sin(phase)


def experimental_frequency_response(u, y, sample_period):
    """Compute G_hat from the DFT definitions in the problem statement."""
    n_samples = len(u)
    sample_indices = np.arange(n_samples)
    frequencies = np.fft.rfftfreq(n_samples, sample_period)
    angular_frequencies = 2 * np.pi * frequencies

    # The matrix multiplication is the direct sum over n in the equations.
    dft_kernel = np.exp(
        -1j * np.outer(angular_frequencies, sample_indices * sample_period)
    )
    u_spectrum = dft_kernel @ u
    y_spectrum = dft_kernel @ y
    suu = u_spectrum * np.conj(u_spectrum) / n_samples
    syu = y_spectrum * np.conj(u_spectrum) / n_samples

    response = np.full_like(syu, np.nan + 1j * np.nan)
    valid_bins = suu > np.finfo(float).eps
    response[valid_bins] = syu[valid_bins] / suu[valid_bins]
    return angular_frequencies, response


def main():
    sample_rate = 1000
    duration = 20
    n_samples = int(sample_rate * duration)
    time = np.arange(n_samples) / sample_rate

    natural_frequency = 2 * np.pi * 5.0
    damping_ratio = 0.08
    initial_frequency = 0.5
    final_frequency = 300

    system = control.TransferFunction(
        [natural_frequency**2],
        [1, 2 * damping_ratio * natural_frequency, natural_frequency**2],
    )

    rng = np.random.default_rng(7)
    input_signal = chirp(time, initial_frequency, final_frequency, duration)
    simulation = control.forced_response(system, T=time, U=input_signal)
    clean_output = np.asarray(simulation.outputs).squeeze()
    output_signal = clean_output + rng.normal(0.0, 0.0, size=n_samples)

    angular_frequencies, experimental_response = experimental_frequency_response(
        input_signal, output_signal, 1 / sample_rate
    )
    analytical_magnitude, analytical_phase, bode_frequencies = control.bode_plot(
        system, omega=angular_frequencies, plot=False
    )
    analytical_magnitude = np.asarray(analytical_magnitude).squeeze()
    analytical_phase = np.asarray(analytical_phase).squeeze()
    bode_frequencies = np.asarray(bode_frequencies).squeeze()

    positive_frequency = angular_frequencies > 0
    frequency_hz = bode_frequencies[positive_frequency] / (2 * np.pi)
    experimental_response = experimental_response[positive_frequency]
    analytical_magnitude = analytical_magnitude[positive_frequency]
    analytical_phase = analytical_phase[positive_frequency]

    figure, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True, dpi=200)
    axes[0].semilogx(
        frequency_hz,
        20 * np.log10(np.abs(experimental_response)),
        label="Experimental",
        alpha=0.8,
    )
    axes[0].semilogx(
        frequency_hz,
        20 * np.log10(analytical_magnitude),
        "--",
        label="Analytical",
    )
    axes[0].set_ylabel("Magnitude (dB)")
    axes[0].grid(True, which="both")
    axes[0].legend()

    axes[1].semilogx(
        frequency_hz,
        np.unwrap(np.angle(experimental_response)) * 180 / np.pi,
        label="Experimental",
        alpha=0.8,
    )
    axes[1].semilogx(
        frequency_hz,
        np.unwrap(analytical_phase) * 180 / np.pi,
        "--",
        label="Analytical",
    )
    axes[1].set_xlabel("Frequency (Hz)")
    axes[1].set_ylabel("Phase (degrees)")
    axes[1].grid(True, which="both")
    axes[1].legend()

    figure.suptitle("Experimental and analytical frequency responses")
    figure.tight_layout()

    plots_directory = Path(__file__).resolve().parent.parent / "plots"
    plots_directory.mkdir(exist_ok=True)
    figure.savefig(plots_directory / "experimental_frequency_response.svg")

    signal_figure, signal_axes = plt.subplots(
        2, 1, figsize=(8, 6), sharex=True, dpi=200
    )
    signal_axes[0].plot(time, input_signal, label="Input chirp")
    signal_axes[0].set_ylabel("Amplitude")
    signal_axes[0].set_title("Input signal")
    signal_axes[0].grid(True)
    signal_axes[0].legend()

    signal_axes[1].plot(time, output_signal, label="Noisy output")
    signal_axes[1].set_xlabel("Time (s)")
    signal_axes[1].set_ylabel("Amplitude")
    signal_axes[1].set_title("Output signal")
    signal_axes[1].grid(True)
    signal_axes[1].legend()

    signal_figure.suptitle("Signals used in the frequency-response calculation")
    signal_figure.tight_layout()
    signal_figure.savefig(plots_directory / "input_output_signals.svg")
    plt.show()


if __name__ == "__main__":
    main()
