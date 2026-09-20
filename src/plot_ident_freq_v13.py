from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter


def _format_decimal_comma(value, _):
    return f"{value:g}".replace(".", ",")


def main():
    data_dir = Path(__file__).resolve().parent.parent / "data"
    plots_dir = data_dir.parent / "plots"
    plots_dir.mkdir(exist_ok=True)
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

    frequency_data = np.loadtxt(data_dir / "ident_freq_v13_freq.txt", delimiter=",")
    time_data = np.loadtxt(data_dir / "ident_freq_v13_time.txt", delimiter=",")

    w = frequency_data[:, 0]
    magnitude_data = frequency_data[:, 1]
    phase_data = frequency_data[:, 2]
    magnitude_sys = frequency_data[:, 3]
    phase_sys = frequency_data[:, 4]

    figure, axes = plt.subplots(2, 1, sharex=True, figsize=(8, 6))
    axes[0].semilogx(w, 20 * np.log10(magnitude_data), label="Experimento")
    axes[0].semilogx(w, 20 * np.log10(magnitude_sys), label="Modelo")
    axes[0].set_ylabel("Magnitude (dB)")
    axes[0].xaxis.set_major_formatter(FuncFormatter(_format_decimal_comma))
    axes[0].yaxis.set_major_formatter(FuncFormatter(_format_decimal_comma))
    axes[0].grid(True, which="both")
    axes[0].legend()

    axes[1].semilogx(w, phase_data, label="Experimento")
    axes[1].semilogx(w, phase_sys - 360, label="Modelo")
    axes[1].set_xlabel("Frequência (rad/s)")
    axes[1].set_ylabel("Fase (graus)")
    axes[1].xaxis.set_major_formatter(FuncFormatter(_format_decimal_comma))
    axes[1].yaxis.set_major_formatter(FuncFormatter(_format_decimal_comma))
    axes[1].grid(True, which="both")

    figure.suptitle("Diagrama de Bode")
    figure.tight_layout()
    figure.savefig(plots_dir / "ident_freq_v13_bode.svg", format="svg")

    t = time_data[:, 0]
    y = time_data[:, 1]
    y_real = time_data[:, 2]
    i = time_data[:, 3]
    i_real = time_data[:, 4]

    min_time = 200.2
    max_time = 201
    y = y[(t > min_time) & (t < max_time)]
    y_real = y_real[(t > min_time) & (t < max_time)]
    i = i[(t > min_time) & (t < max_time)]
    i_real = i_real[(t > min_time) & (t < max_time)]
    t = t[(t > min_time) & (t < max_time)]

    figure, axes = plt.subplots(2, 1, sharex=True, figsize=(8, 6))

    axes[0].plot(t, y, label="Modelo")
    axes[0].plot(t, y_real, label="Experimento")
    axes[0].set_ylabel(r"Deslocamento ($\mu$m)")
    axes[0].xaxis.set_major_formatter(FuncFormatter(_format_decimal_comma))
    axes[0].yaxis.set_major_formatter(FuncFormatter(_format_decimal_comma))
    axes[0].grid(True)
    axes[0].legend()

    axes[1].plot(t, i, label="Modelo")
    axes[1].plot(t, i_real, label="Experimento")
    axes[1].set_xlabel("Tempo (s)")
    axes[1].set_ylabel("Corrente (A)")
    axes[1].xaxis.set_major_formatter(FuncFormatter(_format_decimal_comma))
    axes[1].yaxis.set_major_formatter(FuncFormatter(_format_decimal_comma))
    axes[1].grid(True)
    axes[1].legend()

    figure.tight_layout()
    figure.savefig(plots_dir / "ident_freq_v13_time.svg", format="svg")

    plt.show()


if __name__ == "__main__":
    main()
