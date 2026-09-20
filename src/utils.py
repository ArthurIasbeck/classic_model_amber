import numpy as np


def interpolate_signals_mimo(x, data):
    n_plot_points = 1_000_000
    x = np.asarray(x)
    data = np.asarray(data)
    if data.ndim != 2:
        raise ValueError("data must be a two-dimensional array")
    if data.shape[1] != x.size:
        raise ValueError("the independent axis and samples must have matching sizes")

    x_plot = np.linspace(x[0], x[-1], n_plot_points)
    interpolated_data = np.vstack([np.interp(x_plot, x, signal) for signal in data])

    return x_plot, interpolated_data


def interpolate_signal_siso(x, data):
    n_plot_points = 1_000_000
    x = np.asarray(x)
    data = np.asarray(data)
    if data.ndim != 1:
        raise ValueError("data must be a one-dimensional array")
    if data.size != x.size:
        raise ValueError("the independent axis and samples must have matching sizes")

    x_plot = np.linspace(x[0], x[-1], n_plot_points)
    interpolated_data = np.interp(x_plot, x, data)

    return x_plot, interpolated_data


def interpolate_signals(x, data):
    data = np.asarray(data)
    if data.ndim == 1:
        return interpolate_signal_siso(x, data)
    if data.ndim == 2:
        return interpolate_signals_mimo(x, data)
    raise ValueError("data must be a one- or two-dimensional array")
