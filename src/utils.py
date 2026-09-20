"""Funções auxiliares para interpolar sinais SISO e MIMO para plotagem."""

import numpy as np


def interpolate_signals_mimo(x, data):
    """Interpola cada sinal de uma matriz usando uma malha densa comum.

    Args:
        x: Eixo independente com uma amostra por coluna de ``data``.
        data: Matriz cujas linhas representam os sinais a interpolar.

    Returns:
        O eixo interpolado e a matriz de sinais interpolados.

    Raises:
        ValueError: Se ``data`` não for bidimensional ou tiver tamanho
            incompatível com ``x``.
    """
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
    """Interpola um sinal unidimensional usando uma malha densa.

    Args:
        x: Eixo independente do sinal.
        data: Sinal unidimensional a interpolar.

    Returns:
        O eixo interpolado e o sinal interpolado.

    Raises:
        ValueError: Se ``data`` não for unidimensional ou tiver tamanho
            incompatível com ``x``.
    """
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
    """Seleciona a interpolação SISO ou MIMO conforme a dimensão dos dados.

    Args:
        x: Eixo independente dos sinais.
        data: Sinal unidimensional ou matriz de sinais.

    Returns:
        O eixo interpolado e os dados interpolados.

    Raises:
        ValueError: Se ``data`` não for uni ou bidimensional.
    """
    data = np.asarray(data)
    if data.ndim == 1:
        return interpolate_signal_siso(x, data)
    if data.ndim == 2:
        return interpolate_signals_mimo(x, data)
    raise ValueError("data must be a one- or two-dimensional array")
