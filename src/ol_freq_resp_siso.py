"""Entrypoints para calcular respostas experimentais em malha aberta."""

from matplotlib import pyplot as plt

from experimental_open_loop_frequency_response_siso import (
    ExperimentalOpenLoopFrequencyResponseSiso,
)


def main_v13_ol():
    """Calcula e analisa a resposta em malha aberta da direção V13."""
    open_loop_freq_response_siso = ExperimentalOpenLoopFrequencyResponseSiso("T_v13")
    open_loop_freq_response_siso.compute_open_loop_freq_resp(plot=True)
    open_loop_freq_response_siso.remove_delay(plot=True)


def main_w13_ol():
    """Calcula e analisa a resposta em malha aberta da direção W13."""
    open_loop_freq_response_siso = ExperimentalOpenLoopFrequencyResponseSiso("T_w13")
    open_loop_freq_response_siso.compute_open_loop_freq_resp(plot=True)
    open_loop_freq_response_siso.remove_delay(plot=True)


def main():
    """Executa as análises de malha aberta para V13 e W13."""
    main_v13_ol()
    main_w13_ol()


if __name__ == "__main__":
    main()
    plt.show()
