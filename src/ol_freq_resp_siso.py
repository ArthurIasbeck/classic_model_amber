from matplotlib import pyplot as plt

from experimental_open_loop_frequency_response_siso import (
    ExperimentalOpenLoopFrequencyResponseSiso,
)


def main_v13():
    open_loop_freq_response_siso = ExperimentalOpenLoopFrequencyResponseSiso("T_v13")
    open_loop_freq_response_siso.compute_open_loop_freq_resp(plot=True)
    open_loop_freq_response_siso.remove_delay(plot=True)


def main_w13():
    open_loop_freq_response_siso = ExperimentalOpenLoopFrequencyResponseSiso("T_w13")
    open_loop_freq_response_siso.compute_open_loop_freq_resp(plot=True)
    open_loop_freq_response_siso.remove_delay(plot=True)


def main():
    main_v13()
    main_w13()


if __name__ == "__main__":
    main()
    plt.show()
