from time import perf_counter

from matplotlib import pyplot as plt

from dataset import Dataset
from experimental_frequency_response import ExperimentalFrequencyResponse


def main():

    signal_name = "chirp"
    file_extension = ".txt"

    load_data_v13 = Dataset(file_path=f"../data/{signal_name}_v13{file_extension}")
    t_v13, i_v13, d_v13, y_v13 = load_data_v13.load()

    load_data_w13 = Dataset(file_path=f"../data/{signal_name}_w13{file_extension}")
    t_w13, i_w13, d_w13, y_w13 = load_data_w13.load()

    load_data_v24 = Dataset(file_path=f"../data/{signal_name}_v24{file_extension}")
    t_v24, i_v24, d_v24, y_v24 = load_data_v24.load()

    load_data_w24 = Dataset(file_path=f"../data/{signal_name}_w24{file_extension}")
    t_w24, i_w24, d_w24, y_w24 = load_data_w24.load()

    t = [
        t_v13,
        t_w13,
        t_v24,
        t_w24,
    ]
    d = [
        d_v13,
        d_w13,
        d_v24,
        d_w24,
    ]
    y = [
        y_v13,
        y_w13,
        y_v24,
        y_w24,
    ]

    experimental_frequency_response = ExperimentalFrequencyResponse(t, d, y)
    # experimental_frequency_response.plot_exp_data()
    # experimental_frequency_response.filter_output(cutoff_frequency=200)
    # experimental_frequency_response.plot_filter_output()
    experimental_frequency_response.compute()
    experimental_frequency_response.plot_freq_resp(max_freq=300)


if __name__ == "__main__":
    main()
    plt.show()
