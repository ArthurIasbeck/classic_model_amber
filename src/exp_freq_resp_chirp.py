import numpy as np
from matplotlib import pyplot as plt

from dataset import Dataset
from experimental_frequency_response import ExperimentalFrequencyResponse


def main():
    load_data_v13 = Dataset(file_path="../data/chirp_v13.txt")
    load_data_w13 = Dataset(file_path="../data/chirp_w13.txt")
    load_data_v24 = Dataset(file_path="../data/chirp_v24.txt")
    load_data_w24 = Dataset(file_path="../data/chirp_w24.txt")

    t_v13, i_v13, d_v13, y_v13 = load_data_v13.load()
    t_w13, i_w13, d_w13, y_w13 = load_data_w13.load()
    t_v24, i_v24, d_v24, y_v24 = load_data_v24.load()
    t_w24, i_w24, d_w24, y_w24 = load_data_w24.load()

    t = [
        np.transpose(t_v13),
        np.transpose(t_w13),
        np.transpose(t_v24),
        np.transpose(t_w24),
    ]
    i = [
        np.transpose(d_v13),
        np.transpose(d_w13),
        np.transpose(d_v24),
        np.transpose(d_w24),
    ]
    y = [
        np.transpose(y_v13),
        np.transpose(y_w13),
        np.transpose(y_v24),
        np.transpose(y_w24),
    ]

    experimental_frequency_response = ExperimentalFrequencyResponse(t, i, y)
    # experimental_frequency_response.plot_exp_data()

    experimental_frequency_response.filter_output()
    experimental_frequency_response.plot_filter_output()

    experimental_frequency_response.compute()
    experimental_frequency_response.plot_freq_resp(max_freq=1500)


if __name__ == "__main__":
    main()
    plt.show()
