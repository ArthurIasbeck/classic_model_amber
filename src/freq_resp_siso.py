from dataset import Dataset
from experimental_frequency_response_siso import ExperimentalFrequencyResponseSiso
from loguru import logger
from matplotlib import pyplot as plt


def main():
    dataset_0 = Dataset(file_path="../data/chirp_v13_0.csv")
    dataset_1 = Dataset(file_path="../data/chirp_v13_1.txt")
    dataset_2 = Dataset(file_path="../data/chirp_v13_2.txt")
    dataset_3 = Dataset(file_path="../data/chirp_v13_3.txt")
    dataset_4 = Dataset(file_path="../data/chirp_v13_4.txt")

    logger.info("Iniciando carregamento dos dados...")
    t_0, i_0, d_0, y_0 = dataset_0.load()
    t_1, i_1, d_1, y_1 = dataset_1.load()
    t_2, i_2, d_2, y_2 = dataset_2.load()
    t_3, i_3, d_3, y_3 = dataset_3.load()
    t_4, i_4, d_4, y_4 = dataset_4.load()
    logger.info("Carregamento dos dados concluído.")

    t = [t_0, t_1, t_2, t_3, t_4]
    d = [d_0[0, :], d_1[0, :], d_2[0, :], d_3[0, :], d_4[0, :]]
    y = [y_0[0, :], y_1[0, :], y_2[0, :], y_3[0, :], y_4[0, :]]

    experimental_freq_resp = ExperimentalFrequencyResponseSiso("T_v13", t, d, y)

    logger.info("Iniciando filtragem dos dados...")
    experimental_freq_resp.filter_output(cutoff_frequency=600)
    experimental_freq_resp.plot_filter_output()
    logger.info("Filtragem dos dados concluída.")

    logger.info("Iniciando recorte dos dados...")
    experimental_freq_resp.clip_data(plot=False)
    experimental_freq_resp.plot_exp_data()
    logger.info("Recorte dos dados concluído.")

    logger.info("Iniciando computação da resposta em frequência...")
    experimental_freq_resp.compute()
    experimental_freq_resp.plot_freq_resp(min_freq=1, max_freq=600)
    logger.info("Computação da resposta em frequência concluída.")


if __name__ == "__main__":
    main()
    plt.show()
