"""Entrypoints para calcular respostas em frequência SISO experimentais."""

from dataset import Dataset
from experimental_frequency_response_siso import ExperimentalFrequencyResponseSiso
from loguru import logger
from matplotlib import pyplot as plt


def main_v13():
    """Processa os experimentos da direção V13 e calcula sua resposta SISO."""
    dataset_0 = Dataset(file_path="../data/chirp_v13_0.csv")
    dataset_1 = Dataset(file_path="../data/chirp_v13_1.txt")
    dataset_2 = Dataset(file_path="../data/chirp_v13_2.txt")
    dataset_3 = Dataset(file_path="../data/chirp_v13_3.txt")
    dataset_4 = Dataset(file_path="../data/chirp_v13_4.txt")

    logger.info("Iniciando carregamento dos dados...")
    t_0, _, d_0, y_0 = dataset_0.load()
    t_1, _, d_1, y_1 = dataset_1.load()
    t_2, _, d_2, y_2 = dataset_2.load()
    t_3, _, d_3, y_3 = dataset_3.load()
    t_4, _, d_4, y_4 = dataset_4.load()
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


def main_w13():
    """Processa os experimentos da direção W13 e calcula sua resposta SISO."""
    dataset_1 = Dataset(file_path="../data/chirp_w13_1.txt")
    dataset_2 = Dataset(file_path="../data/chirp_w13_2.txt")
    dataset_3 = Dataset(file_path="../data/chirp_w13_3.txt")
    dataset_4 = Dataset(file_path="../data/chirp_w13_4.txt")

    logger.info("Iniciando carregamento dos dados...")
    t_1, _, d_1, y_1 = dataset_1.load()
    t_2, _, d_2, y_2 = dataset_2.load()
    t_3, _, d_3, y_3 = dataset_3.load()
    t_4, _, d_4, y_4 = dataset_4.load()
    logger.info("Carregamento dos dados concluído.")

    t = [t_1, t_2, t_3, t_4]
    d = [d_1[1, :], d_2[1, :], d_3[1, :], d_4[1, :]]
    y = [y_1[1, :], y_2[1, :], y_3[1, :], y_4[1, :]]

    experimental_freq_resp = ExperimentalFrequencyResponseSiso("T_w13", t, d, y)

    logger.info("Iniciando filtragem dos dados...")
    experimental_freq_resp.filter_output(cutoff_frequency=600)
    experimental_freq_resp.plot_filter_output()
    logger.info("Filtragem dos dados concluída.")

    logger.info("Iniciando recorte dos dados...")
    experimental_freq_resp.clip_data(plot=True)
    experimental_freq_resp.plot_exp_data()
    logger.info("Recorte dos dados concluído.")

    logger.info("Iniciando computação da resposta em frequência...")
    experimental_freq_resp.compute()
    experimental_freq_resp.plot_freq_resp(min_freq=1, max_freq=600)
    logger.info("Computação da resposta em frequência concluída.")


if __name__ == "__main__":
    main_v13()
    main_w13()
    plt.show()
