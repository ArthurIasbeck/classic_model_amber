import csv
from pathlib import Path
import numpy as np


def read_dspace_csv(filename):
    filename = Path(filename)

    with filename.open("r", encoding="cp1252", newline="") as file:
        reader = csv.reader(file)

        path_row = None
        trace_values_start = None

        for row in reader:

            if not row:
                continue

            if row[0] == "path":
                path_row = row

            elif row[0] == "trace_values":
                trace_values_start = row
                break

        if path_row is None:
            raise ValueError("Linha 'path' não encontrada no arquivo.")

        if trace_values_start is None:
            raise ValueError("Linha 'trace_values' não encontrada no arquivo.")

        paths = path_row[2:]
        paths = [path.strip() for path in paths]
        number_of_traces = len(paths)

        if number_of_traces == 0:
            raise ValueError("Nenhum trace foi encontrado.")

        data = [[float(value) for value in trace_values_start[1:]]]

        for row in reader:
            if not row:
                continue

            data.append([float(value) for value in row[1:]])

    data = np.asarray(data, dtype=np.float64)

    time = data[:, 0]
    values = data[:, 1:]

    if values.shape[1] != number_of_traces:
        raise ValueError(
            f"Número de traces inconsistente: "
            f"{number_of_traces} paths encontrados, "
            f"mas {values.shape[1]} colunas de dados."
        )

    trace_values = {path: values[:, i] for i, path in enumerate(paths)}
    return time, paths, trace_values


def main():
    filename = Path(__file__).resolve().parents[1] / "data" / "random_v13.csv"
    time, paths, trace_values = read_dspace_csv(filename)
    print(paths)


if __name__ == "__main__":
    main()
