# Classic Model AMBER

Ferramentas de pesquisa para carregar, visualizar e analisar dados experimentais
de um sistema MIMO com quatro entradas, quatro perturbacoes e quatro saidas. O
repositorio tambem inclui uma estimativa de resposta em frequencia a partir dos
espectros de entrada e saida e exemplos para comparacao com modelos analiticos.

## Funcionalidades

- leitura de experimentos armazenados em arquivos de texto ou NumPy;
- separacao automatica dos sinais de entrada, perturbacao e saida;
- geracao de graficos dos sinais experimentais;
- estimativa da matriz de resposta em frequencia por FFT;
- graficos de magnitude e fase para cada relacao entrada-saida;
- analise de autocorrelacao e correlacao cruzada com decimacao configuravel;
- graficos das FFTs de perturbacoes e saidas;
- filtragem passa-baixas e analise de resposta em frequencia SISO;
- processamento separado dos ensaios chirp `V13` e `W13`;
- calculo da resposta de malha aberta SISO a partir da resposta de malha fechada;
- identificacao de respostas em frequencia e comparacao de modelos MIMO do rotor;
- construcao de modelos ROSS com diferentes malhas, suspensoes e configuracoes de FRF;
- scripts MATLAB para exportar os dados originais para o formato utilizado em
  Python;
- exemplo visual que compara a estimativa experimental com uma funcao de
  transferencia conhecida.

## Estrutura

```text
.
|-- plots/                              # Graficos gerados em SVG
|-- src/
|   |-- dataset.py                     # Carga e visualizacao dos dados
|   |-- experimental_frequency_response.py
|   |-- experimental_frequency_response_siso.py
|   |-- freq_resp.py                   # Analise MIMO com dados experimentais
|   |-- freq_resp_siso.py              # Analise SISO com dados experimentais
|   |-- experimental_open_loop_frequency_response_siso.py
|   |                                      # Classe para analise de malha aberta
|   |-- ol_freq_resp_siso.py             # Execucao das analises V13 e W13
|   |-- build_mimo_model_from_frequency_response.py
|   |                                      # Ajuste de modelo MIMO reduzido
|   |-- utils.py                       # Interpolacao para visualizacao
|   `-- utils/                         # Conversao, modelos e dados MATLAB
|-- tests/
|   |-- 00_compute_exp_freq_resp.py    # Comparacao visual da resposta
|   `-- 01_read_csv.py                 # Leitura de CSV do dSPACE
|-- requirements.txt
`-- README.md
```

O diretorio `data/` nao e versionado porque os arquivos experimentais podem
ocupar varios gigabytes.

## Instalacao

Crie um ambiente virtual e instale as dependencias fixadas no repositorio:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## Formato dos dados

Os arquivos de texto devem ser separados por virgula e conter 13 colunas nesta
ordem:

| Colunas | Conteudo |
| --- | --- |
| 1 | tempo |
| 2 a 5 | quatro entradas de corrente |
| 6 a 9 | quatro perturbacoes |
| 10 a 13 | quatro saidas |

Coloque os arquivos em `data/`. Ao receber, por exemplo,
`../data/chirp.txt`, `Dataset.load()` procura primeiro por
`../data/chirp.npz`. Se o cache nao existir, o arquivo de texto inteiro e
carregado e um arquivo NPZ com as chaves `t`, `i`, `d` e `y` e criado no mesmo
diretorio.

Os scripts em `src/utils/` convertem resultados MATLAB para essa organizacao de
colunas. Eles devem ser executados com o diretorio de trabalho ajustado para que
os caminhos relativos apontem para `data/`.

Os CSVs exportados pelo dSPACE podem conter cabecalhos e uma linha
`trace_values`; `Dataset.load()` localiza essa linha, reorganiza as colunas de
perturbacao e saida para a convencao do projeto e cria um cache `.npz` ao lado
do CSV. Os arquivos brutos e os caches permanecem fora do versionamento.

## Uso

### Visualizacao de um conjunto de dados

O exemplo em `dataset.py` usa `data/chirp.txt`. Execute-o a partir de `src/`,
pois o caminho do arquivo de entrada e relativo ao diretorio de trabalho:

```bash
cd src
python dataset.py
```

Para ambientes sem interface grafica:

```bash
cd src
MPLBACKEND=Agg python dataset.py
```

Os graficos de entrada, perturbacao e saida sao gravados em `plots/`.

`Dataset.plot_crosscorrelation()` e `Dataset.plot_autocorrelation()` geram
correlacoes com os sinais decimados. `Dataset.compute_fft()` gera as magnitudes
das FFTs de perturbacoes e saidas. Esses metodos exigem que o conjunto tenha
sido carregado e gravam os SVGs em `plots/`.

### Resposta em frequencia

`ExperimentalFrequencyResponse` recebe uma lista de experimentos por entrada
excitada. Cada matriz de entrada ou saida deve ter os sinais nas linhas e as
amostras nas colunas. O metodo `compute()` calcula, para cada par entrada-saida,

```text
S_uu = U * conj(U) / N
S_yu = Y * conj(U) / N
G = S_yu / S_uu
```

O `Dataset` segue o mesmo padrao: `i`, `d` e `y` sao armazenados como
`(n_sinais, n_amostras)`. Os arquivos CSV e TXT permanecem em formato tabular,
com uma amostra por linha; a transposicao ocorre durante o carregamento.

O exemplo principal cria um sistema MIMO sintetico de duas entradas e duas
saidas, excita uma entrada por experimento e gera os diagramas de magnitude e
fase:

```bash
cd src
python experimental_frequency_response.py
```

O argumento `max_freq` de `plot_freq_resp()` e os valores armazenados em
`angular_frequencies` usam radianos por segundo.

### Analise SISO e malha aberta

`freq_resp_siso.py` possui fluxos separados para os ensaios `V13` e `W13`. O
script `ol_freq_resp_siso.py` calcula as respostas de malha aberta dos dois
ensaios a partir das respostas de malha fechada salvas em `data/`:

```bash
cd src
MPLBACKEND=Agg python ol_freq_resp_siso.py
```

Os resultados sao salvos com o identificador do ensaio nos nomes dos arquivos
NPZ e dos graficos, evitando que os resultados de `V13` e `W13` sejam
sobrescritos.

Exemplos dos resultados versionados:

![Magnitude da resposta em frequencia](plots/experimental_frequency_response_magnitude.svg)

![Fase da resposta em frequencia](plots/experimental_frequency_response_phase.svg)

## Verificacao

Verifique a formatacao do codigo Python com:

```bash
python -m black --check src tests
```

O arquivo `tests/00_compute_exp_freq_resp.py` e um script de comparacao visual,
nao um teste automatizado com assercoes. Ele pode ser executado a partir da raiz
do repositorio:

```bash
MPLBACKEND=Agg python tests/00_compute_exp_freq_resp.py
```

Esse script constroi uma matriz DFT densa. Com os 20.000 pontos configurados no
exemplo, somente essa matriz pode consumir cerca de 3,2 GB de memoria. Evite
executa-lo em maquinas com memoria limitada.

Os scripts `src/freq_resp.py` e `src/freq_resp_siso.py` executam analises sobre
arquivos experimentais reais. O segundo aplica filtragem passa-baixas, recorta
os experimentos e calcula a resposta SISO media. Eles dependem dos dados locais
em `data/` e podem exigir varios gigabytes de memoria; nao sao testes
automatizados.

### Resposta de malha aberta

`experimental_open_loop_frequency_response_siso.py` implementa a classe que
carrega a resposta SISO de malha fechada salva em `data/`, calcula a resposta de
malha aberta usando o controlador definido no script e grava o resultado como
um novo arquivo NPZ. O script `ol_freq_resp_siso.py` e o ponto de entrada
executavel para os ensaios `V13` e `W13`; ele tambem gera graficos da resposta do
controlador, da resposta de malha aberta e da funcao objetivo usada para avaliar
a remocao de atraso.

O modulo `utils/amber_models.py` constroi modelos de rotor com a biblioteca
ROSS, permite comparar discretizacoes e calcula FRFs numericas para comparacao
com ensaios experimentais. O script `build_mimo_model_from_frequency_response.py`
tambem oferece um fluxo de ajuste de ganhos de um modelo MIMO reduzido a partir
das respostas de malha aberta salvas.
