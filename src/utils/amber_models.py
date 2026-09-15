"""Modelo final do rotor e comparação com as FRFs experimentais.

As laminações dos dois rotores são representadas por elementos de eixo
anulares em paralelo com o eixo interno. Todas as opções de uso
estão reunidas na seção CONFIGURAÇÕES DO USUÁRIO.
"""

from dataclasses import dataclass
from pathlib import Path
import numpy as np
import plotly.colors as pcolors
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import ross as rs
from scipy.linalg import eigh

# =============================================================================
# CONFIGURAÇÕES DO USUÁRIO
# =============================================================================

# Construção e visualização do rotor
TIPO_MALHA = "refinada"  # mínima/intermediária/reduzida/refinada/completa
VALIDAR_REDUCAO_MALHA = False  # compara todas as malhas com a completa
INCLUIR_SUSPENSAO = False  # aproximação livre-livre adotada no modelo final
REPRESENTACAO_SUSPENSAO = "mola"  # "mola" ou "mancal_equivalente"
AMORTECIMENTO_SUSPENSAO = 0.0  # N.s/m por cordão
# Para o rotor sem cordões: "anterior", "compromisso" ou "f1_exato".
AJUSTE_RIGIDEZ_SEM_SUSPENSAO = "compromisso"
INCLUIR_AMORTECIMENTO = True  # amortecimento estrutural
# Instrumentação: os valores abaixo devem ser confirmados por certificado/pesagem.
INCLUIR_MASSAS_ACELEROMETROS = False
MASSA_ACELEROMETRO_DE = 0.004  # kg
MASSA_ACELEROMETRO_NDE = 0.004  # kg
PLOTAR_ROTOR = True
MOSTRAR_NOS_ROTOR = True
MALHAS_PARA_PLOTAR = None
# None plota apenas TIPO_MALHA. Para comparar várias geometrias, use por
# exemplo: ("intermediária", "reduzida", "refinada", "completa").

# Análise de FRF
CALCULAR_FRF_NUMERICA = False
COMPARAR_COM_EXPERIMENTAL = False
LADO_FRF = "NDE"  # opções: "DE" ou "NDE"
DIRECAO_FRF = "y"  # impacto e acelerômetros na parte superior
EXPERIMENTOS_FRF = None  # ("01", "02")  # use None para plotar todos os ensaios
FREQUENCIA_MAX_FRF = 1000.0  # Hz
NUM_PONTOS_FRF = 2561
# Razões obtidas por meia-potência nos cinco ensaios (valores adimensionais).
# Mantém o ajuste anterior como padrão; a opção modal fica para diagnóstico.
USAR_AMORTECIMENTO_MODAL_EXPERIMENTAL = False
RAZOES_AMORTECIMENTO_MODAL = np.array([0.005367, 0.001439, 0.000883])
# Corrige somente a convenção de sinal/fase; a amplitude não é alterada.
POLARIDADE_FRF_NUMERICA = -1.0
AVALIAR_ERRO_COMPLEXO_FRF = False
COERENCIA_MINIMA_ERRO_FRF = 0.80

# Saída dos gráficos
BIBLIOTECA_GRAFICA = "plotly"  # opções: "plotly" ou "matplotlib"
MOSTRAR_GRAFICOS = False
SALVAR_GRAFICOS = False
FORMATO_MATPLOTLIB = "pdf"  # opções usuais: "pdf", "png" ou "svg"
DPI_MATPLOTLIB = 300  # usado principalmente para arquivos PNG
PASTA_DADOS = Path(__file__).resolve().parent / "dados"
ARQUIVO_ROTOR = Path(__file__).resolve().parent / "rotor_my_ross.html"

# Nós de REFERÊNCIA DA MALHA COMPLETA. A função frf_numerica_modal converte
# automaticamente estes planos para a numeração local da malha selecionada.
NO_IMPACTO = 28  # disco central na malha completa
NOS_RESPOSTA_FRF = {"DE": 8, "NDE": 48}  # malha completa


# =============================================================================
# DADOS DO ROTOR E PARÂMETROS IDENTIFICADOS
# =============================================================================

MASSA_ALVO = 5.89  # kg
MASSA_SEM_ANEL = 4.39  # kg, desenho SKF
IP_ALVO = 0.004748  # kg.m²
IT_ALVO = 0.179010  # kg.m², diagnóstico; não imposto na identificação
LARGURA_ANEL = 0.025  # m
DIAMETRO_ANEL = 0.120  # m, usado na representação gráfica em Matplotlib
NO_DISCO_RESIDUAL = 28

FREQUENCIAS_ALVO = np.array([107.03125, 400.78125, 769.140625])  # Hz
E_EIXO_FRF = 190.46001209e9  # Pa
E_LAMINACAO_FRF = 200.0e9  # Pa
# Identificações específicas para a aproximação livre-livre sem cordões.
E_EIXO_LIVRE_COMPROMISSO = 192.32093468e9  # Pa, minimiza o maior erro dos 3 modos
E_EIXO_LIVRE_F1_EXATO = 194.37837736e9  # Pa, ajusta somente o modo de 107 Hz
KY_CORDAO_FRF = 60_265.12385  # N/m por apoio
KX_CORDAO_FRF = 0.0  # livre na direção horizontal
# Referências DE e NDE da malha completa. Também coincidem com os planos de
# atuação magnética preservados nas discretizações atuais.
NOS_CORDAO_FRF = (12, 44)  # posições efetivas identificadas
ALPHA_RAYLEIGH_FRF = 1.83039074  # 1/s
BETA_RAYLEIGH_FRF = 7.76478583e-8  # s

# Cores próximas, mas distintas, para ressaltar os dois materiais.
COR_EIXO = "#69757C"  # cinza aço
COR_LAMINACAO = "#db4b78"  # azul acinzentado
COR_DISCO = "#9c0c16"
COR_SUSPENSAO = "#B08A3E"


@dataclass(frozen=True)
class ElementoSuspensao:
    """Mola viscosa translacional entre um nó do rotor e o solo."""

    n: int
    kxx: float
    kyy: float
    cxx: float = 0.0
    cyy: float = 0.0
    tag: str = "Cordão equivalente"
    color: str = COR_SUSPENSAO


class RotorComSuspensao(rs.Rotor):
    """Rotor ROSS com molas ao solo sem usar ``BearingElement``."""

    def __init__(self, *args, elementos_suspensao=None, **kwargs):
        self.elementos_suspensao = list(elementos_suspensao or [])
        super().__init__(*args, **kwargs)

    def _adicionar_suspensao(self, matriz, atributo_x, atributo_y):
        matriz = matriz.copy()
        for elemento in self.elementos_suspensao:
            dof_x = elemento.n * self.number_dof
            dof_y = dof_x + 1
            matriz[dof_x, dof_x] += getattr(elemento, atributo_x)
            matriz[dof_y, dof_y] += getattr(elemento, atributo_y)
        return matriz

    def K(self, frequency):
        matriz = super().K(frequency)
        return self._adicionar_suspensao(matriz, "kxx", "kyy")

    def C(self, frequency):
        matriz = super().C(frequency)
        return self._adicionar_suspensao(matriz, "cxx", "cyy")


# Posições axiais e diâmetros do eixo principal.
POSICOES_NODAIS = np.array(
    [
        0.0,
        0.012,
        0.032,
        0.052,
        0.072,
        0.092,
        0.112,
        0.1208,
        0.1272,
        0.1348,
        0.1405,
        0.1469,
        0.1530,
        0.1592,
        0.1653,
        0.1804,
        0.1905,
        0.2063,
        0.2221,
        0.2379,
        0.2537,
        0.2695,
        0.2853,
        0.3011,
        0.3169,
        0.3243,
        0.3363,
        0.3580,
        0.3640,
        0.3705,
        0.3825,
        0.3986,
        0.4147,
        0.4308,
        0.4469,
        0.4630,
        0.4791,
        0.4952,
        0.5113,
        0.5274,
        0.5356,
        0.5457,
        0.5607,
        0.5669,
        0.5731,
        0.5792,
        0.5856,
        0.5913,
        0.5989,
        0.6053,
        0.6141,
        0.6341,
        0.6461,
    ]
)

# Cada valor corresponde ao elemento entre os nós i e i+1.
DIAMETROS_EIXO_MM = np.zeros(len(POSICOES_NODAIS) - 1)
DIAMETROS_EIXO_MM[0] = 6.35
DIAMETROS_EIXO_MM[1:5] = 32.0
DIAMETROS_EIXO_MM[5:14] = 34.8
DIAMETROS_EIXO_MM[14:16] = 49.9
DIAMETROS_EIXO_MM[16:27] = 19.05
DIAMETROS_EIXO_MM[27:29] = 54.0
DIAMETROS_EIXO_MM[29:40] = 19.05
DIAMETROS_EIXO_MM[40:42] = 49.9
DIAMETROS_EIXO_MM[42:51] = 34.8
DIAMETROS_EIXO_MM[51] = 6.35

# Elementos anulares que completam o diâmetro externo das laminações.
NOS_LAMINACAO = tuple(range(6, 14)) + tuple(range(42, 50))
DIAMETRO_INTERNO_LAMINACAO = 0.0348  # m
DIAMETRO_EXTERNO_LAMINACAO = 0.0498  # m

# A malha reduzida conserva transições geométricas e planos funcionais. Os
# números abaixo se referem à numeração da malha completa.
INDICES_MALHA_MINIMA = (
    0,
    1,
    5,
    6,
    8,
    12,
    14,
    16,
    25,
    27,
    28,
    29,
    40,
    42,
    44,
    48,
    50,
    51,
    52,
)
INDICES_MALHA_INTERMEDIARIA = (
    0,
    1,
    5,
    6,
    8,
    12,
    14,
    16,
    21,
    25,
    27,
    28,
    29,
    34,
    40,
    42,
    44,
    48,
    50,
    51,
    52,
)
INDICES_MALHA_REDUZIDA = (
    0,
    1,
    5,
    6,
    8,
    12,
    14,
    16,
    18,
    21,
    24,
    25,
    26,
    27,
    28,
    29,
    31,
    34,
    37,
    40,
    42,
    44,
    48,
    50,
    51,
    52,
)
INDICES_MALHA_REFINADA = (
    0,
    1,
    3,
    5,
    6,
    7,
    8,
    10,
    12,
    13,
    14,
    15,
    16,
    17,
    18,
    19,
    21,
    22,
    24,
    25,
    26,
    27,
    28,
    29,
    30,
    31,
    32,
    34,
    35,
    37,
    38,
    40,
    41,
    42,
    43,
    44,
    46,
    48,
    49,
    50,
    51,
    52,
)
MALHAS_CONVERGENCIA = {
    "mínima": INDICES_MALHA_MINIMA,
    "intermediária": INDICES_MALHA_INTERMEDIARIA,
    "reduzida": INDICES_MALHA_REDUZIDA,
    "refinada": INDICES_MALHA_REFINADA,
    "completa": tuple(range(len(POSICOES_NODAIS))),
}
PLANOS_FUNCIONAIS_COMPLETA = {
    "sensor_DE": NOS_RESPOSTA_FRF["DE"],
    "atuador_DE": NOS_CORDAO_FRF[0],
    "CG": 25,
    "disco": NO_IMPACTO,
    "atuador_NDE": NOS_CORDAO_FRF[1],
    "sensor_NDE": NOS_RESPOSTA_FRF["NDE"],
}


def _normalizar_tipo_malha(tipo_malha):
    tipo_malha = tipo_malha.lower().strip()
    aliases = {"minima": "mínima", "intermediaria": "intermediária"}
    tipo_malha = aliases.get(tipo_malha, tipo_malha)
    if tipo_malha not in MALHAS_CONVERGENCIA:
        opcoes = ", ".join(MALHAS_CONVERGENCIA)
        raise ValueError(f"tipo_malha deve ser uma destas opções: {opcoes}.")
    return tipo_malha


def _normalizar_representacao_suspensao(representacao):
    representacao = str(representacao).lower().strip()
    aliases = {
        "mola": "mola",
        "cordao": "mola",
        "cordão": "mola",
        "mancal": "mancal_equivalente",
        "mancal_equivalente": "mancal_equivalente",
    }
    try:
        return aliases[representacao]
    except KeyError as erro:
        raise ValueError(
            "representacao_suspensao deve ser 'mola' ou " "'mancal_equivalente'."
        ) from erro


def _modulo_eixo_identificado(incluir_suspensao, criterio=None):
    """Seleciona o módulo efetivo coerente com a condição de suspensão."""
    if incluir_suspensao:
        return E_EIXO_FRF
    criterio = AJUSTE_RIGIDEZ_SEM_SUSPENSAO if criterio is None else criterio
    criterio = str(criterio).lower().strip()
    valores = {
        "anterior": E_EIXO_FRF,
        "compromisso": E_EIXO_LIVRE_COMPROMISSO,
        "f1_exato": E_EIXO_LIVRE_F1_EXATO,
    }
    try:
        return valores[criterio]
    except KeyError as erro:
        raise ValueError(
            "AJUSTE_RIGIDEZ_SEM_SUSPENSAO deve ser 'anterior', "
            "'compromisso' ou 'f1_exato'."
        ) from erro


def _dados_malha(tipo_malha):
    """Retorna posições, diâmetros e intervalos laminados da malha escolhida."""
    tipo_malha = _normalizar_tipo_malha(tipo_malha)
    indices = np.asarray(MALHAS_CONVERGENCIA[tipo_malha])

    posicoes = POSICOES_NODAIS[indices]
    diametros = DIAMETROS_EIXO_MM[indices[:-1]] * 1e-3
    laminados = []
    for indice_i, indice_f in zip(indices[:-1], indices[1:]):
        laminados.append(
            (6 <= indice_i and indice_f <= 14) or (42 <= indice_i and indice_f <= 50)
        )

    nos_funcionais = {}
    for nome, no_completo in PLANOS_FUNCIONAIS_COMPLETA.items():
        encontrados = np.flatnonzero(np.isclose(posicoes, POSICOES_NODAIS[no_completo]))
        if len(encontrados) != 1:
            raise ValueError(f"O plano funcional {nome} não foi preservado.")
        nos_funcionais[nome] = int(encontrados[0])
    return posicoes, diametros, np.asarray(laminados), nos_funcionais


def _criar_elementos_eixo(
    E_eixo,
    E_laminacao,
    rho_eixo,
    rho_laminacao,
    alpha,
    beta,
    tipo_malha,
):
    """Cria o eixo principal e os anéis laminados sobrepostos."""
    material_eixo = rs.Material(
        name="eixo_aco",
        rho=rho_eixo,
        E=E_eixo,
        Poisson=0.30,
        color=COR_EIXO,
    )
    material_laminacao = rs.Material(
        name="laminacao_M20",
        rho=rho_laminacao,
        E=E_laminacao,
        Poisson=0.31,
        color=COR_LAMINACAO,
    )

    posicoes, diametros, intervalos_laminados, nos_funcionais = _dados_malha(tipo_malha)
    comprimentos = np.diff(posicoes)
    eixo = [
        rs.ShaftElement(
            n=no,
            L=comprimentos[no],
            idl=0.0,
            odl=diametros[no],
            material=material_eixo,
            shear_effects=True,
            rotary_inertia=True,
            gyroscopic=True,
            alpha=alpha,
            beta=beta,
            tag=f"Eixo {no}",
        )
        for no in range(len(comprimentos))
    ]
    laminacoes = [
        rs.ShaftElement(
            n=no,
            L=comprimentos[no],
            idl=DIAMETRO_INTERNO_LAMINACAO,
            odl=DIAMETRO_EXTERNO_LAMINACAO,
            material=material_laminacao,
            shear_effects=True,
            rotary_inertia=True,
            gyroscopic=True,
            alpha=alpha,
            beta=beta,
            tag=f"Laminação {no}",
        )
        for no in np.flatnonzero(intervalos_laminados)
    ]
    return eixo + laminacoes, posicoes, nos_funcionais


def montar_rotor(
    E_eixo=E_EIXO_FRF,
    E_laminacao=E_LAMINACAO_FRF,
    rho_eixo=7850.0,
    rho_laminacao=7850.0,
    incluir_suspensao=False,
    k_suspensao=(KX_CORDAO_FRF, KY_CORDAO_FRF),
    nos_suspensao=NOS_CORDAO_FRF,
    c_suspensao=AMORTECIMENTO_SUSPENSAO,
    alpha=0.0,
    beta=0.0,
    normalizar_massa_sem_anel=True,
    tipo_malha=TIPO_MALHA,
    massas_acelerometros=None,
    representacao_suspensao=REPRESENTACAO_SUSPENSAO,
):
    """Monta o rotor com ou sem a rigidez equivalente dos cordões.

    Quando ``normalizar_massa_sem_anel`` é verdadeiro, as densidades do eixo e
    das laminações recebem o mesmo fator para que sua massa conjunta seja
    4,39 kg. O disco central recebe a massa e o momento polar restantes. As
    massas opcionais dos acelerômetros são acrescentadas depois dessas
    restrições e, portanto, representam a instrumentação externa ao rotor de
    5,89 kg.
    """
    tipo_malha = _normalizar_tipo_malha(tipo_malha)
    representacao_suspensao = _normalizar_representacao_suspensao(
        representacao_suspensao
    )
    elementos, posicoes, nos_funcionais = _criar_elementos_eixo(
        E_eixo, E_laminacao, rho_eixo, rho_laminacao, alpha, beta, tipo_malha
    )

    fator_massa = 1.0
    if normalizar_massa_sem_anel:
        fator_massa = MASSA_SEM_ANEL / sum(elm.m for elm in elementos)
        elementos, posicoes, nos_funcionais = _criar_elementos_eixo(
            E_eixo,
            E_laminacao,
            rho_eixo * fator_massa,
            rho_laminacao * fator_massa,
            alpha,
            beta,
            tipo_malha,
        )

    massa_sem_anel = sum(elm.m for elm in elementos)
    ip_sem_anel = sum(elm.Im for elm in elementos)
    massa_anel = MASSA_ALVO - massa_sem_anel
    ip_anel = IP_ALVO - ip_sem_anel
    if massa_anel <= 0.0 or ip_anel <= 0.0:
        raise ValueError("Eixo e laminações excedem a massa ou o momento polar total.")

    id_anel = ip_anel / 2 + massa_anel * LARGURA_ANEL**2 / 12
    disco = rs.DiskElement(
        n=nos_funcionais["disco"],
        m=massa_anel,
        Ip=ip_anel,
        Id=id_anel,
        tag="Anel de massa removível",
        color=COR_DISCO,
    )

    apoios = []
    elementos_suspensao = []
    nos_suspensao_malha = []
    if incluir_suspensao:
        if np.isscalar(k_suspensao):
            kx = ky = float(k_suspensao)
        else:
            kx, ky = map(float, k_suspensao)
        for no_completo in nos_suspensao:
            encontrados = np.flatnonzero(
                np.isclose(posicoes, POSICOES_NODAIS[no_completo])
            )
            if len(encontrados) != 1:
                raise ValueError(
                    f"O nó de suspensão {no_completo} não existe na malha."
                )
            nos_suspensao_malha.append(int(encontrados[0]))
        elementos_suspensao = [
            ElementoSuspensao(
                n=no,
                kxx=kx,
                kyy=ky,
                cxx=c_suspensao,
                cyy=c_suspensao,
                tag=f"Cordão equivalente - nó {no}",
                color=COR_SUSPENSAO,
            )
            for no in nos_suspensao_malha
        ]
        if representacao_suspensao == "mancal_equivalente":
            apoios = [
                rs.BearingElement(
                    n=elemento.n,
                    kxx=elemento.kxx,
                    kyy=elemento.kyy,
                    cxx=elemento.cxx,
                    cyy=elemento.cyy,
                    tag=elemento.tag,
                    color=elemento.color,
                )
                for elemento in elementos_suspensao
            ]

    massas_acelerometros = massas_acelerometros or {}
    massas_sensores = {
        lado: float(massas_acelerometros.get(lado, 0.0)) for lado in ("DE", "NDE")
    }
    if any(massa < 0.0 for massa in massas_sensores.values()):
        raise ValueError("As massas dos acelerômetros não podem ser negativas.")
    discos_sensores = [
        rs.DiskElement(
            n=nos_funcionais[f"sensor_{lado}"],
            m=massa,
            Id=0.0,
            Ip=0.0,
            tag=f"Acelerômetro {lado}",
            color="#303030",
        )
        for lado, massa in massas_sensores.items()
        if massa > 0.0
    ]

    argumentos_rotor = {
        "shaft_elements": elementos,
        "disk_elements": [disco, *discos_sensores],
        "bearing_elements": apoios or None,
    }
    if incluir_suspensao and representacao_suspensao == "mola":
        rotor = RotorComSuspensao(
            **argumentos_rotor,
            elementos_suspensao=elementos_suspensao,
        )
    else:
        rotor = rs.Rotor(**argumentos_rotor)
    rotor.parametros_modelo = {
        "tipo_malha": tipo_malha,
        "numero_nos": len(posicoes),
        "nos_funcionais": nos_funcionais,
        "indices_referencia_completa": PLANOS_FUNCIONAIS_COMPLETA.copy(),
        "incluir_suspensao": incluir_suspensao,
        "representacao_suspensao": representacao_suspensao,
        "c_suspensao": float(c_suspensao),
        "nos_suspensao": tuple(nos_suspensao_malha),
        "nos_suspensao_referencia": tuple(nos_suspensao),
        "massas_acelerometros": massas_sensores,
        "fator_massa": fator_massa,
        "rho_eixo_efetivo": rho_eixo * fator_massa,
        "rho_laminacao_efetiva": rho_laminacao * fator_massa,
        "E_eixo": float(E_eixo),
        "E_laminacao": float(E_laminacao),
        "massa_sem_anel": massa_sem_anel,
        "massa_anel": massa_anel,
        "ip_anel": ip_anel,
    }
    return rotor


def rotor_identificado_frf(
    incluir_suspensao=True,
    incluir_amortecimento=True,
    tipo_malha=TIPO_MALHA,
    massas_acelerometros=None,
    representacao_suspensao=REPRESENTACAO_SUSPENSAO,
    c_suspensao=AMORTECIMENTO_SUSPENSAO,
    criterio_rigidez_sem_suspensao=None,
):
    """Retorna o modelo identificado para o ensaio de impacto vertical."""
    return montar_rotor(
        E_eixo=_modulo_eixo_identificado(
            incluir_suspensao, criterio_rigidez_sem_suspensao
        ),
        incluir_suspensao=incluir_suspensao,
        alpha=ALPHA_RAYLEIGH_FRF if incluir_amortecimento else 0.0,
        beta=BETA_RAYLEIGH_FRF if incluir_amortecimento else 0.0,
        tipo_malha=tipo_malha,
        massas_acelerometros=massas_acelerometros,
        representacao_suspensao=representacao_suspensao,
        c_suspensao=c_suspensao,
    )


def propriedades_de_massa(rotor):
    """Calcula massa, CG, Ip e momento transversal em torno do CG."""
    momento_estatico = 0.0
    it_origem = 0.0
    for elm in rotor.shaft_elements:
        x = elm.axial_cg_pos
        it_local = elm.Im / 2 + elm.m * elm.L**2 / 12
        momento_estatico += elm.m * x
        it_origem += it_local + elm.m * x**2
    for elm in rotor.disk_elements:
        x = rotor.nodes_pos[elm.n]
        momento_estatico += elm.m * x
        it_origem += elm.Id + elm.m * x**2
    for elm in rotor.point_mass_elements:
        x = rotor.nodes_pos[elm.n]
        momento_estatico += elm.m * x
        it_origem += elm.m * x**2

    cg = momento_estatico / rotor.m
    return {
        "massa": rotor.m,
        "CG": cg,
        "Ip": rotor.Ip,
        "It": it_origem - rotor.m * cg**2,
    }


def frequencias_flexiveis(rotor, direcao="y", alvos=FREQUENCIAS_ALVO):
    """Seleciona os modos flexíveis associados à direção solicitada."""
    direcao = direcao.lower()
    offset = {"x": 0, "y": 1}.get(direcao)
    if offset is None:
        raise ValueError("direcao deve ser 'x' ou 'y'.")

    autovalores, modos = eigh(rotor.K(0), rotor.M(0), check_finite=False, driver="gvd")
    positivos = autovalores > (2 * np.pi * 1.0) ** 2
    frequencias = np.sqrt(autovalores[positivos]) / (2 * np.pi)
    modos = modos[:, positivos]
    dofs = np.arange(offset, rotor.ndof, rotor.number_dof)
    participacao = np.sum(np.abs(modos[dofs, :]) ** 2, axis=0)

    selecionadas = []
    for alvo in np.asarray(alvos):
        candidatos = np.argsort(np.abs(frequencias - alvo))[:2]
        indice = candidatos[np.argmax(participacao[candidatos])]
        selecionadas.append(frequencias[indice])
    return np.asarray(selecionadas)


def frequencias_flexiveis_y(rotor, alvos=FREQUENCIAS_ALVO):
    """Atalho mantido para compatibilidade com as versões anteriores."""
    return frequencias_flexiveis(rotor, direcao="y", alvos=alvos)


def _normalizar_biblioteca(biblioteca=None):
    biblioteca = BIBLIOTECA_GRAFICA if biblioteca is None else biblioteca
    biblioteca = biblioteca.lower().strip()
    if biblioteca not in {"plotly", "matplotlib"}:
        raise ValueError("biblioteca deve ser 'plotly' ou 'matplotlib'.")
    return biblioteca


def _caminho_grafico(caminho_base, biblioteca=None):
    """Aplica a extensão adequada ao backend escolhido."""
    biblioteca = _normalizar_biblioteca(biblioteca)
    caminho_base = Path(caminho_base)
    extensao = ".html" if biblioteca == "plotly" else f".{FORMATO_MATPLOTLIB}"
    return caminho_base.with_suffix(extensao)


def _plotar_rotor_plotly(rotor, mostrar, mostrar_nos, arquivo_saida):
    fig = rotor.plot_rotor(nodes=1)
    if not mostrar_nos:
        fig.data[0].visible = False
    for no_apoio in rotor.parametros_modelo["nos_suspensao"]:
        fig.add_vline(
            x=float(rotor.nodes_pos[no_apoio]),
            line_color=COR_SUSPENSAO,
            line_width=1.2,
            line_dash="dash",
        )
    fig.update_xaxes(title_text="Posição axial [m]")
    fig.update_yaxes(title_text="Raio [m]")
    tipo_malha = rotor.parametros_modelo["tipo_malha"]
    fig.update_layout(
        title=f"Rotor - malha {tipo_malha}",
        height=620,
        showlegend=False,
    )
    if arquivo_saida is not None:
        fig.write_html(str(arquivo_saida), include_plotlyjs=True)
        print(f"Gráfico do rotor salvo em: {arquivo_saida}")
    if mostrar:
        fig.show()
    return fig


def _plotar_rotor_matplotlib(rotor, mostrar, mostrar_nos, arquivo_saida):
    import matplotlib

    if not mostrar:
        matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    fig, ax = plt.subplots(figsize=(12, 4.2), constrained_layout=True)
    ax.axhline(0.0, color="0.25", linewidth=0.8, linestyle="-.", zorder=0)

    # Elementos sólidos e anulares são desenhados separadamente para que o
    # eixo interno permaneça visível dentro dos pacotes laminados.
    for elm in rotor.shaft_elements:
        x0 = float(rotor.nodes_pos[elm.n])
        comprimento = float(elm.L)
        raio_externo = float(elm.odl) / 2
        raio_interno = float(elm.idl) / 2
        cor = elm.material.color
        zorder = 3 if raio_interno > 0 else 2
        if raio_interno == 0:
            ax.add_patch(
                Rectangle(
                    (x0, -raio_externo),
                    comprimento,
                    2 * raio_externo,
                    facecolor=cor,
                    edgecolor="black",
                    linewidth=0.55,
                    alpha=0.9,
                    zorder=zorder,
                )
            )
        else:
            espessura = raio_externo - raio_interno
            for y0 in (raio_interno, -raio_externo):
                ax.add_patch(
                    Rectangle(
                        (x0, y0),
                        comprimento,
                        espessura,
                        facecolor=cor,
                        edgecolor="black",
                        linewidth=0.55,
                        alpha=0.95,
                        zorder=zorder,
                    )
                )

    no_disco = rotor.parametros_modelo["nos_funcionais"]["disco"]
    x_disco = float(rotor.nodes_pos[no_disco])
    ax.add_patch(
        Rectangle(
            (x_disco - LARGURA_ANEL / 2, -DIAMETRO_ANEL / 2),
            LARGURA_ANEL,
            DIAMETRO_ANEL,
            facecolor=COR_DISCO,
            edgecolor="black",
            linewidth=0.8,
            alpha=0.9,
            zorder=4,
        )
    )

    for no_apoio in rotor.parametros_modelo["nos_suspensao"]:
        x_apoio = float(rotor.nodes_pos[no_apoio])
        ax.axvline(
            x_apoio,
            color=COR_SUSPENSAO,
            linewidth=1.2,
            linestyle="--",
            alpha=0.9,
            zorder=1,
        )

    if mostrar_nos:
        y_nos = -DIAMETRO_ANEL / 2 - 0.008
        ax.scatter(
            rotor.nodes_pos,
            np.full(len(rotor.nodes_pos), y_nos),
            s=7,
            color="black",
            zorder=5,
        )
        for no, x_no in enumerate(rotor.nodes_pos):
            ax.text(
                x_no,
                y_nos - 0.003,
                str(no),
                rotation=90,
                ha="center",
                va="top",
                fontsize=5.5,
            )

    tipo_malha = rotor.parametros_modelo["tipo_malha"]
    ax.set_title(f"Rotor - malha {tipo_malha}")
    ax.set_xlabel("Posição axial [m]")
    ax.set_ylabel("Raio [m]")
    ax.set_xlim(POSICOES_NODAIS[0] - 0.015, POSICOES_NODAIS[-1] + 0.015)
    limite_y = DIAMETRO_ANEL / 2 + (0.025 if mostrar_nos else 0.012)
    ax.set_ylim(-limite_y, limite_y)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(axis="x", color="0.88", linewidth=0.6)
    ax.spines[["top", "right"]].set_visible(False)

    if arquivo_saida is not None:
        fig.savefig(
            arquivo_saida,
            dpi=DPI_MATPLOTLIB,
            bbox_inches="tight",
            facecolor="white",
        )
        print(f"Gráfico do rotor salvo em: {arquivo_saida}")
    if mostrar:
        plt.show()
    return fig


def plotar_rotor(
    rotor,
    mostrar=True,
    mostrar_nos=True,
    arquivo_saida=None,
    biblioteca=None,
):
    """Plota o rotor completo com Plotly ou Matplotlib."""
    biblioteca = _normalizar_biblioteca(biblioteca)
    if biblioteca == "plotly":
        return _plotar_rotor_plotly(rotor, mostrar, mostrar_nos, arquivo_saida)
    return _plotar_rotor_matplotlib(rotor, mostrar, mostrar_nos, arquivo_saida)


def plotar_rotores_malhas(
    tipos_malha=None,
    incluir_suspensao=INCLUIR_SUSPENSAO,
    representacao_suspensao=REPRESENTACAO_SUSPENSAO,
    incluir_amortecimento=INCLUIR_AMORTECIMENTO,
    mostrar=True,
    mostrar_nos=True,
    salvar=False,
    biblioteca=None,
    caminho_base=ARQUIVO_ROTOR,
    rotor_atual=None,
    massas_acelerometros=None,
    c_suspensao=AMORTECIMENTO_SUSPENSAO,
):
    """Plota uma ou várias discretizações mantendo as mesmas condições físicas.

    ``tipos_malha=None`` plota somente ``TIPO_MALHA``. Uma string seleciona uma
    única malha e uma sequência permite gerar várias figuras na mesma execução.
    """
    if tipos_malha is None:
        tipos_malha = (TIPO_MALHA,)
    elif isinstance(tipos_malha, str):
        tipos_malha = (tipos_malha,)

    # Normaliza, valida e remove repetições sem alterar a ordem solicitada.
    tipos_normalizados = tuple(
        dict.fromkeys(_normalizar_tipo_malha(tipo) for tipo in tipos_malha)
    )
    if not tipos_normalizados:
        raise ValueError("tipos_malha deve conter pelo menos uma malha.")

    biblioteca = _normalizar_biblioteca(biblioteca)
    representacao_suspensao = _normalizar_representacao_suspensao(
        representacao_suspensao
    )
    caminho_base = Path(caminho_base)
    massas_sensores = {
        lado: float((massas_acelerometros or {}).get(lado, 0.0))
        for lado in ("DE", "NDE")
    }
    figuras = {}
    for tipo in tipos_normalizados:
        reutilizar_atual = (
            rotor_atual is not None
            and rotor_atual.parametros_modelo["tipo_malha"] == tipo
            and rotor_atual.parametros_modelo["incluir_suspensao"] == incluir_suspensao
            and rotor_atual.parametros_modelo["representacao_suspensao"]
            == representacao_suspensao
            and rotor_atual.parametros_modelo["massas_acelerometros"] == massas_sensores
            and np.isclose(rotor_atual.parametros_modelo["c_suspensao"], c_suspensao)
        )
        rotor_malha = (
            rotor_atual
            if reutilizar_atual
            else rotor_identificado_frf(
                incluir_suspensao=incluir_suspensao,
                incluir_amortecimento=incluir_amortecimento,
                tipo_malha=tipo,
                massas_acelerometros=massas_sensores,
                representacao_suspensao=representacao_suspensao,
                c_suspensao=c_suspensao,
            )
        )

        arquivo_saida = None
        if salvar:
            base_malha = caminho_base
            if len(tipos_normalizados) > 1:
                nome_malha = tipo.replace("í", "i").replace("á", "a")
                base_malha = caminho_base.with_name(f"{caminho_base.stem}_{nome_malha}")
            arquivo_saida = _caminho_grafico(base_malha, biblioteca)

        figuras[tipo] = plotar_rotor(
            rotor_malha,
            mostrar=False,
            mostrar_nos=mostrar_nos,
            arquivo_saida=arquivo_saida,
            biblioteca=biblioteca,
        )

    if mostrar:
        if biblioteca == "plotly":
            for figura in figuras.values():
                figura.show()
        else:
            import matplotlib.pyplot as plt

            plt.show()
    return figuras


def _normalizar_lado(lado):
    lado = lado.upper().strip()
    if lado not in NOS_RESPOSTA_FRF:
        raise ValueError("lado deve ser 'DE' ou 'NDE'.")
    return lado


def frf_numerica_modal(
    rotor,
    lado="DE",
    direcao="y",
    no_impacto=None,
    no_resposta=None,
    freq_max_hz=FREQUENCIA_MAX_FRF,
    num_pontos=NUM_PONTOS_FRF,
    razoes_amortecimento=None,
    polaridade=1.0,
):
    """Calcula a acelerância entre o disco e o lado DE ou NDE."""
    lado = _normalizar_lado(lado)
    direcao = direcao.lower()
    offset = {"x": 0, "y": 1}.get(direcao)
    if offset is None:
        raise ValueError("direcao deve ser 'x' ou 'y'.")
    if no_impacto is None:
        no_impacto = rotor.parametros_modelo["nos_funcionais"]["disco"]
    if no_resposta is None:
        no_resposta = rotor.parametros_modelo["nos_funcionais"][f"sensor_{lado}"]

    frequencias = np.linspace(0.0, freq_max_hz, num_pontos)
    omega = 2 * np.pi * frequencias
    M = rotor.M(0)
    K = rotor.K(0)
    C = rotor.C(0)
    autovalores, modos = eigh(K, M, check_finite=False)
    limite_rigido = (2 * np.pi * 0.1) ** 2
    rigidos = np.abs(autovalores) <= limite_rigido
    if np.any(autovalores < -limite_rigido):
        raise ValueError("O modelo possui autovalores negativos não rígidos.")
    # O resíduo dos modos rígidos deve permanecer na acelerância livre-livre.
    # Apenas seus pequenos autovalores numéricos são zerados.
    autovalores = autovalores.copy()
    autovalores[rigidos] = 0.0
    amortecimento_modal = np.real(np.diag(modos.T @ C @ modos)).copy()

    if razoes_amortecimento is not None:
        razoes_amortecimento = np.asarray(razoes_amortecimento, dtype=float)
        if len(razoes_amortecimento) != len(FREQUENCIAS_ALVO):
            raise ValueError(
                "razoes_amortecimento deve ter um valor para cada frequência alvo."
            )
        if np.any(razoes_amortecimento < 0.0):
            raise ValueError("As razões de amortecimento não podem ser negativas.")
        positivos = autovalores > limite_rigido
        frequencias_naturais = np.sqrt(autovalores[positivos]) / (2 * np.pi)
        indices_positivos = np.flatnonzero(positivos)
        dofs_direcao = np.arange(offset, rotor.ndof, rotor.number_dof)
        participacao_direcional = np.sum(
            np.abs(modos[dofs_direcao, :][:, positivos]) ** 2, axis=0
        )
        for alvo, zeta in zip(FREQUENCIAS_ALVO, razoes_amortecimento):
            candidatos = np.argsort(np.abs(frequencias_naturais - alvo))[:2]
            indice_local = candidatos[np.argmax(participacao_direcional[candidatos])]
            indice = indices_positivos[indice_local]
            amortecimento_modal[indice] = 2 * zeta * np.sqrt(autovalores[indice])

    dof_impacto = no_impacto * rotor.number_dof + offset
    dof_resposta = no_resposta * rotor.number_dof + offset
    participacao = modos[dof_resposta, :] * modos[dof_impacto, :]
    acelerancia = np.zeros_like(omega, dtype=complex)
    frequencias_nao_nulas = omega > 0.0
    omega_ativos = omega[frequencias_nao_nulas]
    denominador = (
        autovalores[:, None]
        - omega_ativos[None, :] ** 2
        + 1j * amortecimento_modal[:, None] * omega_ativos[None, :]
    )
    receptancia = np.sum(participacao[:, None] / denominador, axis=0)
    acelerancia[frequencias_nao_nulas] = float(polaridade) * (
        -(omega_ativos**2) * receptancia
    )

    return {
        "frequencia": frequencias,
        "amplitude": np.abs(acelerancia) / 9.80665,
        "fase": np.angle(acelerancia, deg=True),
        "lado": lado,
        "direcao": direcao,
        "no_impacto": no_impacto,
        "no_resposta": no_resposta,
        "no_impacto_referencia": NO_IMPACTO,
        "no_resposta_referencia": NOS_RESPOSTA_FRF[lado],
        "tipo_malha": rotor.parametros_modelo["tipo_malha"],
        "razoes_amortecimento": (
            None if razoes_amortecimento is None else razoes_amortecimento.copy()
        ),
        "polaridade": float(polaridade),
        "numero_modos_rigidos": int(np.count_nonzero(rigidos)),
    }


def _ler_arquivo_colunas(caminho, numero_colunas):
    linhas = []
    with open(caminho, "r", encoding="utf-8") as arquivo:
        for linha in arquivo:
            partes = linha.strip().split()
            if len(partes) != numero_colunas:
                continue
            try:
                linhas.append([float(v.replace(",", ".")) for v in partes])
            except ValueError:
                continue
    return np.asarray(linhas, dtype=float)


def carregar_frf_experimental(
    pasta=PASTA_DADOS,
    lado="DE",
    experimentos=None,
):
    """Carrega amplitude, fase e coerência somente do lado escolhido."""
    pasta = Path(pasta)
    lado = _normalizar_lado(lado)
    filtro = None if experimentos is None else {str(i) for i in experimentos}
    experimentos = {}
    for caminho_amp in sorted(pasta.glob(f"*_amp_phase_{lado.lower()}.txt")):
        identificador = caminho_amp.stem.split("_")[0]
        if filtro is not None and identificador not in filtro:
            continue
        caminho_coe = pasta / f"{identificador}_coe_{lado.lower()}.txt"
        amp_fase = _ler_arquivo_colunas(caminho_amp, 3)
        coerencia = (
            _ler_arquivo_colunas(caminho_coe, 2)
            if caminho_coe.exists()
            else np.empty((0, 2))
        )
        if amp_fase.size:
            experimentos[identificador] = {
                "frequencia": amp_fase[:, 0],
                "amplitude": amp_fase[:, 1],
                "fase": amp_fase[:, 2],
                "frequencia_coerencia": (
                    coerencia[:, 0] if coerencia.size else np.array([])
                ),
                "coerencia": coerencia[:, 1] if coerencia.size else np.array([]),
            }
    if not experimentos:
        raise FileNotFoundError(
            f"Nenhuma FRF experimental do lado {lado} encontrada em {pasta}."
        )
    return experimentos


def avaliar_erro_complexo_frf(
    dados_numericos,
    dados_experimentais,
    coerencia_minima=0.80,
    frequencia_minima=1.0,
    frequencia_maxima=None,
):
    """Calcula erro complexo normalizado usando a coerência como peso."""
    if not 0.0 <= coerencia_minima <= 1.0:
        raise ValueError("coerencia_minima deve estar entre zero e um.")

    freq_num = dados_numericos["frequencia"]
    h_num = dados_numericos["amplitude"] * np.exp(
        1j * np.deg2rad(dados_numericos["fase"])
    )
    limite_superior = (
        float(freq_num[-1]) if frequencia_maxima is None else float(frequencia_maxima)
    )

    resultados = {}
    for identificador, dados in dados_experimentais.items():
        freq_exp = dados["frequencia"]
        h_exp = dados["amplitude"] * np.exp(1j * np.deg2rad(dados["fase"]))
        h_num_exp = np.interp(freq_exp, freq_num, h_num.real) + 1j * np.interp(
            freq_exp, freq_num, h_num.imag
        )
        if dados["coerencia"].size:
            coerencia = np.interp(
                freq_exp,
                dados["frequencia_coerencia"],
                dados["coerencia"],
            )
        else:
            coerencia = np.ones_like(freq_exp)

        mascara = (
            (freq_exp >= frequencia_minima)
            & (freq_exp <= limite_superior)
            & (coerencia >= coerencia_minima)
            & np.isfinite(h_exp)
            & np.isfinite(h_num_exp)
        )
        if not np.any(mascara):
            erro = np.nan
        else:
            numerador = np.sum(
                coerencia[mascara] * np.abs(h_num_exp[mascara] - h_exp[mascara]) ** 2
            )
            denominador = np.sum(coerencia[mascara] * np.abs(h_exp[mascara]) ** 2)
            erro = 100 * np.sqrt(numerador / denominador)
        resultados[identificador] = {
            "erro_complexo_percentual": float(erro),
            "numero_pontos": int(np.count_nonzero(mascara)),
        }

    erros_validos = [
        dados["erro_complexo_percentual"]
        for dados in resultados.values()
        if np.isfinite(dados["erro_complexo_percentual"])
    ]
    resultados["mediana_percentual"] = (
        float(np.median(erros_validos)) if erros_validos else np.nan
    )
    return resultados


def _plotar_frf_plotly(
    dados_numericos,
    dados_experimentais=None,
    mostrar=True,
    arquivo_saida=None,
):
    """Implementação interativa da FRF em Plotly."""
    lado = dados_numericos["lado"]
    fig = make_subplots(rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.07)

    if dados_experimentais:
        cores = pcolors.qualitative.Plotly
        for indice, identificador in enumerate(sorted(dados_experimentais)):
            dados = dados_experimentais[identificador]
            cor = cores[indice % len(cores)]
            grupo = f"experimental-{identificador}-{lado}"
            fig.add_trace(
                go.Scatter(
                    x=dados["frequencia"],
                    y=dados["amplitude"],
                    name=f"Exp. {identificador} {lado}",
                    legendgroup=grupo,
                    line=dict(color=cor),
                ),
                row=1,
                col=1,
            )
            fig.add_trace(
                go.Scatter(
                    x=dados["frequencia"],
                    y=dados["fase"],
                    showlegend=False,
                    legendgroup=grupo,
                    line=dict(color=cor),
                ),
                row=2,
                col=1,
            )
            if dados["coerencia"].size:
                fig.add_trace(
                    go.Scatter(
                        x=dados["frequencia_coerencia"],
                        y=dados["coerencia"],
                        showlegend=False,
                        legendgroup=grupo,
                        line=dict(color=cor),
                    ),
                    row=3,
                    col=1,
                )

    grupo_numerico = f"numerico-{lado}"
    nome_numerico = (
        f"Numérico {lado}: disco → sensor {lado} " f"({dados_numericos['tipo_malha']})"
    )
    fig.add_trace(
        go.Scatter(
            x=dados_numericos["frequencia"],
            y=dados_numericos["amplitude"],
            name=nome_numerico,
            legendgroup=grupo_numerico,
            line=dict(color="black", dash="dot", width=3),
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=dados_numericos["frequencia"],
            y=dados_numericos["fase"],
            showlegend=False,
            legendgroup=grupo_numerico,
            line=dict(color="black", dash="dot", width=3),
        ),
        row=2,
        col=1,
    )

    titulo = f"FRF numérica - lado {lado}"
    if dados_experimentais:
        titulo = f"FRF numérica e experimental - lado {lado}"
    fig.update_layout(title=titulo, height=850, legend_title="Curvas")
    fig.update_yaxes(title_text="Amplitude [g/N]", type="log", row=1, col=1)
    fig.update_yaxes(title_text="Fase [graus]", row=2, col=1)
    fig.update_yaxes(title_text="Coerência", range=[0, 1.05], row=3, col=1)
    fig.update_xaxes(
        title_text="Frequência [Hz]",
        range=[0, float(dados_numericos["frequencia"][-1])],
        row=3,
        col=1,
    )

    if arquivo_saida is not None:
        fig.write_html(str(arquivo_saida), include_plotlyjs=True)
        print(f"Gráfico de FRF salvo em: {arquivo_saida}")
    if mostrar:
        fig.show()
    return fig


def _plotar_frf_matplotlib(
    dados_numericos,
    dados_experimentais=None,
    mostrar=True,
    arquivo_saida=None,
):
    """Implementação da FRF em Matplotlib para publicação."""
    import matplotlib

    if not mostrar:
        matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    lado = dados_numericos["lado"]
    fig, eixos = plt.subplots(
        3,
        1,
        figsize=(8.0, 8.5),
        sharex=True,
        gridspec_kw={"height_ratios": [1.4, 1.0, 0.8]},
        constrained_layout=True,
    )
    ax_amp, ax_fase, ax_coe = eixos

    if dados_experimentais:
        mapa_cores = plt.get_cmap("tab10")
        for indice, identificador in enumerate(sorted(dados_experimentais)):
            dados = dados_experimentais[identificador]
            cor = mapa_cores(indice % 10)
            ax_amp.semilogy(
                dados["frequencia"],
                dados["amplitude"],
                color=cor,
                linewidth=1.0,
                label=f"Exp. {identificador} {lado}",
            )
            ax_fase.plot(
                dados["frequencia"],
                dados["fase"],
                color=cor,
                linewidth=0.9,
            )
            if dados["coerencia"].size:
                ax_coe.plot(
                    dados["frequencia_coerencia"],
                    dados["coerencia"],
                    color=cor,
                    linewidth=0.9,
                )

    rotulo_numerico = (
        f"Numérico {lado}: disco → sensor {lado} " f"({dados_numericos['tipo_malha']})"
    )
    ax_amp.semilogy(
        dados_numericos["frequencia"],
        dados_numericos["amplitude"],
        color="black",
        linestyle="--",
        linewidth=1.6,
        label=rotulo_numerico,
    )
    ax_fase.plot(
        dados_numericos["frequencia"],
        dados_numericos["fase"],
        color="black",
        linestyle="--",
        linewidth=1.4,
    )

    titulo = f"FRF numérica - lado {lado}"
    if dados_experimentais:
        titulo = f"FRF numérica e experimental - lado {lado}"
    fig.suptitle(titulo)
    ax_amp.set_ylabel("Amplitude [g/N]")
    ax_fase.set_ylabel("Fase [graus]")
    ax_coe.set_ylabel("Coerência")
    ax_coe.set_xlabel("Frequência [Hz]")
    ax_coe.set_ylim(0.0, 1.05)
    ax_coe.set_xlim(0.0, float(dados_numericos["frequencia"][-1]))
    ax_amp.legend(loc="best", frameon=False, fontsize=8)
    for eixo in eixos:
        eixo.grid(True, which="both", color="0.88", linewidth=0.55)
        eixo.spines[["top", "right"]].set_visible(False)

    if arquivo_saida is not None:
        fig.savefig(
            arquivo_saida,
            dpi=DPI_MATPLOTLIB,
            bbox_inches="tight",
            facecolor="white",
        )
        print(f"Gráfico de FRF salvo em: {arquivo_saida}")
    if mostrar:
        plt.show()
    return fig


def plotar_frf(
    dados_numericos,
    dados_experimentais=None,
    mostrar=True,
    arquivo_saida=None,
    biblioteca=None,
):
    """Plota a FRF com Plotly ou Matplotlib."""
    biblioteca = _normalizar_biblioteca(biblioteca)
    if biblioteca == "plotly":
        return _plotar_frf_plotly(
            dados_numericos,
            dados_experimentais,
            mostrar,
            arquivo_saida,
        )
    return _plotar_frf_matplotlib(
        dados_numericos,
        dados_experimentais,
        mostrar,
        arquivo_saida,
    )


def comparar_malhas(
    incluir_suspensao=True,
    tipos_malha=None,
    massas_acelerometros=None,
    representacao_suspensao=REPRESENTACAO_SUSPENSAO,
    c_suspensao=AMORTECIMENTO_SUSPENSAO,
):
    """Compara propriedades e modos das malhas com a discretização completa."""
    if tipos_malha is None:
        tipos_malha = tuple(MALHAS_CONVERGENCIA)
    elif isinstance(tipos_malha, str):
        tipos_malha = (tipos_malha,)

    tipos_malha = tuple(
        dict.fromkeys(_normalizar_tipo_malha(tipo) for tipo in tipos_malha)
    )
    if "completa" not in tipos_malha:
        tipos_malha += ("completa",)

    rotores = {
        tipo: rotor_identificado_frf(
            incluir_suspensao=incluir_suspensao,
            incluir_amortecimento=False,
            tipo_malha=tipo,
            massas_acelerometros=massas_acelerometros,
            representacao_suspensao=representacao_suspensao,
            c_suspensao=c_suspensao,
        )
        for tipo in tipos_malha
    }
    resultados = {}
    for tipo, rotor in rotores.items():
        resultados[tipo] = {
            "rotor": rotor,
            "propriedades": propriedades_de_massa(rotor),
            "frequencias": frequencias_flexiveis_y(rotor),
            "numero_nos": len(rotor.nodes_pos),
            "numero_elementos_eixo": len(rotor.shaft_elements),
            "ndof": rotor.ndof,
            "nos_funcionais": rotor.parametros_modelo["nos_funcionais"].copy(),
        }

    frequencias_referencia = resultados["completa"]["frequencias"]
    ndof_referencia = resultados["completa"]["ndof"]
    for tipo in tipos_malha:
        erro = 100 * (resultados[tipo]["frequencias"] / frequencias_referencia - 1)
        resultados[tipo]["erro_frequencias_percentual"] = erro
        resultados[tipo]["erro_modal_maximo_percentual"] = float(np.max(np.abs(erro)))
        resultados[tipo]["reducao_gdl_percentual"] = float(
            100 * (1 - resultados[tipo]["ndof"] / ndof_referencia)
        )

    # Mantém compatibilidade com chamadas anteriores que consultavam somente
    # o erro da malha reduzida diretamente na raiz do dicionário.
    if "reduzida" in resultados:
        resultados["erro_frequencias_percentual"] = resultados["reduzida"][
            "erro_frequencias_percentual"
        ]
    return resultados


def imprimir_resumo(rotor):
    """Imprime as propriedades impostas e as frequências identificadas."""
    props = propriedades_de_massa(rotor)
    frequencias = frequencias_flexiveis(rotor, DIRECAO_FRF)
    erros = 100 * (frequencias / FREQUENCIAS_ALVO - 1)
    print("\nModelo do rotor")
    print(f"  Malha:              {rotor.parametros_modelo['tipo_malha']}")
    print(f"  Nós / GDL:          {len(rotor.nodes_pos)} / {rotor.ndof}")
    print(f"  Suspensão incluída: {rotor.parametros_modelo['incluir_suspensao']}")
    print(
        "  Representação:      " f"{rotor.parametros_modelo['representacao_suspensao']}"
    )
    print(
        "  Amort. por cordão:  " f"{rotor.parametros_modelo['c_suspensao']:.6f} N.s/m"
    )
    print("  E efetivo do eixo:  " f"{rotor.parametros_modelo['E_eixo'] / 1e9:.9f} GPa")
    print(f"  Massa total:        {props['massa']:.9f} kg")
    massas_sensores = rotor.parametros_modelo["massas_acelerometros"]
    print(
        "  Massas sensores:   "
        f"DE={massas_sensores['DE']:.6f} kg / "
        f"NDE={massas_sensores['NDE']:.6f} kg"
    )
    print(f"  Massa sem anel:     {rotor.parametros_modelo['massa_sem_anel']:.9f} kg")
    print(f"  Massa do anel:      {rotor.parametros_modelo['massa_anel']:.9f} kg")
    print(f"  Ip total:           {props['Ip']:.9f} kg.m²")
    print(f"  CG:                 {props['CG']:.9f} m")
    print(f"  It (diagnóstico):   {props['It']:.9f} kg.m²")
    print(f"  Modos em {DIRECAO_FRF}:          {frequencias} Hz")
    print(f"  Erros relativos:    {erros} %")


def executar():
    """Executa as opções selecionadas na seção de configurações."""
    if COMPARAR_COM_EXPERIMENTAL and not INCLUIR_SUSPENSAO:
        print(
            "Aviso: os dados experimentais foram obtidos com o rotor suspenso, "
            "mas INCLUIR_SUSPENSAO está definido como False."
        )
    massas_acelerometros = (
        {
            "DE": MASSA_ACELEROMETRO_DE,
            "NDE": MASSA_ACELEROMETRO_NDE,
        }
        if INCLUIR_MASSAS_ACELEROMETROS
        else None
    )
    rotor = rotor_identificado_frf(
        incluir_suspensao=INCLUIR_SUSPENSAO,
        incluir_amortecimento=INCLUIR_AMORTECIMENTO,
        tipo_malha=TIPO_MALHA,
        massas_acelerometros=massas_acelerometros,
        representacao_suspensao=REPRESENTACAO_SUSPENSAO,
        c_suspensao=AMORTECIMENTO_SUSPENSAO,
    )
    imprimir_resumo(rotor)
    if VALIDAR_REDUCAO_MALHA:
        validacao = comparar_malhas(
            INCLUIR_SUSPENSAO,
            massas_acelerometros=massas_acelerometros,
            representacao_suspensao=REPRESENTACAO_SUSPENSAO,
            c_suspensao=AMORTECIMENTO_SUSPENSAO,
        )
        print("\nValidação modal das malhas (referência: completa)")
        print(
            "  Malha          Nós  Elem.  GDL      f1 [Hz]      f2 [Hz]"
            "      f3 [Hz]   Erro máx. [%]"
        )
        for tipo in MALHAS_CONVERGENCIA:
            dados = validacao[tipo]
            f1, f2, f3 = dados["frequencias"]
            print(
                f"  {tipo:14s} {dados['numero_nos']:3d}  "
                f"{dados['numero_elementos_eixo']:5d}  {dados['ndof']:3d}  "
                f"{f1:11.5f}  {f2:11.5f}  {f3:11.5f}  "
                f"{dados['erro_modal_maximo_percentual']:13.5f}"
            )

        print("\nMapeamento dos nós em cada malha")
        print(
            "  Malha          Sensor DE  Atuador DE  CG  Disco  "
            "Atuador NDE  Sensor NDE"
        )
        for tipo in MALHAS_CONVERGENCIA:
            nos = validacao[tipo]["nos_funcionais"]
            print(
                f"  {tipo:14s} {nos['sensor_DE']:9d}  "
                f"{nos['atuador_DE']:10d}  {nos['CG']:2d}  "
                f"{nos['disco']:5d}  {nos['atuador_NDE']:11d}  "
                f"{nos['sensor_NDE']:10d}"
            )
    biblioteca = _normalizar_biblioteca(BIBLIOTECA_GRAFICA)

    if PLOTAR_ROTOR:
        plotar_rotores_malhas(
            tipos_malha=MALHAS_PARA_PLOTAR,
            incluir_suspensao=INCLUIR_SUSPENSAO,
            representacao_suspensao=REPRESENTACAO_SUSPENSAO,
            incluir_amortecimento=INCLUIR_AMORTECIMENTO,
            mostrar=MOSTRAR_GRAFICOS,
            mostrar_nos=MOSTRAR_NOS_ROTOR,
            salvar=SALVAR_GRAFICOS,
            biblioteca=biblioteca,
            caminho_base=ARQUIVO_ROTOR,
            rotor_atual=rotor,
            massas_acelerometros=massas_acelerometros,
            c_suspensao=AMORTECIMENTO_SUSPENSAO,
        )

    if CALCULAR_FRF_NUMERICA or COMPARAR_COM_EXPERIMENTAL or AVALIAR_ERRO_COMPLEXO_FRF:
        razoes_amortecimento = (
            RAZOES_AMORTECIMENTO_MODAL
            if (INCLUIR_AMORTECIMENTO and USAR_AMORTECIMENTO_MODAL_EXPERIMENTAL)
            else None
        )
        dados_numericos = frf_numerica_modal(
            rotor,
            lado=LADO_FRF,
            direcao=DIRECAO_FRF,
            razoes_amortecimento=razoes_amortecimento,
            polaridade=POLARIDADE_FRF_NUMERICA,
        )
        dados_experimentais = None
        if COMPARAR_COM_EXPERIMENTAL or AVALIAR_ERRO_COMPLEXO_FRF:
            dados_experimentais = carregar_frf_experimental(
                PASTA_DADOS,
                lado=LADO_FRF,
                experimentos=EXPERIMENTOS_FRF,
            )
        if AVALIAR_ERRO_COMPLEXO_FRF:
            avaliacao = avaliar_erro_complexo_frf(
                dados_numericos,
                dados_experimentais,
                coerencia_minima=COERENCIA_MINIMA_ERRO_FRF,
                frequencia_maxima=FREQUENCIA_MAX_FRF,
            )
            print(
                "\nErro da FRF complexa ponderado pela coerência "
                f"(mínima={COERENCIA_MINIMA_ERRO_FRF:.2f})"
            )
            for identificador in sorted(dados_experimentais):
                resultado = avaliacao[identificador]
                print(
                    f"  Ensaio {identificador}: "
                    f"{resultado['erro_complexo_percentual']:.2f}% "
                    f"({resultado['numero_pontos']} pontos)"
                )
            print("  Mediana: " f"{avaliacao['mediana_percentual']:.2f}%")
        caminho_frf = None
        if SALVAR_GRAFICOS:
            sufixo = "comparacao" if COMPARAR_COM_EXPERIMENTAL else "numerica"
            caminho_base = (
                Path(__file__).resolve().parent / f"frf_{sufixo}_{LADO_FRF.lower()}"
            )
            caminho_frf = _caminho_grafico(caminho_base, biblioteca)
        if CALCULAR_FRF_NUMERICA or COMPARAR_COM_EXPERIMENTAL:
            plotar_frf(
                dados_numericos,
                dados_experimentais=(
                    dados_experimentais if COMPARAR_COM_EXPERIMENTAL else None
                ),
                mostrar=MOSTRAR_GRAFICOS,
                arquivo_saida=caminho_frf,
                biblioteca=biblioteca,
            )
    return rotor


if __name__ == "__main__":
    rotor = executar()
