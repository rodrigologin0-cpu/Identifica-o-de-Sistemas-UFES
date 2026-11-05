#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SIMULACAO ACADEMICA | FONTE BUCK-BOOST IDEAL | CINCO CONFIGURACOES

Objetivo: comparar o seguimento de tensao, a rejeicao de perturbacoes,
          o esforco de controle e a estimacao dos estados.

Controladores:
    1. PI: apenas realimentacao da tensao medida.
    2. LQR: corrente e tensao medidas, sem acao integral.
    3. LQR + Luenberger: corrente e tensao estimadas pelo observador.
    4. LQI: LQR com estado integral, corrente e tensao medidas.
    5. LQI + Luenberger: LQI com estados estimados pelo observador.

Diferenca conceitual:
    PI, LQR e LQI SAO controladores.
    Luenberger e um OBSERVADOR, que fornece estados ao LQR/LQI.

Modelo: medio ideal, CCM, sem representar transistor, comutacao PWM,
        ripple, resistencias parasitas, atrasos ou nao idealidades.

Preserva todos os parametros, PI e LQR do primeiro experimento.
O LQI tem pesos novos e explicitamente identificados (Q_LQI/R_LQI).

Instalacao: py -3.13 -m pip install numpy scipy matplotlib
Execucao:   py -3.13 .\buckboost_comparacao_5_controladores.py

Graficos e CSV sao salvos na pasta deste arquivo.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # Gera PNG mesmo em terminais sem ambiente grafico.
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import solve_continuous_are
from scipy.signal import place_poles


# ==========================================================================
# 1. PARAMETROS ORIGINAIS (NAO ALTERADOS EM RELACAO AO PRIMEIRO ENSAIO)
# ==========================================================================
VIN_NOM = 12.0
VREF = 24.0  # Magnitude da tensao da Buck-Boost inversora [V]
R_NOM = 12.0
INDUTANCIA = 1.0e-3
CAPACITANCIA = 470.0e-6
T_FINAL = 0.45
DT = 50.0e-6
D_MIN, D_MAX = 0.05, 0.90
SEMENTE = 2026
DESVIO_RUIDO_V = 0.05  # [V]
DESVIO_RUIDO_I = 0.02  # [A], para os casos com sensor de corrente

# Mudancas na planta: identicas para os cinco casos.
T_CARGA = 0.12
R_APOS_CARGA = 9.0
T_ENTRADA = 0.28
VIN_APOS_ENTRADA = 11.0

# Ganhos PI mantidos.
KP_PI = 0.0025
KI_PI = 0.70

# Pesos LQR mantidos.
Q_LQR = np.diag([0.02, 1.0])
R_LQR = np.array([[1000.0]])

# Polos originais do observador mantidos [1/s].
POLOS_OBSERVADOR = [-350.0, -550.0]
OBSERVADOR_CONHECE_CARGA = True

# ==========================================================================
# 2. PARAMETROS NOVOS: PROJETO LQI COM ESTADO INTEGRAL
# ==========================================================================
# z(t) = integral(VREF - Vout) dt. A matriz aumentada inclui z.
# Q_LQI pondera [erro de corrente, erro de tensao, integral de tensao].
# O peso integral abaixo foi escolhido como ponto de partida para estudo.
# Nao foi imposto que PI e LQI tenham exatamente o mesmo tempo de resposta.
Q_LQI = np.diag([0.02, 1.0, 10000.0])
R_LQI = np.array([[1000.0]])

# O integrador de LQI usa a tensao MEDIDA em ambas as configuracoes.
# O estimador Luenberger fornece corrente e tensao para o ganho de estados.
# A integraçao utiliza anti-windup condicional, como no PI.

CASOS = (
    "PI",
    "LQR",
    "LQR + Luenberger",
    "LQI",
    "LQI + Luenberger",
)
PASTA_SAIDA = Path(__file__).resolve().parent


# ==========================================================================
# 3. MODELO DA BUCK-BOOST EM CONDUCAO CONTINUA (TENSAO EM MODULO)
# ==========================================================================
def dinamica(x: np.ndarray, duty: float, vin: float, resistencia: float) -> np.ndarray:
    """x=[iL (A), vC (V)]; equacoes medias nao lineares."""
    i_l, v_c = x
    return np.array([
        (duty * vin - (1.0 - duty) * v_c) / INDUTANCIA,
        ((1.0 - duty) * i_l - v_c / resistencia) / CAPACITANCIA,
    ])


def passo_rk4(x: np.ndarray, duty: float, vin: float, resistencia: float) -> np.ndarray:
    """Integra a planta por Runge-Kutta de quarta ordem em cada DT."""
    k1 = dinamica(x, duty, vin, resistencia)
    k2 = dinamica(x + 0.5 * DT * k1, duty, vin, resistencia)
    k3 = dinamica(x + 0.5 * DT * k2, duty, vin, resistencia)
    k4 = dinamica(x + DT * k3, duty, vin, resistencia)
    return x + DT * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0


def condicoes_externas(t: float) -> tuple[float, float]:
    """Mantem as mesmas duas perturbacoes do primeiro ensaio."""
    resistencia = R_NOM if t < T_CARGA else R_APOS_CARGA
    vin = VIN_NOM if t < T_ENTRADA else VIN_APOS_ENTRADA
    return vin, resistencia


def matrizes_linearizadas():
    """Jacobianas A, B, C no equilibrio nominal; entrada = duty cycle."""
    d0 = VREF / (VIN_NOM + VREF)
    i0 = VREF / (R_NOM * (1.0 - d0))
    x0 = np.array([i0, VREF])
    A = np.array([
        [0.0, -(1.0 - d0) / INDUTANCIA],
        [(1.0 - d0) / CAPACITANCIA, -1.0 / (R_NOM * CAPACITANCIA)],
    ])
    B = np.array([[(VIN_NOM + VREF) / INDUTANCIA], [-i0 / CAPACITANCIA]])
    C = np.array([[0.0, 1.0]])
    return A, B, C, x0, d0


def sintetizar_ganhos(A: np.ndarray, B: np.ndarray, C: np.ndarray):
    """Sintetiza K_LQR, K_LQI e L_observador via ARE e alocacao de polos."""
    # Controle LQR: delta_d = -K_LQR * delta_x.
    P_lqr = solve_continuous_are(A, B, Q_LQR, R_LQR)
    K_lqr = np.linalg.solve(R_LQR, B.T @ P_lqr).ravel()

    # z_ponto = VREF - vC, logo a linha adicional e [-C, 0].
    A_aug = np.block([
        [A, np.zeros((2, 1))],
        [-C, np.zeros((1, 1))],
    ])
    B_aug = np.vstack((B, np.zeros((1, 1))))
    # delta_d = -K_LQI * [delta_i, delta_v, z].
    # Com z_ponto = VREF - vC, K_LQI[2] e NEGATIVO.
    # Assim, a parcela integral efetiva no duty e POSITIVA.
    P_lqi = solve_continuous_are(A_aug, B_aug, Q_LQI, R_LQI)
    K_lqi = np.linalg.solve(R_LQI, B_aug.T @ P_lqi).ravel()

    # Ganho do observador (mantido igual ao da primeira versao).
    L_obs = place_poles(A.T, C.T, POLOS_OBSERVADOR).gain_matrix.T.ravel()
    return K_lqr, K_lqi, L_obs, A_aug, B_aug


# ==========================================================================
# 4. SIMULACAO DAS CINCO CONFIGURACOES COM RUIDO IDENTICO
# ==========================================================================
def simular(tipo: str, K_lqr: np.ndarray, K_lqi: np.ndarray, L_obs: np.ndarray,
            x0: np.ndarray, d0: float, ruido_v: np.ndarray, ruido_i: np.ndarray) -> dict:
    """Uma planta por caso; mesma condicao inicial e mesma realizacao de ruido."""
    n = len(ruido_v)
    tempo = np.arange(n) * DT
    x = x0.copy()
    xhat = x0 + np.array([-1.0, 0.5])  # Estimativa inicialmente incorreta.
    integral_erro = 0.0

    tensao = np.empty(n)
    corrente = np.empty(n)
    tensao_medida = np.empty(n)
    duty_hist = np.empty(n)
    corrente_est = np.full(n, np.nan)
    tensao_est = np.full(n, np.nan)
    integral_hist = np.full(n, np.nan)
    resíduo_hist = np.full(n, np.nan)

    usa_observador = tipo in ("LQR + Luenberger", "LQI + Luenberger")
    usa_integral = tipo in ("PI", "LQI", "LQI + Luenberger")

    for k, t in enumerate(tempo):
        vin, resistencia = condicoes_externas(t)
        y_medido = x[1] + ruido_v[k]
        tensao[k], corrente[k], tensao_medida[k] = x[1], x[0], y_medido
        erro = VREF - y_medido

        if tipo == "PI":
            # PI original sem qualquer alteracao de ganhos.
            comando = d0 + KP_PI * erro + KI_PI * integral_erro
        elif tipo == "LQR":
            # LQR original: ambos os estados medidos, sujeitos a ruido.
            x_medido = np.array([x[0] + ruido_i[k], y_medido])
            comando = d0 - K_lqr @ (x_medido - x0)
        elif tipo == "LQR + Luenberger":
            # LQR original: estado estimado por observador de Luenberger.
            comando = d0 - K_lqr @ (xhat - x0)
        elif tipo == "LQI":
            # LQI com sensor de corrente e tensao; z e erro integral da tensao.
            x_medido = np.array([x[0] + ruido_i[k], y_medido])
            vetor_aumentado = np.array([*(x_medido - x0), integral_erro])
            comando = d0 - K_lqi @ vetor_aumentado
        elif tipo == "LQI + Luenberger":
            # LQI: corrente estimada. Integral da tensao MEDIDA para eliminar offset.
            vetor_aumentado = np.array([*(xhat - x0), integral_erro])
            comando = d0 - K_lqi @ vetor_aumentado
        else:
            raise ValueError(f"Configuracao desconhecida: {tipo}")

        duty = float(np.clip(comando, D_MIN, D_MAX))

        # Anti-windup condicional preserva logica original do PI.
        # O sinal de K_lqi[2] e negativo; erro > 0 tende a aumentar duty.
        if usa_integral:
            saturacao_alta = comando >= D_MAX and erro > 0
            saturacao_baixa = comando <= D_MIN and erro < 0
            if not (saturacao_alta or saturacao_baixa):
                integral_erro += erro * DT
            integral_hist[k] = integral_erro

        if usa_observador:
            corrente_est[k], tensao_est[k] = xhat
            R_modelo = resistencia if OBSERVADOR_CONHECE_CARGA else R_NOM
            residuo = y_medido - xhat[1]
            resíduo_hist[k] = residuo
            # Mesma integracao por Euler e polos da primeira versao.
            xhat = xhat + DT * (dinamica(xhat, duty, vin, R_modelo) + L_obs * residuo)

        duty_hist[k] = duty
        x = passo_rk4(x, duty, vin, resistencia)
        if not np.all(np.isfinite(x)):
            raise RuntimeError(f"Divergencia de {tipo} em t={t:.6f} s")

    return dict(t=tempo, v=tensao, i=corrente, v_medido=tensao_medida,
                duty=duty_hist, i_hat=corrente_est, v_hat=tensao_est,
                integral=integral_hist, residuo=resíduo_hist)


# ==========================================================================
# 5. METRICAS OBJETIVAS E EXPORTACAO EM CSV
# ==========================================================================
def metricas(sim: dict) -> dict:
    """Erro entre 120 e 450 ms; erro final na ultima janela de 25 ms."""
    t, v, d = sim["t"], sim["v"], sim["duty"]
    mascara = t >= T_CARGA
    erro = VREF - v[mascara]
    janela_final = t >= (T_FINAL - 0.025)
    mascara_vin = (t >= T_ENTRADA) & (t < T_FINAL)
    # Desvio padrao de duty na janela final (indicador simples de ruido).
    return {
        "RMSE [V]": float(np.sqrt(np.mean(erro**2))),
        "MAE [V]": float(np.mean(np.abs(erro))),
        "Maior erro [V]": float(np.max(np.abs(erro))),
        "Erro medio final [V]": float(np.mean(VREF - v[janela_final])),
        "RMSE apos Vin [V]": float(np.sqrt(np.mean((VREF - v[mascara_vin])**2))),
        "Duty saturado [%]": float(100 * np.mean((d <= D_MIN + 1e-8) | (d >= D_MAX - 1e-8))),
        "Duty medio final": float(np.mean(d[janela_final])),
        "Desvio duty final": float(np.std(d[janela_final])),
    }


def salvar_metricas(resultados: dict) -> Path:
    """Gera CSV com ponto-e-virgula para importacao no Excel em PT-BR."""
    arquivo = PASTA_SAIDA / "buckboost_5_controladores_metricas.csv"
    with arquivo.open("w", newline="", encoding="utf-8-sig") as f:
        escritor = csv.writer(f, delimiter=";")
        escritor.writerow(["Controlador", *metricas(next(iter(resultados.values()))).keys()])
        for nome, sim in resultados.items():
            m = metricas(sim)
            escritor.writerow([nome, *(f"{v:.6f}".replace(".", ",") for v in m.values())])
    return arquivo


# ==========================================================================
# 6. GRAFICOS
# ==========================================================================
def grafico_comparativo(resultados: dict) -> Path:
    """Tres paines: tensao, corrente do indutor e duty cycle."""
    fig, eixos = plt.subplots(3, 1, figsize=(13.3, 9.6), sharex=True, layout="constrained")
    cores = {
        "PI": "tab:blue", "LQR": "tab:orange", "LQR + Luenberger": "tab:green",
        "LQI": "tab:red", "LQI + Luenberger": "tab:purple",
    }
    estilos = {
        "PI": "-", "LQR": "-", "LQR + Luenberger": "--",
        "LQI": "-", "LQI + Luenberger": "--",
    }
    for nome, sim in resultados.items():
        t_ms = 1000.0 * sim["t"]
        for ax, campo in zip(eixos, ("v", "i", "duty")):
            ax.plot(t_ms, sim[campo], label=nome, color=cores[nome], ls=estilos[nome], lw=1.45,
                    alpha=0.95 if "Luenberger" not in nome else 0.90)
    eixos[0].axhline(VREF, color="black", ls=":", lw=1.4, label="Referencia 24 V")
    for ax in eixos:
        ax.axvline(1000.0 * T_CARGA, ls=":", color="gray", lw=1.0)
        ax.axvline(1000.0 * T_ENTRADA, ls=":", color="gray", lw=1.0)
        ax.grid(True, alpha=0.25)
    eixos[0].set_ylabel("Tensao de saida [V]")
    eixos[1].set_ylabel("Corrente no indutor [A]")
    eixos[2].set_ylabel("Duty cycle")
    eixos[2].set_xlabel("Tempo [ms]")
    eixos[0].legend(loc="best", fontsize=9, ncol=3)
    eixos[0].set_title("Buck-Boost | PI, LQR e LQI | com/sem observador")
    arquivo = PASTA_SAIDA / "buckboost_5_controladores.png"
    fig.savefig(arquivo, dpi=160)
    plt.close(fig)
    return arquivo


def grafico_observadores(resultados: dict) -> Path:
    """Compara correntes reais/estimadas nos dois casos com Luenberger."""
    fig, eixos = plt.subplots(3, 1, figsize=(13.3, 9.4), sharex=True, layout="constrained")
    cores = {"LQR + Luenberger": "tab:green", "LQI + Luenberger": "tab:purple"}
    for nome, cor in cores.items():
        sim = resultados[nome]
        t_ms = 1000.0 * sim["t"]
        eixos[0].plot(t_ms, sim["i"], color=cor, lw=1.5, label=f"Corrente real — {nome}")
        eixos[0].plot(t_ms, sim["i_hat"], color=cor, ls="--", lw=1.1, label=f"Estimada — {nome}")
        eixos[1].plot(t_ms, sim["i"] - sim["i_hat"], color=cor, lw=1.1, label=nome)
        eixos[2].plot(t_ms, sim["residuo"], color=cor, lw=0.7, alpha=0.7, label=nome)
    eixos[0].set_ylabel("Corrente [A]")
    eixos[1].set_ylabel("Erro i real - estimada [A]")
    eixos[2].set_ylabel("Inovacao v medida - estimada [V]")
    eixos[2].set_xlabel("Tempo [ms]")
    eixos[1].axhline(0, color="black", lw=0.8, ls="--")
    for ax in eixos:
        ax.axvline(1000.0 * T_CARGA, ls=":", color="gray", lw=1.0)
        ax.axvline(1000.0 * T_ENTRADA, ls=":", color="gray", lw=1.0)
        ax.grid(True, alpha=0.25)
        ax.legend(loc="best", fontsize=9, ncol=2)
    fig.suptitle("Observador de Luenberger | corrente, erro e ruido de medicao")
    arquivo = PASTA_SAIDA / "buckboost_5_observadores.png"
    fig.savefig(arquivo, dpi=160)
    plt.close(fig)
    return arquivo


# ==========================================================================
# 7. EXECUCAO E VALIDACAO
# ==========================================================================
def main() -> None:
    if not (0.0 < T_CARGA < T_ENTRADA < T_FINAL):
        raise ValueError("Tempos de perturbacoes invalidos")
    if DT <= 0 or not (0 < D_MIN < D_MAX < 1):
        raise ValueError("Passo DT ou limites de duty invalidos")

    A, B, C, x0, d0 = matrizes_linearizadas()
    O = np.vstack((C, C @ A))
    rank_obs = np.linalg.matrix_rank(O)
    if rank_obs != 2:
        raise ValueError("O sistema nao e observavel a partir de vC no ponto nominal")
    K_lqr, K_lqi, L_obs, A_aug, B_aug = sintetizar_ganhos(A, B, C)

    polos_lqr = np.linalg.eigvals(A - B @ K_lqr[None, :])
    polos_lqi = np.linalg.eigvals(A_aug - B_aug @ K_lqi[None, :])
    polos_obs = np.linalg.eigvals(A - L_obs[:, None] @ C)
    for nome, polos in (("LQR", polos_lqr), ("LQI", polos_lqi), ("Luenberger", polos_obs)):
        if not np.all(np.real(polos) < 0):
            raise RuntimeError(f"Autovalores instaveis para {nome}: {polos}")

    print("\n=== SIMULACAO BUCK-BOOST | CINCO CONFIGURACOES ===")
    print("Modelo medio ideal (CCM) | output em modulo | dt=50 us")
    print(f"Equilibrio nominal: Vin={VIN_NOM:.1f} V, Vref={VREF:.1f} V, R={R_NOM:.1f} ohms")
    print(f"iL0={x0[0]:.4f} A | duty0={d0:.5f}")
    print("A =\n", A, "\nB =\n", B, "\nC =\n", C)
    print(f"Posto da observabilidade = {rank_obs} de 2")
    print(f"K_LQR = {np.round(K_lqr, 7)}")
    print(f"K_LQI [K_iL, K_vC, K_z] = {np.round(K_lqi, 7)}")
    print(f"L_observador = {np.round(L_obs, 5)}")
    print(f"Polos LQR: {np.round(polos_lqr, 2)}")
    print(f"Polos LQI: {np.round(polos_lqi, 2)}")
    print(f"Polos observador: {np.round(polos_obs, 2)}")
    print(f"Q_LQI = diag({Q_LQI[0,0]:g}, {Q_LQI[1,1]:g}, {Q_LQI[2,2]:g})")
    print(f"R_LQI = {R_LQI[0,0]:g}")
    print(f"Ruido de tensao: {DESVIO_RUIDO_V:g} V (1 sigma)")
    print(f"Observador conhece R atualizado: {OBSERVADOR_CONHECE_CARGA}")
    print(f"Perturbacoes: R={R_NOM:g}->{R_APOS_CARGA:g} ohms em {T_CARGA*1000:g} ms; "
          f"Vin={VIN_NOM:g}->{VIN_APOS_ENTRADA:g} V em {T_ENTRADA*1000:g} ms")

    n = round(T_FINAL / DT) + 1
    rng = np.random.default_rng(SEMENTE)
    ruido_v = rng.normal(0.0, DESVIO_RUIDO_V, size=n)
    ruido_i = rng.normal(0.0, DESVIO_RUIDO_I, size=n)
    resultados = {nome: simular(nome, K_lqr, K_lqi, L_obs, x0, d0, ruido_v, ruido_i)
                  for nome in CASOS}

    print("\n=== METRICAS APOS PRIMEIRA PERTURBACAO (120 a 450 ms) ===")
    cabecalho = f"{'Configuracao':<24} {'RMSE V':>9} {'MAE V':>9} {'Max V':>9} {'Erro final V':>13} {'Ruido duty':>12}"
    print(cabecalho)
    print("-" * len(cabecalho))
    for nome, sim in resultados.items():
        m = metricas(sim)
        print(f"{nome:<24} {m['RMSE [V]']:>9.4f} {m['MAE [V]']:>9.4f} "
              f"{m['Maior erro [V]']:>9.4f} {m['Erro medio final [V]']:>13.4f} "
              f"{m['Desvio duty final']:>12.6f}")

    print("\n=== ERRO DA ESTIMATIVA DE CORRENTE (APOS 40 ms) ===")
    for nome in ("LQR + Luenberger", "LQI + Luenberger"):
        sim = resultados[nome]
        mascara = sim["t"] >= 0.04
        e_i = sim["i"][mascara] - sim["i_hat"][mascara]
        print(f"{nome:<24} MAE={np.mean(np.abs(e_i)):.5f} A, "
              f"RMSE={np.sqrt(np.mean(e_i**2)):.5f} A")

    arquivo_png = grafico_comparativo(resultados)
    arquivo_obs = grafico_observadores(resultados)
    arquivo_csv = salvar_metricas(resultados)
    print("\nArquivos gerados:")
    for arquivo in (arquivo_png, arquivo_obs, arquivo_csv):
        print(f"  {arquivo}")
    print("\nIMPORTANTE: LQI incorpora integrador; LQR puro nao.")
    print("IMPORTANTE: ganhos nao foram otimizados para igualar tempos de resposta.")
    print("IMPORTANTE: convergencia de tensao != estabilidade do estado.")


if __name__ == "__main__":
    main()