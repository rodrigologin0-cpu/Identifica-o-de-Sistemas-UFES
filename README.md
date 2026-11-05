# Simulação Buck-Boost — PI, LQR, LQI e observador de Luenberger

Estudo acadêmico em Python de um conversor CC-CC Buck-Boost ideal, usando modelo médio em condução contínua (CCM) e representação em espaço de estados.

## Estratégias comparadas

1. **PI** — realimentação da tensão medida.
2. **LQR** — realimentação da tensão e corrente medidas, sem ação integral.
3. **LQR + Luenberger** — realimentação dos estados estimados por observador.
4. **LQI** — LQR com estado integral, usando tensão e corrente medidas.
5. **LQI + Luenberger** — LQI com estados estimados.

O observador de Luenberger não é um controlador independente: ele estima estados que são utilizados pelo LQR ou LQI.

## Executar

Requer Python 3, NumPy, SciPy e Matplotlib.

```bash
python -m pip install -r requirements.txt
python buckboost_comparacao_5_controladores.py
```

A simulação apresenta métricas no terminal e gera gráficos comparativos e um CSV com resultados.

## Cenário e limitações

- Referência de saída: 24 V.
- Entrada nominal: 12 V.
- Perturbação de carga em 120 ms.
- Alteração da tensão de entrada em 280 ms.
- Ruído de medição adicionado ao teste.
- Modelo médio idealizado: não representa a comutação PWM individual, ripple, resistências parasitas ou outros efeitos não ideais.

**Interpretação:** estabilidade da malha não equivale necessariamente à eliminação do erro estacionário. O LQI inclui ação integral; o LQR clássico deste experimento não a inclui. A proximidade das respostas com e sem observador depende da qualidade do modelo e da estimação.

Os ganhos não foram otimizados para igualar tempos de resposta ou esforço de controle entre métodos. Os resultados são ilustrativos e não constituem validação para implementação em hardware.

## Arquivo principal

[`buckboost_comparacao_5_controladores.py`](buckboost_comparacao_5_controladores.py)

## Uso

Material de estudo acadêmico em teoria de controle, sistemas lineares e observabilidade.
