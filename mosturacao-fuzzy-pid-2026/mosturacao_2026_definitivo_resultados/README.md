# Mosturação 2026.4 — versão final computacional

## Objetivo
Comparar PID fixo, PID Gain Scheduling por faixas, Fuzzy-PID Gain Scheduling e Fuzzy puro incremental nas **curvas completas** das receitas 01 (62, 70,5, 78 °C) e 02 (44, 52, 63, 70,5, 78 °C) do trabalho IFES de 2018.

## Como executar no Windows (Python 3.11+ recomendado)

Instale `numpy`, `scipy` e `matplotlib`; `numba` é recomendado para desempenho durante a otimização.

```powershell
python -m pip install numpy scipy matplotlib numba
python mosturacao_2026_definitivo.py
```

O comando executa e exporta os gráficos usando as sintonias já otimizadas no arquivo `mosturacao_2026_definitivo_resultados/parametros_definitivos.json`.

Para repetir **exatamente a configuração de busca** utilizada nesta entrega:

```powershell
python mosturacao_2026_definitivo.py --otimizar --maxiter 36 --popsize 6
```

A busca pode produzir resultados ligeiramente diferentes conforme as versões de NumPy/SciPy/Numba. Para obter os mesmos números do relatório, utilize o JSON fornecido e a primeira forma de execução.

## Desenho experimental

- Dois perfis térmicos: Receita 01 (3 patamares) e Receita 02 (5 patamares).
- Primeira etapa: 1.800 s; seguintes: 1.200 s. Tempos de 2026, **não extraídos dos testes originais**.
- Temperatura inicial: 26 °C.
- Banda: ±0,2 °C e permanência mínima de 600 s consecutivos.
- Sobressinal: penalidade prioritária, com **tolerância numérica de 0,02 °C na escala hipotética**. Não é garantia física de sobresinal zero.
- Amostragem: 1 s.
- MV assumida: 0 a 100%.
- Treinamento: 5 cenários por receita = 40 patamares por controlador.
- Validação: 4 cenários adicionais por receita = 32 patamares por controlador, nunca acessados no ajuste.
- Comparação usa os mesmos cenários para todos os quatro controladores.
- A otimização utiliza diferencial evolutivo com a mesma quantidade de **gerações** e tamanho relativo de população. Entretanto, o número de avaliações difere devido à quantidade de parâmetros de cada controlador; não é uma comparação de custo computacional idêntico.

## Limitações fundamentais

1. Os valores GP1–GP4 foram publicados em 2018; a **associação dos GP às temperaturas** (centros 42, 52, 63 e 70,5 °C) foi suposta em 2026. A interpolação dinâmica entre eles também é uma hipótese.
2. O ganho estático foi interpretado como °C por ponto percentual de MV para construir o modelo; isso **não foi validado** contra os registros brutos. Unidades apresentadas nos gráficos dependem dessa convenção.
3. O modelo não representa energia térmica armazenada no metal da panela e na resistência, evaporação, variação de composição do mosto, atraso e ruído do sensor, nem bomba de recirculação explicitamente.
4. As sintonias e regras Fuzzy são **implementações 2026**, não reprodução fiel do controlador embarcado em 2018. A Fuzzy-PID usa consequentes numéricos ponderados (singleton), não a defuzzificação Mamdani original.
5. As condições reservadas para teste são variações do **mesmo modelo hipotético**, não ensaios físicos independentes.
6. Não é uma recomendação de parâmetros para conectar a uma resistência real sem identificação, proteções de segurança e validação experimental.

## Conteúdo

- `mosturacao_2026_definitivo.py`: código completo e autocontido da simulação/otimização.
- `parametros_definitivos.json`: parâmetros finais e configuração numérica da busca.
- `resumo_comparativo.csv`: resultado global de 40 patamares de treino + 32 de validação.
- `resumo_por_receita.csv`: resumo por receita e controlador.
- `metricas_por_patamar.csv`: indicadores individuais.
- `series_temporais_nominais.csv`: SP, PV, MV e ganhos por amostra (cenário nominal).
- `receita_*_curva.png`: sobreposição de SP, PV e MV para cada receita e cenários selecionados.
- `receita_*_erro.png`: erros em relação ao setpoint.
- `sintese_tempo_assentamento.png` / `sintese_validacao.png`: gráficos-resumo para a publicação.
- `ARTIGO_LINKEDIN_MOSTURACAO_2026.md`: texto pronto para editar e publicar no LinkedIn.

## Interpretação dos indicadores

- `Ts_s` só é computado para o patamar se a temperatura permaneceu dentro da banda até o fim do patamar, com pelo menos 600 s disponíveis.
- A mediana de Ts considera apenas os patamares acomodados; se uma estratégia não estabilizou alguns patamares, não comparar a mediana isoladamente.
- `atingiu_sem_ultrapassagem` requer acomodação e sobressinal ≤ 0,02 na escala hipotética.
- A duração do patamar é fixa para todos os controladores e não é adaptada conforme a resposta. Se um falha, o SP muda na hora programada; isso garante uma comparação justa da mesma receita.

Referência histórica: R. S. Lima, monografia/artigo IFES 2018, *Controle de Temperatura Aplicado ao Processo de Fabricação de Cerveja Artesanal com Uso de PID Fuzzy Gain Scheduling*.