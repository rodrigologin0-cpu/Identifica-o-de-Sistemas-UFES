# Mosturação 2026 — PID e alternativas para múltiplos pontos de operação

Este projeto reúne código Python e **dados sintéticos de simulação** usados na comparação de quatro estratégias: PID fixo, PID Gain Scheduling por faixas, Fuzzy-PID Gain Scheduling e Fuzzy puro incremental.

## Reproduzir os resultados

1. Instale Python 3.11+.
2. Execute `pip install -r requirements.txt`.
3. Execute `python mosturacao_2026_definitivo.py`.

Os parâmetros utilizados na publicação estão em `mosturacao_2026_definitivo_resultados/parametros_definitivos.json`. O comando acima lê essas sintonias e regenera os arquivos CSV e os gráficos originais com Matplotlib. Para reexecutar a busca: `python mosturacao_2026_definitivo.py --otimizar --maxiter 36 --popsize 6`.

## Datasets e plots

- `metricas_por_patamar.csv`: indicadores individuais, incluindo sobressinal e acomodação.
- `resumo_comparativo.csv` e `resumo_por_receita.csv`: resumos por controlador e receita.
- `series_temporais_nominais.csv`: sinais sintéticos de SP, PV, MV e ganhos nas duas receitas e quatro controladores (quando incluído).
- `*_curva.png`, `*_erro.png`: plots Matplotlib gerados a partir da simulação, ao executar o script.

A primeira receita segue 62 / 70,5 / 78 °C e a segunda 44 / 52 / 63 / 70,5 / 78 °C. Os quatro modelos FOPDT foram extraídos da identificação do estudo de 2018.

**Limitações:** associação dos modelos GP1–GP4 às temperaturas, interpolação, escala de MV e duração dos patamares são hipóteses de 2026. Não há validação física destas sintonias; os CSVs publicados são simulados, não os dados brutos experimentais de 2018. O modelo não representa toda a inércia térmica da panela. **Não usar os parâmetros diretamente em equipamento real.**

Artigo: [texto para publicação no LinkedIn](ARTIGO_LINKEDIN_MOSTURACAO_2026.md). Mais detalhes: [documentação técnica](mosturacao_2026_definitivo_resultados/README.md).

Autor do estudo original: R. S. Lima, IFES (2018). Código e análise computacional: Versão 2026.4.
