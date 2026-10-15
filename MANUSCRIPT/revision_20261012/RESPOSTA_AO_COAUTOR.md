# Resposta técnica à revisão do coautor

Data: 12/10/2026. Manuscrito: `manuscript3_EAAI.tex`. Suplemento: `supplementary_information.tex`.

## Resultado e limite desta revisão

A versão foi substancialmente revisada. Os comentários foram tratados individualmente com correção, nova análise ou delimitação explícita. **Isso não significa que todos os experimentos sugeridos tenham sido executados ou que o artigo esteja incondicionalmente pronto para submissão.** Não foram executadas novas simulações MCNP/LS-DYNA, otimizações comparativas completas, CV aninhada ou ensaios físicos. Não foi feito depósito público nem presumida liberação institucional.

O PDF anterior não continha as últimas alterações do abstract. O texto-fonte atualizado pelo autor foi usado como ponto de partida. Título, autoria, afiliações, financiamento e CRediT foram preservados. A abertura solicitada para o abstract foi mantida. Uma cópia anterior está em `manuscript_before_revision.tex`.

## Resposta aos 25 itens da tabela C

### 1. Modelos incorretos na figura — corrigido

As figuras antigas não são usadas na nova compilação. O fluxograma identifica Ridge para dose e XGBoost para penetração. A figura de desempenho mostra o melhor aprendiz dentro de cada representação, diretamente dos CSVs congelados. Ver Métodos 2.1 e Resultados 3.1.

### 2. Duas versus três representações, S23 e seleção — corrigido com limitação metodológica declarada

RAW, PHYSICS e PHYSICS-S23 estão diferenciadas. A tabela de features enumera sete valores de resistência, sete produtos resistência–espessura, um total, um centroide, seis saltos e uma complacência de cisalhamento, totalizando 23 acréscimos. S23 não é a dimensão total do pipeline.

A seleção de representação/hiperparâmetros reutilizou folds de desenvolvimento. **Não foi CV aninhada.** Isso foi escrito expressamente em Métodos e Limitações. Nenhum experimento aninhado foi inventado. A evidência de teste continua condicionada ao modelo e à partição congelados.

### 3. Faixas de erro relativo MCNP — corrigido

Recalculadas dos dados finais: DoE 1,12–1,53%, validação 1,14–1,62%. A faixa 1,04–1,47% foi removida. `evidence_audit.json` identifica as fontes.

### 4. Pre-registered versus pre-specified — corrigido

Usado “pre-specified”. O manifesto local registra 2026-10-03 10:30:06 UTC, regra de seleção e SHA256 da frente fonte `610884374c37c9ddf636c6765f15f9b62d33023633b9747191c8e4d2ac2b1230`. Isso é proveniência local, não registro externo ou prova de preregistro público.

### 5. Restrição não ativa — nova análise e reenquadramento, ablação completa pendente

No DoE, previsão máxima 38,315 mm. Na união retida, 21,568 mm. Nenhum ponto desses conjuntos cai em (56,0666, 60] mm. A tabela de sensibilidade mostra contagens sem/com margem para 15, 20, 30 e 60 mm.

A 20 mm, a margem elimina 195 candidatos do DoE e 47 da união que seriam aceitos pela previsão pontual. **É filtragem dos conjuntos existentes, não reotimização.** O artigo não afirma que a margem melhorou segurança/Pareto a 60 mm. Reotimizações pareadas e simulações próximas ao novo limite continuam sendo novos experimentos.

### 6. Margem constante — reformulado

O artigo apresenta a equivalência exata `p_pred + 3,933441 <= 60` e `p_pred <= 56,066559`. É um back-off aditivo calibrado, não incerteza dependente do projeto. A calibração bilateral usada como limite superior é declarada.

Foi identificado ainda um detalhe de implementação: o quantil NumPy “higher” seleciona a ordem 287, não 286, para n=300. Os valores congelados são ligeiramente mais conservadores e foram preservados. Os valores de ambas as ordens constam no artigo e em `conformal_implementation_audit.csv`. Não há correção silenciosa pós-validação.

CQR, scores normalizados e ponderação sob shift são discutidos como alternativas não executadas. Trocar agora o método e avaliar nos mesmos 30 casos não produziria uma nova validação independente.

### 7. Viés da dose — nova análise concluída

Erro = solver − previsão. Média +0,08023 mSv, 28/30 erros positivos, MAE 0,08832 mSv e RMSE 0,10297 mSv. MAE é 1,913 vezes o do teste. Desvio MC médio 0,08314 mSv. Cinco casos têm |erro/σ_MC| > 1,96. Essa padronização não inclui incerteza do surrogate e não é tratada como teste z formal.

Configurações próximas dos limites de espessura representam 35,35% do DoE, 88,43% da união e 83,33% dos 30 casos. Só oito permutações dos 30 foram vistas no DoE. Esses indicadores sustentam uma discussão de deslocamento de distribuição, sem atribuir causa única ao viés.

### 8. Hipervolume e orçamento — corrigido, comparação igual-orçamento pendente

O resultado principal passa a ser o ganho médio por execução: 0,685769 → 0,724699, +5,68%. A união 0,7446/+8,6% fica secundária. O texto informa que o orçamento do DoE não é equivalente e não reivindica superioridade sobre busca aleatória, NSGA-III ou MOEA/D. Um benchmark equivalente ainda requer novas execuções de otimização.

### 9. Dominância em alta fidelidade — nova análise concluída

Usados dose MCNP e os mesmos cálculos determinísticos de massa/custo para ambos os conjuntos. Dos 1.996 casos elegíveis, 295 são não dominados em três objetivos. Nenhum dos 30 candidatos é dominado pelo DoE. Vinte e nove permanecem não dominados na união DoE+30. Os 30 dominam, em conjunto, 73 dos 295 membros da frente inicial.

A figura explicita que a projeção dose–massa não é o teste completo de dominância em três objetivos. As comparações usam estimativas pontuais MCNP, não dominância probabilística certificada. Ver `all_30_designs_and_HF_dominance.csv`.

### 10. Portabilidade — corrigido

Removido o enquadramento como blindagem pessoal/vestível ou comprovadamente carregável à mão. Massa 95,15–256,21 kg, altura 70 cm, raios 35–41 cm e áreas curvas 0,770–0,902 m² são explícitos. Interpretação: conceito de barreira substancial que exigiria apoio/manuseio assistido. Nenhum sistema de rodas ou suporte foi declarado como já projetado.

### 11. Custo — corrigido

Equações de massa e custo foram incluídas com sete coeficientes. O custo é índice comparativo congelado por kg, não preço de mercado validado. Misturas e precisão das densidades são distinguidas do arredondamento dos inputs. Não foram inventadas cotações comerciais.

### 12. Física, malha, FSP e fonte — especificação incluída, validação física pendente

Tabela nova extraída dos inputs: atividade 3,7 TBq, espectro Co-60, fonte cilíndrica, distâncias, detector, MODE P, DE/DF, FM=7400, 300 s, 10^8 histórias. FSP cúbico de 1,1 g a 650 m/s, 232,375 J, gap 10 mm, malha, constitutivos, erosão, contato, horizonte e timestep.

Não foi presumido que esse cubo equivale a um FSP normativo. A tabela identifica o contato SOFT=2 da Fase 5 sem afirmar que todo o histórico de recuperação foi idêntico. **Não há novo estudo de convergência de malha nem benchmark experimental.** Inputs completos e validação das propriedades constitutivas continuam importantes para eventual qualificação física.

### 13. “Physical authority” — corrigido

Substituído por modelo de referência numérica. Comparação surrogate–solver não é certificação experimental nem proteção garantida em campo.

### 14. Dez diagnósticos MCNP — auditoria e sensibilidade concluídas

Lidos os 2.000 outputs, com dez flags extraídas em todos. Falhas: check 1 (mean behavior), 11 casos; check 9 (FOM behavior), 36 casos; demais oito, zero. São 43 casos distintos, quatro falham ambos. O inventário inicial de mensagens de aviso não substitui essa leitura da tabela de dez testes.

Sensibilidade: retirados 32 casos afetados apenas do treino, hiperparâmetros Ridge congelados, teste original inalterado. RMSE 0,05730 → 0,05753 mSv, MAE 0,04616 → 0,04637. O modelo oficial não foi substituído. Evidência: `MCNP_all_ten_checks.csv` e `extra_diagnostics.json`.

### 15. Quatro exclusões LS-DYNA — descritas, imputação não executada

IDs 733, 832, 1277 e 1567. Desvios máximos registrados 30,65%, 13,39%, 13,26% e 6,98%, acima do gate de 5%. Todos possuem Al externo com espessuras 2,61–2,89 mm. Portanto, o risco de missingness informativa foi explicitado, não descartado com base nas previsões.

Os registros completos e as configurações foram preservados no suplemento/CSV. Imputar 60 mm como se fosse resultado simulado seria incorreto. Um experimento contrafactual de treinamento pode ser planejado separadamente, mas não foi executado nem misturado aos dados reais.

### 16. Horizonte temporal — limites declarados, extensão ampla pendente

9016 e 9026 a 0,60 ms produzem 3,961 e 4,042 mm. O caso 9016 ainda avança. Isso foi mantido explicitamente. Dois casos não demonstram estabilização universal. O artigo define penetração como máximo no horizonte observado e não alega parada definitiva em todos os casos.

### 17. RAW e baselines físicos — novas comparações concluídas

Melhor RAW dose: XGBoost, RMSE CV 0,32536. PHYSICS/Ridge: 0,05880, redução 81,9%. Melhor RAW penetração: HistGB, 2,94744. PHYSICS-S23/XGBoost: 2,31533, redução 21,4%.

Baselines pós-revisão exclusivamente no treino: Beer–Lambert log-linear 0,17139 mSv, transmissão exponencial afim 0,17062 mSv, RAW Ridge/log-dose 0,71014 mSv na melhor penalização inspecionada. Nenhum desses é um transporte completo com buildup. O nome e alcance estão delimitados, e nenhuma seleção oficial foi alterada.

### 18. XGBoost versus HistGB — pares reportados e linguagem corrigida

Diferença média pareada de RMSE −0,01743 mm. XGBoost perde em três dos cinco folds e ganha em dois. O texto não reivindica significância nem equivalência formal. Mantém a escolha operacional congelada. IC bootstrap do teste foram acrescentados, condicionais ao modelo. Não foi fabricado teste independente de seleção nem criado ensemble pós-hoc.

### 19. NSGA-II, métricas e dispersão — corrigido

Operadores OX/SBX, probabilidades, mutações, projeção, quantização, dez seeds, orçamento e checkpoints estão especificados. Figura de convergência incluída. SD HV=0,00248931 e IQR=0,00251764: não são iguais, apenas arredondam a 0,0025. S7 fornece os espaços de busca de hiperparâmetros extraídos do código.

### 20. Eficiência — quantificada apenas onde medida

Inferência por caso em lote de 300: Ridge/PHYSICS ~0,409 ms, XGBoost/S23 ~0,776 ms. Dez NSGA-II somam 9.814 s (~2,73 h). Removida promessa de speedup global não medido. Tempos de inferência de todos os modelos permanecem no suplemento. Curva de aprendizado e comparação controlada de solver/treino/otimizadores não foram apresentadas como executadas.

### 21. Literatura e Forrester — corrigido e ampliado

Forrester passou a sustentar metodologia geral, não um estudo específico de blindagem. Acrescentadas fontes primárias verificadas de ML de atenuação (Nahool 2024), otimização de blindagem (Larsen 2024), ParEGO, K-RVEA, CQR, conformal sob covariate shift e feedback shift no projeto (Fannjiang 2022). A tabela de trabalhos relacionados não reivindica revisão sistemática nem prioridade universal. A revisão de surrogates balísticos específicos ainda pode ser ampliada, sem inventar referências.

### 22. Biblioteca obrigatória e B4C — corrigido

Não é alegado benefício neutrônico em MODE P. O uso dos sete materiais é reconhecido como restrição da arquitetura, não resultado de uma seleção ótima irrestrita. As conclusões limitam a validade à biblioteca estudada.

### 23. Intervalos negativos — corrigido

Gráfico regenerado com limite inferior max(0, previsão−q), sem alterar q nem o limite superior. O método matemático e a justificativa para suporte não negativo são declarados.

### 24. Colinearidade e configurações — análise concluída

Spearman na união: dose–massa −0,9953, dose–custo −0,6763, massa–custo +0,6229. Discutida a quase unidimensionalidade dose–massa. Tabela dos 30 projetos contém ordem, espessuras, dose/penetração MCNP/LS-DYNA e objetivos determinísticos. Frequências de posição de materiais foram incluídas, sem interpretação causal de uma restrição não ativa.

### 25. Edição, figuras e disponibilidade — corrigido com pendência administrativa

Removidas notas de produção, “Thesis methodological framework”, duplicação de disponibilidade, figuras desatualizadas e arte pictórica antiga. Esquemas novos são vetoriais e gerados por código a partir da geometria declarada. Fluxograma abre a metodologia. Abstract tem 209 palavras na versão-base desta revisão. CRediT preservado. Suplemento acompanha o PDF.

DOI, depósito público, eventual autorização institucional e identificação de versão pública da tese dependem dos autores. Não foram inventados. A página oficial do Guide for Authors retornou 403 nesta consulta, portanto não foi alegada conferência integral das normas atuais.

## Comentários gerais A–I e perguntas H

- A/B/D/F/G: introdução, contribuições, discussão e conclusões foram reescritas para retirar promessas não demonstradas. A diferença entre referência numérica e validação experimental é explícita. As fases são definidas antes de seu uso. Incerteza constante não foi apresentada como localmente adaptativa.
- B.4/E: tabela crítica de literatura ampliada, sem afirmar busca sistemática. Recomendações de novos métodos são discutidas, não convertidas em resultados inexistentes.
- H.1–2: contagens por limite e máximos estão nos CSVs. O alvo geométrico não é artificialmente saturado em 60 mm. A evidência atual não cobre um regime de perfuração para validar essa extrapolação.
- H.3–5: geometria, massa, fonte, FSP, constitutivos, malha, método de amostragem e concentração nos limites foram acrescentados.
- H.6: dez checagens e sensibilidade de exclusão concluídas. Histórico DE/DF versus campanhas anteriores distinguido.
- H.7–9: dominância em HF, viés, baselines e protocolo de seleção documentados. CV aninhada permanece não executada.
- H.10: tempos medidos apresentados. Custo de obter uma frente comparável sem surrogate e curva de aprendizado não podem ser deduzidos com rigor dos registros existentes.
- I: revisão compilada e auditada deve passar por leitura/autorização dos autores. Não há certificação automática de prontidão científica para submissão.

## O que ainda exige decisão ou trabalho novo

1. **Experimentos apenas computacionais nos surrogates:** comparação de orçamento equivalente, outro MOEA, ablação pareada de margem, CV aninhada/repetida, aprendizado por tamanho de amostra e incerteza local. Podem ser feitos sem novas campanhas físicas, mas consomem tempo e precisam de protocolo congelado antes da execução.
2. **Novas simulações/validação física:** convergência de malha, extensão temporal amostral e casos próximos de um limite relevante. A revisão não autoriza supor que foram feitos. Ensaios experimentais exigem evidência externa.
3. **Ações dos autores:** depósito/DOI, aprovação de divulgação e revisão final de pressupostos físicos, custos e sobreposição com tese/SPHIN.

O artigo pode ser defendido como estudo de caso com essas limitações. Para afirmar que todas as recomendações experimentais do coautor foram satisfeitas, os itens acima precisam ser efetivamente executados. Uma reescrita não substitui essa evidência.
# Atualização em 14/10/2026 — comparações de otimização concluídas

Esta atualização prevalece sobre as pendências de otimização descritas na resposta original abaixo. Foram concluídas 100 execuções, com 60.300 avaliações cada, preservando os modelos e o conjunto original de 30 casos da Fase 5.

- Os itens relativos a comparação de MOEAs/busca aleatória em orçamento igual e ablação da margem agora têm evidência executada. O artigo e a seção S8 do suplemento apresentam as condições, todas as réplicas, métricas, contrastes pareados e limites da interpretação.
- Em 60 mm, HV médio: NSGA-II 0,752552, NSGA-III 0,761918, busca aleatória com arquivo limitado 0,720850. Os dois contrastes de método envolvendo NSGA-II têm p ajustado por Holm 0,01171875.
- A margem não muda as frentes em 60 mm. Em 15/20/30 mm, muda os projetos e rejeita alguns membros das frentes pontuais, mas nenhum contraste de HV com/sem margem é significativo após Holm. Não foi atribuída uma vantagem de segurança a esses resultados.
- A nova definição de HV preserva coordenadas melhores que o ideal do DoE, sem clipping. O artigo distingue os números deste benchmark dos números históricos.
- O menor generalized spread da busca aleatória é reportado junto com seus piores HV/IGD. Não se afirma superioridade uniforme de um método em todas as métricas.
- Os 30 casos reavaliados pelos solvers continuam pertencendo à campanha original NSGA-II. Nenhum caso NSGA-III ou candidato de limite reduzido foi apresentado como reavaliado em alta fidelidade.
- Permanecem fora desta atualização: novos ajustes de modelos, CV aninhada/repetição de splits, convergência de malha e novas extensões temporais. O depósito público aguarda o fechamento da metodologia.

---

