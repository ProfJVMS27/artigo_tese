# Estudo direcionado de convergência numérica

Data: 14/10/2026. Estado: auditoria inicial concluída, novas execuções não iniciadas.

## Escopo e rastreabilidade

Base: Fase 5/VALIDACAO_ALTA_FIDELIDADE/fase4_20261003_043316_MCNP63_DEDF_AP_30/study_30.
Auditoria reproduzível: `auditar_historicos.py`. Evidência com caminhos e SHA256: `auditoria_historicos.json`.
O script lê os 30 outputs MCNP, os 30 nodout originais e os dois nodout de duração estendida. Não executa solvers nem altera entradas ou datasets.
O estudo é direcionado aos candidatos Pareto. Não demonstra convergência de todas as 1.996 configurações de treinamento.

## Resultados existentes

- MCNP: 30/30 outputs apresentam aprovação nos dez testes do bin TFC do tally 4, com 100 milhões de histórias.
- Erro relativo final: 0,0114 a 0,0162. É uma estimativa de um desvio padrão, não erro absoluto conhecido.
- Entre os pontos TFC com N >= 50 milhões, o maior desvio de uma média acumulada em relação à média final é 1,7409%, considerando os 30 casos.
- As médias acumuladas são correlacionadas. Não aplicar teste de duas amostras independentes entre pontos do mesmo TFC.
- sim_9030 preserva o output como output_9030.txt, não .o. A auditoria utiliza o arquivo real, sem renomeá-lo. O CSV consolidado ainda registra o caminho antigo .o.
- LS-DYNA 9016: pico 3,961 mm em 0,027880254 ms. Na execução até 0,59982799 ms, o máximo após 0,30 ms não aumentou na precisão registrada. Velocidade final -2,238533125 m/s, classificação ainda indeterminada pelo parser.
- LS-DYNA 9026: pico 4,042 mm em 0,028972992 ms. Na execução até 0,59994156 ms, o máximo após 0,30 ms não aumentou. Velocidade final -0,489090875 m/s, classificado parado.
- Isso sustenta insensibilidade do máximo observado à extensão de tempo nesses dois casos, não convergência temporal universal ou parada comprovada do 9016.

## Seleção anterior às novas execuções

Malha LS-DYNA: 9012, 9015 e 9021. Regra determinística: mínimo, mediana inferior (posição 15 de 30) e máximo da penetração original, desempate por ID. Valores: 0,879 mm, 1,732 mm e 18,485 mm. Não selecionar pelo erro do surrogate. Documentar a composição dos três, sem alegar que cobrem todas as classes de materiais.

MCNP: 9001, 9014 e 9030. Mínimo, mediana inferior e máximo da dose simulada. O 9001 também possui o maior erro relativo, portanto não é necessário adicionar um quarto caso para esse critério.

## Etapa A: MCNP, sem nova campanha neste momento

Os históricos existentes justificam começar por gráficos de média, R, VOV e FOM versus N, para os três casos, acompanhados da auditoria de todos os 30. Não declarar independência física a partir dos dez testes.

Se for necessária confirmação adicional da precisão, executar somente esses três casos com 400 milhões de histórias e sementes novas documentadas. O fator quatro permitiria esperar aproximadamente metade do erro relativo se o regime de variância finita e amostragem estável for mantido. Não é garantia.

Comparar estimativas independentes com suas incertezas combinadas. Registrar a diferença em mSv e em desvios padrão combinados. Um intervalo aproximado de 95% exige hipóteses estatísticas adequadas e não é garantia de exatidão do modelo. Preservar DE/DF, FM, materiais, geometrias, bibliotecas e normalização por 300 s.

Não usar 2% ou 5% como selo universal. A tolerância deve depender da diferença de dose que se pretende distinguir entre projetos. Se diferenças Pareto forem menores que a incerteza combinada, informar resolução insuficiente para ordenar esses projetos.

## Etapa B: malha LS-DYNA

Planejar três níveis por caso: original e dois refinamentos conformais coerentes na região de impacto e no projétil. Se viáveis, seis novas execuções mais os três resultados originais.

Antes de gerar os decks, há duas verificações obrigatórias:

1. A entrada original inspecionada possui 19.853 nós e 17.473 sólidos. Confirmar o limite efetivo da licença Student instalada e contar os elementos/nós de cada proposta. Não assumir que dois refinamentos cabem na licença.
2. O parser atual fixa os nós FSP em 19846–19853. A malha refinada precisa de seleção dos nós por parte e de uma definição equivalente da penetração e velocidade. Não reutilizar os IDs originais nem trocar silenciosamente média nodal por média de massa. Auditar a equivalência da medida para o projétil rígido/refinado antes de comparar.

Usar os decks efetivamente executados como referência, não confiar apenas no nome `modelo_mecanico_congelado.py`. O gerador presente contém comentários de versão V3, que precisam ser reconciliados com os decks originais.

Manter constantes geometria física, composição, espessuras, massa e velocidade do FSP, leis constitutivas, erosão, formulação sólida, contato, restrições e controle de mass scaling. Gerar exclusivamente cópias de estudo. Não calibrar parâmetros para reduzir diferenças entre malhas.

Registrar tamanho local dos elementos nas três direções, contagens, qualidade e passo de tempo. Preferir razões de refinamento comparáveis, mas só fixá-las após a verificação de licença e topologia. Refinar apenas uma direção não demonstra independência geral de malha.

Comparar penetração máxima, instante do pico, trajetória, velocidade, energias, massa erodida e regime de contenção. Usar o mesmo horizonte de tempo nos três níveis. Verificar resolução temporal da saída em torno do pico. Duração final não é equivalente a estudo de passo de tempo.

Apresentar diferenças absolutas em mm e relativas entre níveis. Como referência de relevância, comparar com RMSE balístico de 2,306 mm e margem conformal de 3,933 mm, sem tratar esses números como tolerância física. Mudanças que alterem contenção, regime de erosão ou interpretação do artigo impedem uma conclusão de independência.

Não calcular ordem de convergência ou GCI automaticamente para resultados não monotônicos/discontínuos por erosão. Se não houver aproximação estável entre níveis, relatar sensibilidade à malha e a limitação correspondente.

## Etapa C: duração

Não repetir automaticamente os 30 casos. Preservar a comparação existente de 9016 e 9026. Se for necessário classificar o estado final do 9016, preparar extensão específica, inicialmente até 1,20 ms, como sensibilidade de duração, sem prometer término da interação. A penetração máxima e a classificação final devem permanecer como conclusões distintas.

## Entrega e limites

Tabelas e gráficos com todos os níveis executados, sem excluir resultados desfavoráveis. Falhas também devem aparecer. Novos estudos não substituem automaticamente os labels de treinamento, não exigem rerodar as fases 3/4 e não autorizam afirmar que estas convergiram. Decidir eventual impacto no artigo somente após os resultados.

Não foi iniciado nenhum solver. A próxima preparação operacional é reconciliar o gerador com os decks originais e definir refinamentos que caibam na licença. Ainda não há comando de simulação seguro para entregar.

## Referência de interpretação MCNP

Manual MCNP 6.3.1, LANL LA-UR-24-24602 Rev.1, discussão de TFC, FOM e dez testes estatísticos:
https://mcnp.lanl.gov/pdf_files/TechReport_2024_LANL_LA-UR-24-24602Rev.1_KuleszaAdamsEtAl.pdf
O solver usado nos outputs auditados é 6.3.0. A referência é usada para interpretação dos diagnósticos, não como identificação da versão executada.
