# Diagnóstico dos extremos nodais e da sensibilidade à malha

Data: 15/10/2026. Escopo: três casos originais e seis refinados. Sem execução de solvers, alteração de inputs, datasets ou manuscrito.

## Achado principal: pico extraído depois da remoção de um nó

No sim_9012_M2, a penetração de 58,068 mm foi calculada a partir do nó 121316, com X=351,932 mm em t=0,29996069 ms. O messag e o d3hsp registram a remoção desse nó em t=0,035713 ms. Seus elementos adjacentes 113351 e 113348 apresentam volume negativo em 0,028598 e 0,035713 ms, respectivamente. Portanto esse extremo posterior não é evidência de penetração de material ativo.

O parser usado no estudo considera todos os nós solicitados no nodout em todos os instantes, sem máscara de atividade. Esta limitação foi material no caso 9012/M2. Não se deve interpretar a tabela preliminar como resultado físico validado. Os valores brutos foram preservados para rastreabilidade.

Reextração diagnóstica: excluir o nó a partir de sua remoção registrada e excluir nós sem nenhum elemento adjacente ainda ativo segundo os eventos registrados. Incluir tanto os eventos `solid element ... failed` como `negative volume in solid element # ...`. A segunda forma não aparece na primeira expressão e deve ser contabilizada separadamente.

## Comparação após filtragem diagnóstica

| Caso | Original (mm) | M1 (mm) | M2, máscara baseada em logs (mm) |
|---|---:|---:|---:|
| 9012 | 0,879 | 9,220 | 9,772 |
| 9015 | 1,732 | 12,590 | 10,235 |
| 9021 | 18,485 | 30,443 | 19,804 |

Esta máscara depende da completude dos logs e de seus tempos arredondados. Não houve leitura binária dos estados de atividade do d3plot nesta auditoria. Portanto a coluna M2 é diagnóstico, não nova verdade de referência nem autorização para substituir labels. As outras oito estimativas não mudaram pela filtragem descrita.

## Falhas por volume negativo do FSP

- Originais: nenhum evento identificado nos três casos, cada um com um sólido de FSP.
- M1: nenhum evento identificado nos três casos, cada um com oito sólidos de FSP.
- 9012/M2: dez de 27 sólidos do FSP com volume negativo.
- 9015/M2: nenhum evento identificado nos 27 sólidos.
- 9021/M2: dois de 27 sólidos com volume negativo.

Ausência de evento no log não prova ausência de qualquer outro problema numérico. Falha por volume negativo não deve ser descrita automaticamente como fragmentação física prevista por critério constitutivo.

## Energia de contato

Energia final de interface deslizante em glstat, em J:

| Caso | Original | M1 | M2 |
|---|---:|---:|---:|
| 9012 | 0,52503 | -70,668 | -27,3739 |
| 9015 | 12,6621 | -14,6842 | 2,88499 |
| 9021 | 2,20529 | -352,308 | -81,8219 |

Como escala, a energia cinética inicial nominal do FSP é 0,5 x 0,0011 x 650² = 232,375 J. A magnitude da energia de interface em 9021/M1 excede essa escala. Isso exige investigação do contato, inclusive possíveis compensações entre termos, mesmo com E/E0 global próximo de um. Não prova, isoladamente, a causa única da divergência.

No 9012/M2, o sleout final localiza as contribuições não nulas nas interfaces 206 e 207: aproximadamente -2,918 e -24,45 J, somando cerca de -27,37 J. A razão energética final de 0,656466 continua inadequada após corrigir a interpretação do extremo nodal. A filtragem não altera o balanço energético.

## Conclusão científica

Há uma limitação identificada no pós-processamento de nós removidos, junto com evidências de sensibilidade numérica real. Corrigir apenas a primeira não estabelece convergência. Não escolher só os oito cantos porque isso produziria valores mais próximos dos originais.

Blindagem e projétil foram refinados simultaneamente. As diferenças não identificam, sozinhas, qual discretização é dominante. O estudo também altera a resolução geométrica das superfícies de contato. Não atribuir tudo à malha da blindagem, à lei de material, ou ao algoritmo de ML.

Próximo passo proposto: auditar visualmente os estados ao redor de 0,0286–0,0357 ms no 9012/M2, conferir atividade/forma dos sólidos do FSP e o contato. Se forem autorizadas novas execuções, separar refinamento da blindagem com FSP original e refinamento do FSP com blindagem original, começando no caso 9012. Isso isola fatores, não constitui por si só prova de convergência. Não mudar contato, erosão ou materiais simultaneamente para obter aprovação.

Antes de adotar pós-processamento com atividade como método oficial, verificar também os 30 casos de validação e a campanha de treinamento para incidência desse mecanismo. Nenhuma alegação de que todos os 1996 labels estejam afetados é sustentada por esta auditoria.

## Reprodutibilidade

Script diagnóstico: diagnosticar_extremos.py. Saída: diagnostico_extremos_20261015.json, com caminhos e SHA256 de messag, nodout, glstat e matsum dos nove casos. Script não invoca LS-DYNA e não sobrescreve saídas dos solvers. O resultado permanece separado do pipeline oficial.
