# Avaliação das seis execuções

Verificação em 15/10/2026. Nenhuma nova simulação iniciada. Originais, datasets e artigo não alterados.

## Evidência conferida

Seis execuções com retorno zero e normal termination nos arquivos messag, histórias nodais até aproximadamente 0,30 ms e ausência de MESH_RUNNING.lock. SHA256 dos seis decks, seis originais e seis nodout conferidos contra os manifestos e postprocess.json. Penetração e razão energética reextraídas dos nodout/glstat pelo parser existente com os IDs específicos de cada malha, sem modificar arquivos de solver. Essa reextração verifica consistência, não constitui um segundo parser independente.

## Penetração geométrica registrada

Valores em mm. M1/M2: menor coordenada X entre todos os nós monitorados do FSP. A coluna de cantos usa apenas os oito cantos correspondentes ao elemento original.

| Caso | Original | M1 todos | M1 cantos | M2 todos | M2 cantos |
|---|---:|---:|---:|---:|---:|
| 9012 | 0,879 | 9,220 | 0,638 | 58,068 | 30,302 |
| 9015 | 1,732 | 12,590 | 12,590 | 10,235 | 0,450 |
| 9021 | 18,485 | 30,443 | 30,443 | 19,804 | 19,523 |

Os valores não demonstram independência de malha. A proximidade isolada entre original e M2 no 9021 não elimina a excursão em M1. A diferença grande entre todos os nós e os cantos exige examinar deformação, contato e eventual condição de elementos/nós removidos antes de interpretar o extremo nodal como trajetória física de fragmento ativo.

## Diagnóstico energético

| Caso | E/E0 final M1 | Maior desvio M1 | E/E0 final M2 | Maior desvio M2 |
|---|---:|---:|---:|---:|
| 9012 | 1,00087 | 0,106% | 0,656466 | 34,3543% |
| 9015 | 1,00022 | 0,058% | 1,00010 | 0,025% |
| 9021 | 1,03949 | 4,841% | 0,946429 | 5,3574% |

Valores da linha total energy / initial energy de glstat. Não são da linha energy ratio w/o eroded energy. Não somar novamente energia erodida ao indicador sem reconstruir sua definição. O 9012/M2 apresenta desvio material e o 9021/M2 ultrapassa o limiar histórico de 5%. Não relaxar o limiar para aprovar o estudo. Causa dos desvios ainda não estabelecida.

No 9012/M2, o mínimo X ocorre no último estado nodal (0,29996069 ms), logo o pico também não está demonstrado como estabilizado. A conclusão normal do solver não valida essa resposta.

## Limites da interpretação e próximo passo

O estudo refinou simultaneamente blindagem e FSP deformável. Portanto ainda não separa o efeito de cada discretização. A definição baseada em extremos nodais também precisa de auditoria em presença de deformação intensa/erosão. Não atribuir as diferenças somente ao material ou à malha da blindagem, nem concluir que a malha mais fina é automaticamente a verdade física.

Próximo passo recomendado, sem solver: inspecionar d3plot, matsum, sleout e eventos de eliminação nos casos críticos, rastrear o nó responsável pelo pico e confrontar o estado ativo dos elementos do FSP. Depois decidir se é necessário estudo desacoplado de refinamento de blindagem e projétil. Não iniciar outra rodada automaticamente.

Não substituir os labels originais nem recalibrar materiais. Não calcular GCI ou ordem assintótica com esses resultados. A adequação física da discretização balística permanece uma limitação material, que não é resolvida por boas métricas do surrogate nem por redação. Nenhuma conclusão de segurança real decorre destas execuções.

## Fontes locais

- mesh_preparation_manifest.json, execution_progress.jsonl e auditoria_historicos.json nesta pasta.
- malhas/sim_9012_M1 e M2, sim_9015_M1 e M2, sim_9021_M1 e M2: postprocess.json, nodout, glstat e messag.
