# Isolamento de discretização — caso 9012

## Objetivo e escopo

Este estudo diagnóstico separa a refinação da blindagem da refinação do FSP no caso 9012. Ele não altera o modelo mecânico congelado, a campanha original ou os respectivos rótulos do dataset. Os resultados são evidência de sensibilidade numérica localizada, não uma certificação de convergência de malha.

## Execuções concluídas

| Configuração | Terminação | Penetração bruta (mm) | Penetração ativa por log (mm) | E/E0 final | Máx. |E/E0-1| | Diagnóstico FSP |
|---|---:|---:|---:|---:|---:|---|
| Original / original | normal | 0.879 | 0.879 | 1.00537 | 0.00537 | sem erosão registrada |
| Blindagem M1 / FSP original | normal | 3.657 | 3.657 | 1.07232 | 0.07375 | FSP íntegro; 52 falhas em sólidos da blindagem |
| Blindagem original / FSP M1 | normal | 128.096 | 0.879 antes da remoção | 0.523659 | 0.476357 | 4/8 elementos do FSP falharam por volume negativo em 0.022394–0.022476 ms |
| Blindagem M1 / FSP M1* | normal | 9.220 | 9.220 | 1.00087 | 0.00087 | FSP íntegro |

\*Resultado proveniente do estudo M1 anterior, incluído apenas para interpretar a interação entre as duas discretizações.

## Interpretação

1. Refinar apenas a blindagem elevou a penetração de 0.879 para 3.657 mm. Embora o FSP tenha permanecido íntegro, o desvio energético máximo foi 7.375%, acima da tolerância operacional de 5% adotada para a campanha. Portanto, esse resultado não certifica convergência.
2. Refinar apenas o FSP produziu erosão precoce do próprio projétil e desequilíbrio energético de 47.636%. O pico bruto de 128.096 mm pertence a nó removido e não tem significado físico. O valor de 0.879 mm filtrado por atividade é apenas a maior posição anterior à falha, não um resultado convergido.
3. O resultado M1/M1 é energeticamente estável, mas é aproximadamente dez vezes maior que o da malha original. Com o nível M2 já diagnosticado como instável por erosão e desequilíbrio energético, não existe sequência de malhas válida e monotônica que permita estimar erro de discretização ou aplicar extrapolação.

## Conclusão auditável

As quatro combinações demonstram que a resposta do caso 9012 é sensível à interação entre discretização da blindagem, discretização do FSP, contato e erosão. A configuração congelada não fornece, para este caso, uma demonstração formal de independência de malha. Não se deve substituir nenhum rótulo da campanha original pelos resultados refinados nem alegar convergência de malha no artigo.

## Implicação para o artigo

É defensável relatar que: (i) as soluções de validação de alta fidelidade foram executadas no modelo numérico pré-especificado; (ii) a análise de sensibilidade temporal disponível foi caso-específica; e (iii) uma futura verificação de independência de malha requer uma formulação de refinamento do FSP que preserve estabilidade de contato, erosão e balanço energético antes de qualquer reexecução de campanha. O texto não deve sugerir validação experimental nem convergência numérica geral.
