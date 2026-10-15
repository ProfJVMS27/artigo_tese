# Execução do estudo de malha

As seis entradas estão em `malhas`. Originais não modificados. Manifesto: `mesh_preparation_manifest.json`.

| Nível | Nós | Sólidos | Subdivisões radiais por camada | Passo angular central | Passo vertical central | FSP |
|---|---:|---:|---:|---:|---:|---|
| Original | 19853 | 17473 | 2 | 0,5 grau | 5 mm | 1 sólido |
| M1 | 56193 | 51416 | 3 | 0,25 grau | 2,5 mm | 8 sólidos |
| M2 | 121371 | 113371 | 4 | 0,166667 grau | 1,666667 mm | 27 sólidos |

Casos: 9012, 9015 e 9021. Ordem: M1 e M2 de cada caso. Tempo físico 0,30 ms, igual ao original. Sem mass scaling. Elementos/nós únicos, conectividade conformal, Jacobianos positivos nos oito pontos de Gauss e volumes verificados. Os cartões não relacionados à malha são idênticos aos originais. A reconstrução usa a topologia dos decks executados, não o gerador congelado.

O limite publicado é de 128 mil nós/elementos: https://ansys.synopsys.com/academic/students/ansys-ls-dyna-student
A aceitação efetiva pelo executável instalado só será confirmada na execução. Não foi iniciado um solver para esse fim. M2 está próximo do limite, embora abaixo dele.

O refinamento angular altera ligeiramente a aproximação poligonal da superfície curva. A mudança de volume por parte é registrada no manifesto e é menor que 0,1%. Não é uma alteração intencional de densidade, massa física ou espessura. O volume geométrico do FSP permanece equivalente dentro da precisão do arquivo.

O FSP é deformável (*MAT_JOHNSON_COOK), não rígido. Seu refinamento pode revelar deformações que o sólido único não representa. A penetração obtida do menor X de todos os nós e aquela dos oito cantos originais serão apresentadas separadamente. A média nodal de velocidade NÃO é velocidade do centro de massa. Não substituir automaticamente os targets do dataset com resultados novos.

## Comando PowerShell

```powershell
& "C:\Users\Admin\anaconda3\envs\spyder_dl8\python.exe" -u "C:\Users\Admin\Desktop\Scripts Tese 1\Fase 5\CONVERGENCIA_20261014\executar_malhas.py" --run
```

Console em tempo real, um solver por vez, três CPUs por padrão. São seis execuções novas, não uma repetição da campanha. Não há estimativa confiável de duração antes do primeiro caso refinado. A malha fina pode exigir muito mais tempo e espaço em disco que a original.

Sem `--run`, o programa só verifica hashes e mostra o plano. Após cada caso concluído, salva `completed.json`, `postprocess.json` e uma linha em `execution_progress.jsonl`. O checkpoint é por caso concluído, não retomada intrassimulação.

Ao repetir o comando, casos com conclusão verificada são ignorados. Saídas parciais são preservadas e bloqueiam sobrescrita. Se interromper, não apague o lock nem relance enquanto houver solver ativo. A recuperação de uma execução parcial exige revisão, não remoção indiscriminada de arquivos.

Não chamar o summarizer antigo diretamente nas malhas refinadas: ele usa IDs fixos. O executor deste estudo realiza a extração com o manifesto de cada malha.

Conclusão normal, ausência de erro e balanço energético adequado não provam independência de malha. Avaliar os três níveis, inclusive resultados desfavoráveis, após terminar. É um estudo local de sensibilidade/convergência, não refinamento uniforme de todo o domínio ou prova para os 1996 casos.
