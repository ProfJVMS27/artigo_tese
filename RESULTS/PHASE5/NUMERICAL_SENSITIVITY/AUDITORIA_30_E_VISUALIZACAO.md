# Auditoria de atividade dos 30 casos e visualização diagnóstica

Data: 15/10/2026.

## Casos de validação originais

Foram inspecionados os casos 9001 a 9030 em study_30/lsdyna. Em cada caso, a conectividade do FSP foi lida do keyword (parte 99), os eventos de falha e volume negativo foram lidos do messag e os extremos geométricos foram reextraídos do nodout. Foram considerados também os registros explícitos de remoção de nós.

Resultado: 30 casos verificados, zero eventos identificados de falha de elementos do FSP e zero nós do FSP removidos. As 30 penetrações máximas permaneceram idênticas ao excluir nós sem suporte ativo segundo os logs. Os valores e hashes estão em auditoria_30_atividade_20261015.json.

Conclusão limitada: o mecanismo de extremo pós-remoção identificado no 9012/M2 não alterou os 30 resultados originais de validação. Isso não demonstra independência de malha, não valida automaticamente toda a física e não estende o achado aos 1996 casos de treinamento, que não foram auditados nesta etapa.

## Visualização do 9012/M2

Arquivo: fsp_9012_M2_diagnostico.png. Script: visualizar_fsp_nodout.py.

A figura usa coordenadas X/Z do nodout e a conectividade do keyword. Elementos com eventos registrados de falha são mostrados em cinza tracejado a partir de seus tempos registrados. Os demais são azuis. As escalas variam entre painéis e o contorno da blindagem deformada não está representado. A linha X=410 mm é apenas a referência geométrica inicial.

Os estados mostram 27 elementos ativos em aproximadamente 0,026 ms, 22 em 0,030 ms e 17 em 0,036 ms. O estado final apresenta grande extensão das conexões cinza, que não deve ser interpretada como material ativo atravessando a blindagem. A deformação dos elementos ativos também é importante.

Limitação: esta é uma reconstrução diagnóstica, não uma captura do LS-PrePost e não uma leitura do estado de atividade binário do d3plot. O LS-PrePost 4.12 foi localizado instalado, mas a inspeção visual interativa do d3plot não foi executada nesta etapa. Não atribuir a essa figura prova de auto-interseção tridimensional ou de geometria detalhada do contato.

## Próxima decisão técnica

A evolução do FSP e a sensibilidade conjunta justificam isolar fatores, se novas execuções forem autorizadas: no 9012, testar primeiro blindagem M1 com FSP original e blindagem original com FSP M1, mantendo todos os demais cartões e as métricas diagnósticas. Comparar com os dois resultados já existentes (original/original e M1/M1). Esse pequeno arranjo 2x2 separa os efeitos iniciais e sua interação. Não é um teste definitivo de convergência e não deve substituir a documentação da falha M2.

Não foram geradas novas entradas, iniciadas simulações ou alterados resultados de treinamento/validação nesta etapa. Não é necessário rerodar automaticamente as fases 3 e 4 por causa do erro de pós-remoção, pois ele não foi identificado nestes 30 casos originais. O problema de sensibilidade à discretização permanece aberto.
