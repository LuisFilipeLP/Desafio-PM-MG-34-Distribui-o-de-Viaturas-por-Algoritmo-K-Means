# Desafio-PM-MG-34-Distribui-o-de-Viaturas-por-Algoritmo-K-Means
Solução do desafio proposto pela PM-MG 34° na edição de 2026 do Hackatohm, UFCI-MG.

Equipe Faraday

Membros:
Luís Filipe Lima Porto    20243004239
João Vitor Ferreira      20253015521

Nossa equipe desenvolveu um protótipo de sistema de apoio à decisão focado na otimização e alocação dinâmica de policiamento preventivo orientado por dados. A solução proposta transforma o histórico de ocorrências policiais em estratégias operacionais eficientes, equilibrando a dissuasão criminal em áreas críticas com a garantia de cobertura territorial abrangente.
Metodologia e Arquitetura Técnica:
1. Modelagem Espacial e Ponderação: O território é segmentado em uma malha geométrica regular utilizando hexágonos de 200 metros de diâmetro. Cada ocorrência histórica é espacializada e processada individualmente, aplicando-se o Índice de Gravidade Acumulada baseado nos pesos oficiais do edital (de A a F).
2. Otimização por Algoritmo de Agrupamento (Weighted K-Means): Para realizar a distribuição das viaturas, implementamos o algoritmo de K-Means Ponderado (Machine Learning/Clusterização). O número de agrupamentos (K) é parametrizado dinamicamente com base na quantidade de viaturas ativas no momento. O algoritmo processa as coordenadas geográficas dos hexágonos utilizando o Índice de Gravidade como peso (sample weight), calculando centroides matemáticos ideais que representam os focos de calor criminal para aquela faixa horária.
3. Geração Automatizada e Restrita do Cartão-Programa: O sistema analisa as 24 faixas horárias do dia e as divide automaticamente nas Classes 1, 2 e 3 de incidência. As viaturas são alocadas nos centroides de seus respectivos clusters. Para cumprir estritamente a regra de alternância da Classe 1, inserimos uma restrição global no código (através de um filtro de exclusão temporal) que impede que uma mesma quadrícula receba Pontos-Bases em horas consecutivas, direcionando o recurso para a segunda região mais crítica daquele horário caso haja sobreposição.
4. Reação Dinâmica em Tempo Real: O grande diferencial da solução é o seu motor de recálculo instantâneo. Caso uma ou mais viaturas fiquem indisponíveis devido a empenhos operacionais, o sistema atualiza o valor de K, roda o algoritmo novamente e redesenha as áreas de cobertura das viaturas remanescentes em tempo real, mitigando "vazios urbanos".
Interface e Visualização:
O sistema conta com um front-end interativo desenvolvido em HTML, CSS e JavaScript (utilizando a biblioteca Leaflet) conectado a um motor em Python. No mapa, os avaliadores conseguem visualizar as manchas de calor dos hexágonos, a posição atual das viaturas, seus raios de atuação (1km, 2km ou 3km dependendo da classe horária) e os respectivos Cartões-Programa gerados de forma automatizada.
