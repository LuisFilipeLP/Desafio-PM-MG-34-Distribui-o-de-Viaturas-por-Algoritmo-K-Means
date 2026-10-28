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

# Cartão-Programa dinâmico – 39ª CIA PM (Hackathon PMMG)

## Como executar
```
pip install numpy scipy
python3 patrulha.py --k 5 --porta 8000     # abre o servidor; acesse http://localhost:8000
```
O navegador precisa de internet apenas para Leaflet e tiles OpenStreetMap/CARTO. Tudo é código aberto (Python, NumPy, SciPy, Leaflet, OSM).

## Arquitetura
- `algoritmo_kmeans.py` – algoritmo original (K-Means balanceado, Húngaro, trocas), reutilizado sem alterações.
- `patrulha.py` – backend: roda o K-Means por hora, escolhe as quadrículas, gera os Cartões-Programa e serve a API (`/api/base`, `/api/plano?indisp=V2,V4`).
- `index.html` – interface (Leaflet): mapa, animação, cartões, indisponibilização.
- Dados: `hexagonos_200m.geojson` (reconstruído do mapa_hexagonos_39cia.html: peso, nível e nº de ocorrências por hexágono/hora), `faixas_horarias_classes.csv`, `39_cia.geojson`, `ocorrencias.txt`.

## Metodologia
1. Peso por ocorrência (A=4, B=5, C=1, D=1, E=2, F=3); soma por hexágono/hora; níveis 1–5; prioritários = níveis 4–5.
2. K-Means ponderado e balanceado (K = viaturas disponíveis) sobre os prioritários de cada hora; demais hexágonos vão para o centroide mais próximo (cobertura 100%). Rótulos casados com a hora anterior (menor deslocamento).
3. Classes das faixas (CSV): Classe 1 = ponto-base 20 min + patrulha 40 min (raio 1 km); Classe 2 = 4 faixas iguais à Classe 1 e 4 com patrulhamento até 2 km; Classe 3 = patrulhamento até 3 km.
4. Quadrícula em destaque de cada viatura/hora = hexágono de maior peso da sua área que não foi selecionado na hora anterior nem na seguinte (circular, 23h→0h) por nenhuma viatura. Verificado por `assert` a cada plano.
5. O início do ponto-base varia entre :05 e :35 dentro da faixa (por viatura/hora); a patrulha é dividida antes/depois, somando 40 min.
6. Indisponibilidade: a viatura sai do conjunto, o K-Means roda de novo com K-1 e todos os Cartões são regerados; o painel mostra hexágonos reatribuídos (comparação com o cenário base), área e peso por viatura.

## Limitações
- Rotas são linhas retas entre quadrículas (sem malha viária); localização textual é por hexágono + coordenadas (não há base de logradouros offline).
- Replanejamento recalcula as 24 h (a animação suaviza a transição em 4 min); não há ainda "volta automática" após o empenho.
- Tipos de viatura (atendimento/preventiva) são premissa em `TIPOS` no `patrulha.py`; hoje só identificam o cartão.
- Recálculo leva alguns segundos (cache por combinação de viaturas).

## Evoluções
Roteamento por vias (OSRM/OSM), tempo real de empenho, tipos de viatura com regras distintas, tempo de atendimento estimado, métricas de tempo-resposta.
