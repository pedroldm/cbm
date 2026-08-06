instancePath=/home/pedroldm/MSc/cbm/instances/a1

# true = imprime apenas o custo da melhor solucao (um inteiro em stdout), sem o
# relatorio JSON. Para o target-runner do iRace ler direto.
iRace=false

threads=8
blockMovement=RANDOM

maxIterations=1000
# Segundos
maxTime=60
# Segundos (vira TIME_LIMIT no arquivo .par do LKH, que e' em segundos)
lkhMaxTime=5

constructionBias=2.5

neighborBias=1.0
minNeighborBias=0.2

minSegmentScore=10.0
minSegmentScoreLowerBound=2.0

# Segment sizes are fractions of the instance's column count, resolved at load
# time (e.g. on a 1000-column instance: 0.1 -> 100, 0.2 -> 200). PEAK and
# INTERVAL segments never exceed the resolved maxSegmentSize (widened toward
# maxSegmentSizeUpperBound by diversification); MERGE spans at most two of them.
minSegmentSize=5
maxSegmentSize=0.1
maxSegmentSizeUpperBound=0.2

segmentSizeGrowthFactor=1.25
segmentScoreDecayFactor=0.85
neighborBiasDecayFactor=0.90

adaptationInterval=20

# Criterio de parada por histograma: interrompe a trajetoria apos N iteracoes
# consecutivas sem achatar o histograma de blocos por coluna da solucao
# incumbente. 0 desativa (a parada fica so' por maxIterations / maxTime).
histogramStopInterval=0
# Como o achatamento e' medido (menor = mais achatado em todas as medidas):
#   VARIANCE - variancia dos blocos por coluna; achatar = rebaixar os picos
#   PEAK     - altura do maior pico; so' conta quando a pior coluna cai
#   ENTROPY  - entropia de Shannon normalizada (invariante de escala): mede a
#              redistribuicao dos blocos, nao a reducao do total
#   GINI     - coeficiente de Gini (invariante de escala, O(c log c))
histogramStopMeasure=VARIANCE