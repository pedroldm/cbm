instancePath=/home/pedroldm/MSc/cbm/instances/a1

threads=8
blockMovement=RANDOM

maxIterations=1000
# Segundos
maxTime=60
# Minutos
lkhMaxTime=1

constructionBias=2.5

neighborBias=1.0
minNeighborBias=0.2

minSegmentScore=10.0
minSegmentScoreLowerBound=2.0

# Segment sizes are fractions of the instance's column count, resolved at load
# time (e.g. on a 1000-column instance: 0.1 -> 100, 0.2 -> 200).
maxSegmentSize=0.1
maxSegmentSizeUpperBound=0.2

segmentSizeGrowthFactor=1.25
segmentScoreDecayFactor=0.85
neighborBiasDecayFactor=0.90

adaptationInterval=20