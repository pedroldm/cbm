# Builds every solver used by the experiment runner (experiments/run_experiments.py).
#   make            -> ENS/ENS, ILS/ils, src/CBMLKH/main_prd, src/CBMLKH/lkh_standalone, src/LKH3/LKH
# Linkern is not built here: it is an external Concorde binary passed via --linkern-path / LINKERN_PATH.

CC ?= gcc
CXX ?= g++

all: ens ils cbmlkh lkh

ens: ENS/ENS
ENS/ENS: ENS/ENS_for_CBM.c src/common/cbm_seed.h src/common/cbm_spawn.h src/common/cbm_lkh_params.h
	$(CC) -O2 -std=c11 -Wall -o $@ $< -lm

ils: ILS/ils
ILS/ils: ILS/ils.cpp src/common/cbm_seed.h src/common/cbm_spawn.h src/common/cbm_lkh_params.h
	$(CXX) -O2 -std=c++17 -Wall -Wextra -o $@ $<

cbmlkh:
	$(MAKE) -C src/CBMLKH prd lkh_standalone

lkh:
	$(MAKE) -C src/LKH3

clean:
	rm -f ENS/ENS ILS/ils
	$(MAKE) -C src/CBMLKH clean

.PHONY: all ens ils cbmlkh lkh clean
