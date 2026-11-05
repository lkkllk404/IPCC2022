CC = gcc
CFLAGS = -O3 -march=native -Wall -Wextra
LDLIBS = -lm
PYTHON = /opt/anaconda3/bin/python

pivot: pivot.c
	$(CC) $(CFLAGS) $< $(LDLIBS) -o $@

reports/work/baseline: tests/baseline.c
	mkdir -p reports/work
	$(CC) $(CFLAGS) $< $(LDLIBS) -o $@

check: pivot reports/work/baseline
	taskset -c $${BENCH_CPU:-100} $(PYTHON) tests/check.py reports/work/baseline ./pivot

bench: pivot
	bash reports/benchmark.sh ./pivot current 3

.PHONY: check bench
