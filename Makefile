CC = gcc
CFLAGS = -O3 -march=native -Wall -Wextra
LDLIBS = -lm
PYTHON = /opt/anaconda3/bin/python

pivot: pivot.c
	$(CC) $(CFLAGS) $< $(LDLIBS) -o $@

reports/work/baseline: tests/baseline.c
	mkdir -p reports/work
	$(CC) $(CFLAGS) $< $(LDLIBS) -o $@

reports/work/test_score: tests/test_score.c pivot.c
	mkdir -p reports/work
	$(CC) $(CFLAGS) $< $(LDLIBS) -o $@

reports/work/test_sort: tests/test_sort.c pivot.c
	mkdir -p reports/work
	$(CC) $(CFLAGS) $< $(LDLIBS) -o $@

reports/work/test_ranking: tests/test_ranking.c pivot.c
	mkdir -p reports/work
	$(CC) $(CFLAGS) $< $(LDLIBS) -o $@

reports/work/test_bounds: tests/test_bounds.c pivot.c
	mkdir -p reports/work
	$(CC) $(CFLAGS) $< $(LDLIBS) -o $@

check: pivot reports/work/baseline reports/work/test_score reports/work/test_sort reports/work/test_ranking reports/work/test_bounds
	taskset -c $${BENCH_CPU:-100} $(PYTHON) tests/check.py reports/work/baseline ./pivot
	taskset -c $${BENCH_CPU:-100} reports/work/test_score
	taskset -c $${BENCH_CPU:-100} reports/work/test_sort
	taskset -c $${BENCH_CPU:-100} reports/work/test_ranking
	taskset -c $${BENCH_CPU:-100} reports/work/test_bounds
	taskset -c $${BENCH_CPU:-100} $(PYTHON) tests/check_invalid.py ./pivot

bench: pivot
	bash reports/benchmark.sh ./pivot current 3

.PHONY: check bench
