#!/usr/bin/env python3
"""Differential validation. Equal-score ties may differ by rounding order."""
import itertools
import math
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile

base, candidate = [str(Path(p).resolve()) for p in sys.argv[1:3]]
rng = random.Random(2022)
cases = []
for dim, n, k in [(2, 1, 1), (2, 2, 2), (2, 16, 1), (2, 18, 2), (3, 14, 3), (5, 12, 4), (2, 50, 2)]:
    cases.append((f'random-{dim}-{n}-{k}', dim, k,
                  [[rng.random() for _ in range(dim)] for _ in range(n)]))
cases += [('duplicates', 2, 2, [[0., 0.]] * 8),
          ('grid', 2, 2, [[float(i), float(j)] for i in range(4) for j in range(4)]),
          ('line', 2, 2, [[float(i), 0.] for i in range(16)])]
for name, dim, k, coords in cases:
    with tempfile.TemporaryDirectory() as tmp:
        inp = Path(tmp) / 'input.txt'
        inp.write_text(f'{dim} {len(coords)} {k}\n' + '\n'.join(' '.join(map(str, p)) for p in coords))
        results = []
        for binary in [base, candidate]:
            subprocess.run([binary, str(inp)], cwd=tmp, stdout=subprocess.DEVNULL, check=True)
            results.append([tuple(map(int, row.split())) for row in (Path(tmp) / 'result.txt').read_text().splitlines()])
        # Independently evaluate each output combination in double precision.
        distances = [[math.sqrt(sum((a-b)**2 for a,b in zip(p,q))) for p in coords] for q in coords]
        cache = {}
        def score(pivots):
            if pivots not in cache:
                cache[pivots] = 2 * math.fsum(max(abs(distances[p][i]-distances[p][j]) for p in pivots)
                                            for i in range(len(coords)) for j in range(i+1, len(coords)))
            return cache[pivots]
        assert all(len(r) == 2000 for r in results), name
        count = min(1000, math.comb(len(coords), k))
        for section in [0, 1000]:
            # With all-zero input the legacy maximum list retains zero placeholders.
            if name != 'duplicates' or section != 0:
                valid = results[1][section:section+count]
                assert all(len(p) == k and all(0 <= x < len(coords) for x in p)
                           and list(p) == sorted(set(p)) for p in valid), name
                assert len(set(valid)) == len(valid), name
            for i in range(section, section+1000):
                a, b = score(results[0][i]), score(results[1][i])
                assert math.isclose(a, b, rel_tol=2e-12, abs_tol=2e-10), (name, i, a, b)
        print(f'PASS {name}: all 2000 ranks, valid/unique combinations; byte-identical={results[0] == results[1]}')
print('All differential checks passed (tolerance only for floating-point equal-score ties).')
