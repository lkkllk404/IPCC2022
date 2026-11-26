#!/usr/bin/env python3
"""Compare exact ranks and non-timing stdout against serial enumeration."""
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile

binary = str(Path(sys.argv[1]).resolve())
rng = random.Random(2022)
cases = [
    ("wide-and-tail-batches", 2, [(rng.random(), rng.random()) for _ in range(150)]),
    ("rank-cutoff ties", 2, [(i % 3, i % 2) for i in range(110)]),
    ("zero-score placeholders", 2, [(0, 0)] * 72),
    ("one pair", 2, [(0, 0), (1, 1)]),
    ("generic fallback", 3, [(i % 3, i // 3) for i in range(10)]),
]

with tempfile.TemporaryDirectory() as tmp:
    data = Path(tmp) / "input.txt"
    result = Path(tmp) / "result.txt"

    def run(threads, limit=4, check=True):
        env = dict(os.environ, PIVOT_THREADS=str(threads), OMP_DYNAMIC="FALSE",
                   OMP_THREAD_LIMIT=str(limit))
        return subprocess.run([binary, data.name], cwd=tmp, env=env, capture_output=True,
                              timeout=10, check=check)

    def outputs(threads, limit=4):
        stdout = run(threads, limit).stdout.splitlines()
        return result.read_bytes(), [line for line in stdout if not line.startswith(b"Using time : ")]

    for name, k, coords in cases:
        data.write_text(f"2 {len(coords)} {k}\n" + "".join(f"{x} {y}\n" for x, y in coords))
        serial = outputs(1)
        for threads in (2, 4):
            assert outputs(threads) == serial, (name, threads)
        if name == "wide-and-tail-batches":
            # 11175 pairs fill the warm-up and a full 8192-pair batch, then a tail.
            assert outputs(64, limit=2) == serial
        print(f"PASS {name}: ranks and progress match serial byte-for-byte")

    for invalid in ("", "0", "-1", "65", "99999999999999999999999", "2x"):
        proc = run(invalid, check=False)
        assert proc.returncode != 0 and b"PIVOT_THREADS" in proc.stderr, invalid
    print("PASS PIVOT_THREADS validates 1..64; max request tested with OpenMP limited to 2 threads")

    data.write_text("1 64 2\n" + "0\n" * 63 + "1e200\n")
    result.unlink(missing_ok=True)
    proc = run(4, check=False)
    assert proc.returncode == 1 and b"Coordinate range too large" in proc.stderr
    assert not result.exists()
    print("PASS a worker's distance overflow aborts before scoring or writing a result")
