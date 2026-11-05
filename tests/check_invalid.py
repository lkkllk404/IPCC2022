#!/usr/bin/env python3
import subprocess
import sys
import tempfile
from pathlib import Path

binary = str(Path(sys.argv[1]).resolve())
cases = ['', 'x y z\n', '0 2 1\n', '2 -1 1\n', '2 2 0\n', '2 2 3\n',
         '2 2 2\n0 0\n', '2 2 2\nnan 0\n0 0\n', '2 2 2\ninf 0\n0 0\n',
         '2 2 2\n1e200 0\n0 0\n']
for text in cases:
    with tempfile.TemporaryDirectory() as tmp:
        inp = Path(tmp) / 'input.txt'
        inp.write_text(text)
        result = subprocess.run([binary, str(inp)], cwd=tmp, capture_output=True, timeout=5)
        assert result.returncode > 0, (text, result.returncode, result.stderr)
        assert not (Path(tmp) / 'result.txt').exists()
with tempfile.TemporaryDirectory() as tmp:
    inp = Path(tmp) / 'input.txt'
    inp.write_text('2 2 2\n0 0\n1 1\n')
    (Path(tmp) / 'result.txt').mkdir()
    result = subprocess.run([binary, str(inp)], cwd=tmp, capture_output=True, timeout=5)
    assert result.returncode > 0
print('PASS 11 invalid-input/output-error cases (rejected without crash).')
