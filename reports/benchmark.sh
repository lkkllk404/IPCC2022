#!/usr/bin/env bash
# Run from the repository root; one pinned CPU, serial executions, no GPU.
set -euo pipefail
binary=${1:?binary}
label=${2:?label}
repeats=${3:-3}
mkdir -p reports/work
printf 'date=%s\ncpu=%s\n' "$(date --iso-8601=seconds)" "${BENCH_CPU:-100}"
uptime
free -m | head -2
for ((i=1;i<=repeats;i++)); do
    /usr/bin/time -f "elapsed=%e rss_kib=%M" -o "reports/work/$label-$i.time" \
        timeout 180 taskset -c "${BENCH_CPU:-100}" "$binary" uniformvector-2dim-5h.txt \
        > "reports/work/$label-$i.log"
    cat "reports/work/$label-$i.time"
    grep 'Using time' "reports/work/$label-$i.log" | tail -1
    cmp result.txt refer-2dim-5h.txt
    echo 'official reference: byte-identical'
done
