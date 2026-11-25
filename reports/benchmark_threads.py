#!/usr/bin/env python3
"""Run the official thread matrix on hzcu-hpct."""
import argparse
from collections import Counter, defaultdict
import fcntl
import os
from pathlib import Path
import re
from statistics import median
import subprocess
import tempfile

GATE_MASKS = {1: "100", 2: "100-101", 4: "100-103", 8: "100-107",
              16: "100-115", 32: "96-127", 64: "32-63,96-127"}
REPORT_MASKS = {1: "100", 2: "4,100", 4: "4-5,100-101", 8: "4-7,100-103",
                16: "4-11,100-107", 32: "4-19,100-115", 64: "0-31,96-127"}
NUMA_NODE = 1
REPEATS = 5
AFFINITY = re.compile(rb"thread=(\d+) affinity=(\d+)")
TIMING = re.compile(rb"^Using time : (\d+\.\d+) ms$", re.MULTILINE)


def expand_mask(mask):
    cpus = set()
    for group in mask.split(","):
        start, _, end = group.partition("-")
        cpus.update(range(int(start), int(end or start) + 1))
    return cpus


def verify_topology(work, env, cpu_masks, report):
    topology = subprocess.check_output(["lscpu", "-p=CPU,SOCKET,CORE,NODE"], text=True, env=env)
    (work / "topology.csv").write_text(topology)
    cpu_info = {}
    for line in topology.splitlines():
        if not line.startswith("#"):
            cpu, socket, core, node = line.split(",")
            cpu_info[int(cpu)] = (socket, core, node)
    allowed = os.sched_getaffinity(0) & cpu_info.keys()
    placement = {}
    for threads, mask in cpu_masks.items():
        cpus = expand_mask(mask)
        if len(cpus) != threads or not cpus <= allowed:
            raise SystemExit(f"{threads} threads need CPUs {mask} available to this checkout")
        cores = Counter(cpu_info[cpu][:2] for cpu in cpus)
        nodes = Counter(cpu_info[cpu][2] for cpu in cpus)
        expected_nodes = ({"0": threads // 2, "1": threads // 2} if report and threads > 1
                          else {str(NUMA_NODE): threads})
        per_core = 2 if not report and threads == 64 else 1
        if nodes != expected_nodes or set(cores.values()) != {per_core}:
            raise SystemExit(f"{threads} threads: CPUs {mask} do not match the expected NUMA / SMT layout")
        placement[threads] = (",".join(f"{node}:{count}" for node, count in sorted(nodes.items())),
                              len(cores), per_core)
    return placement


def benchmark(root, work, report):
    cpu_masks = REPORT_MASKS if report else GATE_MASKS
    env = dict(os.environ, LC_ALL="C", OMP_DYNAMIC="FALSE", OMP_PROC_BIND="close",
               OMP_PLACES="threads", OMP_WAIT_POLICY="PASSIVE", OMP_DISPLAY_ENV="FALSE",
               OMP_AFFINITY_FORMAT="thread=%n affinity=%A")
    for key in ("GOMP_CPU_AFFINITY", "GOMP_SPINCOUNT"):
        env.pop(key, None)
    with (work / "environment.log").open("w") as log:
        for command in (["date", "--iso-8601=seconds"], ["uptime"], ["free", "-h"],
                        ["mpstat", "-N", "ALL" if report else str(NUMA_NODE), "1", "1"]):
            subprocess.run(command, stdout=log, env=env, check=True)
        log.write(f"layout={'dual-numa-physical-cores' if report else 'single-numa-gate'}\n")
        log.write("OMP_PLACES=threads OMP_PROC_BIND=close OMP_WAIT_POLICY=PASSIVE OMP_DYNAMIC=FALSE\n")
    placement = verify_topology(work, env, cpu_masks, report)

    reference = (root / "refer-2dim-5h.txt").read_bytes()
    result = work / "result.txt"

    def run(threads, tag, affinity=False):
        setting = str(threads)
        run_env = dict(env, PIVOT_THREADS=setting, OMP_NUM_THREADS=setting, OMP_THREAD_LIMIT=setting,
                       OMP_DISPLAY_AFFINITY="TRUE" if affinity else "FALSE")
        command = ["timeout", "-k", "1s", "4s", "taskset", "-c", cpu_masks[threads],
                   str(root / "pivot"), str(root / "uniformvector-2dim-5h.txt")]
        if not affinity:
            command = ["/usr/bin/time", "-f", "%e\t%U\t%S\t%M", "-o", str(work / f"{tag}.time"), *command]
        result.unlink(missing_ok=True)
        proc = subprocess.run(command, env=run_env, cwd=work, capture_output=True)
        (work / f"{tag}.stdout").write_bytes(proc.stdout)
        (work / f"{tag}.stderr").write_bytes(proc.stderr)
        if proc.returncode:
            raise SystemExit(f"{tag} exited {proc.returncode}; see {work / (tag + '.stderr')}")
        if not result.is_file() or result.read_bytes() != reference:
            raise SystemExit(f"{tag}: result.txt differs from the official reference")
        return proc

    with (work / "affinity.tsv").open("w") as log:
        log.write("threads\tlogical_cpus\tnuma_node:thread_count\tphysical_cores\tthreads_per_core\n")
        for threads, mask in cpu_masks.items():
            if threads > 1:  # The serial path does not create an OpenMP team.
                proc = run(threads, f"affinity-{threads}", affinity=True)
                matches = [AFFINITY.fullmatch(line) for line in proc.stderr.splitlines()]
                rows = [tuple(map(int, m.groups())) for m in matches if m is not None]
                if (len(rows) != threads or len(matches) != threads or
                    {thread for thread, _ in rows} != set(range(threads)) or
                    {cpu for _, cpu in rows} != expand_mask(mask)):
                    raise SystemExit(f"{threads} threads: unexpected team affinity; see affinity-{threads}.stderr")
            nodes, cores, per_core = placement[threads]
            log.write(f"{threads}\t{mask}\t{nodes}\t{cores}\t{per_core}\n")

    samples = defaultdict(list)
    with (work / "timing.tsv").open("w", buffering=1) as log:
        log.write("repeat\tthreads\tlogical_cpus\tinternal_ms\twall_s\tuser_s\tsys_s\trss_kib\n")
        for repeat in range(1, REPEATS + 1):
            for threads, mask in cpu_masks.items():
                tag = f"t{threads}-{repeat}"
                proc = run(threads, tag)
                internal = TIMING.findall(proc.stdout)
                if proc.stderr or not internal:
                    raise SystemExit(f"{tag}: unexpected stderr or missing timing; see {work}")
                progress = TIMING.sub(b"", proc.stdout)
                if threads == 1:
                    serial_progress = progress
                elif progress != serial_progress:
                    raise SystemExit(f"{tag}: non-timing stdout differs from the serial run")
                wall, user, sys, rss = (work / f"{tag}.time").read_text().split()
                internal_ms = internal[-1].decode()
                log.write(f"{repeat}\t{threads}\t{mask}\t{internal_ms}\t{wall}\t{user}\t{sys}\t{rss}\n")
                samples[threads].append((float(internal_ms), float(wall), float(user) + float(sys), int(rss)))

    serial = median(row[0] for row in samples[1])
    summary = ["threads  internal_ms  wall_s  speedup_vs_1  cpu_s  peak_rss_KiB"]
    for threads, rows in samples.items():
        internal, wall, cpu, rss = zip(*rows)
        ms = median(internal)
        summary.append(f"{threads:7d}  {ms:11.3f}  {median(wall):6.2f}  {serial / ms:11.2f}x  "
                       f"{median(cpu):5.2f}  {max(rss):12d}")
    summary.append(f"PASS {len(cpu_masks) * REPEATS} official runs; all outputs and progress match, affinity verified.")
    text = "\n".join(summary)
    (work / "summary.log").write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", action="store_true", help="use distinct physical cores across both NUMA nodes")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent.parent
    logdir = os.environ.get("BENCH_THREADS_LOG_DIR")
    if logdir:
        work = Path(logdir).resolve()
    else:
        parent = root / "reports" / "work"
        parent.mkdir(parents=True, exist_ok=True)
        work = Path(tempfile.mkdtemp(prefix="report-threads." if args.report else "threads.", dir=parent))
    print(f"Thread matrix logs: {work}", flush=True)
    with (work.parent / ".thread-matrix.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise SystemExit("Another thread matrix is running in this checkout; retry later")
        benchmark(root, work, args.report)
