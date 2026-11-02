#!/usr/bin/env python3
"""Generate complete legal pointwise source rows; makes no rank assertion."""
from pathlib import Path
import argparse, datetime, fcntl, hashlib, json, math, os, shutil, signal, subprocess, sys, time

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'prior_restored/latest/finite'
EXPECTED = {
    'phi_stream_selected.cpp': 'bc0724e42b707c7a269c70b89299013c0cbeb35751f982ba82b9834fb6cdb8d2',
    'phi_stream_partitioned.cpp': '011f27ee78cc2a8ba105429c5b48ab18e58da6a1e538c1e04ff2df87f2e4d851',
}

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()

def save(path, value):
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2) + '\n')
    temp.replace(path)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--d', type=int, required=True)
    ap.add_argument('--threads', type=int, default=1)
    args = ap.parse_args()
    d = args.d
    assert d >= 167 and all(d % q for q in range(2, math.isqrt(d) + 1))
    assert args.threads in [1, 2]
    folder = ROOT / 'new_pointwise' / f'd{d}'
    folder.mkdir(parents=True, exist_ok=True)
    lock = (folder / 'pipeline.lock').open('a+')
    fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    for name, digest in EXPECTED.items():
        assert sha(SOURCE / name) == digest
        if not (folder / name).exists():
            shutil.copyfile(SOURCE / name, folder / name)
        assert sha(folder / name) == digest
    binding = {'d': d, 'N': 4*d, 'source_code': EXPECTED, 'driver_sha256': sha(__file__),
               'purpose': 'Pointwise source preparation only; no uniform d criterion or finite interval asserted.'}
    bp = folder / 'SOURCE_BINDINGS.json'
    if bp.exists():
        assert json.loads(bp.read_text()) == binding
    else:
        save(bp, binding)
    child = None
    started = time.monotonic()
    def event(stage, **extra):
        record = {'d': d, 'stage': stage, 'utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  'seconds': time.monotonic() - started, 'source_bindings_sha256': sha(bp), **extra}
        save(folder / 'STATE.json', record)
        with (folder / 'EVENTS.jsonl').open('a') as f:
            f.write(json.dumps(record) + '\n')
        print(json.dumps(record), flush=True)
    def interrupted(sig, frame):
        if child is not None and child.poll() is None:
            os.killpg(child.pid, signal.SIGTERM)
        event('INTERRUPTED', signal=sig)
        raise SystemExit(128 + sig)
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    def run(stage, command, log):
        nonlocal child
        event(stage, command=command)
        with (folder / log).open('a') as f:
            child = subprocess.Popen(command, cwd=folder, stdout=f, stderr=subprocess.STDOUT, start_new_session=True)
            rc = child.wait()
        child = None
        event(stage + '_RETURNED', returncode=rc)
        assert rc == 0
    binary = folder / 'phi_stream_partitioned'
    if not binary.exists():
        run('COMPILE', ['g++', '-O3', '-std=c++17', '-fopenmp', 'phi_stream_partitioned.cpp', '-lgmpxx', '-lgmp', '-o', str(binary)], 'compile.log')
        save(folder / 'BINARY.json', {'sha256': sha(binary), 'source_sha256': EXPECTED['phi_stream_partitioned.cpp']})
    assert sha(binary) == json.loads((folder / 'BINARY.json').read_text())['sha256']
    raw = folder / f'd{d}_prefix_rows.jsonl'
    expected = {(k*d, m) for k in range(1, 5) for m in range(k*(d-1)//2) if k < 4 or m <= d}
    assert len(expected) == 4*d-2
    def keys():
        seen = set()
        if raw.exists():
            with raw.open() as f:
                for line in f:
                    assert line.endswith('\n'), 'Preserve and explicitly repair partial last row before resuming.'
                    row = json.loads(line)
                    key = row['n'], row['m']
                    assert key in expected and key not in seen
                    assert row['d'] == d and row['nd'] == 4*d and len(row['coefficients']) == row['n']+1
                    seen.add(key)
        return seen
    if keys() != expected:
        run('GENERATE_COMPLETE_LEGAL_PREFIX', [str(binary), str(d), str(raw), str(args.threads), str(4*d), 'prefix'], 'source_generation.log')
    assert keys() == expected
    event('COMPLETE_SOURCE_ONLY_NO_RANK_CLAIM', rows=len(expected), source_sha256=sha(raw), bytes=raw.stat().st_size)

if __name__ == '__main__':
    main()
