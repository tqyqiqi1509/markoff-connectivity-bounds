#!/usr/bin/env python3
"""Prepare, run, and inspect two root-owned finite pointwise pipelines.

Preparation compiles/copies audited code but executes no determinant. Only an
explicit `run --worker W` runs determinants. A missing source gives a resumable
WAITING_SOURCE state. No command here asserts a complete prime interval.
"""
from pathlib import Path
import argparse, collections, datetime, fcntl, hashlib, json, math, os
import shutil, signal, subprocess, sys, time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OLD = ROOT/'agents/attachment_audit/finite52/batch_000_d139'
DATA = HERE/'CRT57_1E65_INDEPENDENT_CLASSIFICATION.json'
PLAN = HERE/'BATCH_PLAN.json'
EXPECTED = {
    'run_manifest_group_verification.py': '286dac5a325ffe9182fade7ff1732b80cdfeded61ddebf8599df493e237853eb',
    'independent_modular_det.cpp': '5d6ea13aaeaf059fd38f7e953e95a105907a733ecbb9e0414199928c8e0e7d31',
    'phi_stream_selected.cpp': 'bc0724e42b707c7a269c70b89299013c0cbeb35751f982ba82b9834fb6cdb8d2',
    'phi_stream_partitioned.cpp': '011f27ee78cc2a8ba105429c5b48ab18e58da6a1e538c1e04ff2df87f2e4d851',
}


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    temp = path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2)+'\n')
    temp.replace(path)


def immutable(path, value):
    if path.exists():
        assert load(path) == value, ('Immutable record differs', str(path))
    else:
        save(path, value)


def source_ready(d):
    directory = ROOT/'new_pointwise'/f'd{d}'
    state_path = directory/'STATE.json'
    state = load(state_path)
    if state['stage'] != 'COMPLETE_SOURCE_ONLY_NO_RANK_CLAIM':
        return None
    assert state['d'] == d and state['rows'] == 4*d-2
    source = directory/f'd{d}_prefix_rows.jsonl'
    assert sha(source) == state['source_sha256'] and source.stat().st_size == state['bytes']
    bp = directory/'SOURCE_BINDINGS.json'
    assert sha(bp) == state['source_bindings_sha256']
    binding = load(bp)
    assert binding['d'] == d and binding['N'] == 4*d
    for name in ('phi_stream_selected.cpp', 'phi_stream_partitioned.cpp'):
        assert sha(directory/name) == binding['source_code'][name] == EXPECTED[name]
    driver = ROOT/'control/run_pointwise_prefix.py'
    assert sha(driver) == binding['driver_sha256']
    binary_record = load(directory/'BINARY.json')
    binary = directory/'phi_stream_partitioned'
    assert sha(binary) == binary_record['sha256']
    assert binary_record['source_sha256'] == EXPECTED['phi_stream_partitioned.cpp']
    # The legacy wrapper records selected.cpp as the mathematical source;
    # this separate receipt records the actual compiled/executed entry point.
    receipt = {'status': 'COMPLETE_SOURCE_IDENTITY_BOUND_NOT_A_RANK_CLAIM', 'd': d, 'N': 4*d,
               'source_path': str(source.relative_to(ROOT)), 'source_sha256': sha(source),
               'source_bytes': source.stat().st_size,
               'actual_generator': 'phi_stream_partitioned.cpp',
               'generator_sha256': EXPECTED['phi_stream_partitioned.cpp'],
               'legacy_selected_generator_equivalent_prefix_sha256': EXPECTED['phi_stream_selected.cpp'],
               'scope': 'Source identity and completed row-key generation state. Full integer coefficient, moment and legal-row checks are performed by each audited pointwise wrapper; this receipt does not claim a second source generation.',
               'bindings': {str(p.relative_to(ROOT)): sha(p) for p in (state_path, bp, driver, directory/'BINARY.json', binary)}}
    immutable(HERE/f'd{d}_SOURCE_READY.json', receipt)
    return receipt


def prepare():
    assert not PLAN.exists(), 'Plan is immutable; do not regenerate over running/completed jobs.'
    data = load(DATA)
    assert data['status'] == 'CONDITIONAL_PENDING_CRITERIA_AND_POINTWISE_EVIDENCE'
    assert data['candidate_count'] == 2689 and data['pointwise_count'] == 159
    candidates = collections.defaultdict(list)
    for r in data['records']:
        if r['status'] == 'POINTWISE_EVIDENCE_REQUIRED':
            d, n = r['preferred_d'], int(r['n'])
            assert d in (167, 173) and n > 8*d and n % (2*d) not in (1, 2*d-1)
            candidates[d].append(r['n'])
    assert {d: len(v) for d, v in candidates.items()} == {167: 157, 173: 2}
    common = HERE/'common'
    common.mkdir(exist_ok=True)
    for name, expected in EXPECTED.items():
        src = (ROOT/'new_pointwise/d167'/name if name == 'phi_stream_partitioned.cpp' else OLD/name)
        assert sha(src) == expected
        shutil.copyfile(src, common/name)
        assert sha(common/name) == expected
    command = ['g++','-O3','-std=c++17',str(common/'independent_modular_det.cpp'),'-lgmpxx','-lgmp','-o',str(common/'independent_modular_det')]
    with (HERE/'compile_gmp.log').open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    code_sha = {name: sha(common/name) for name in list(EXPECTED)+['independent_modular_det']}
    binary_binding = load(ROOT/'control/FINITE_BINARY_ROOT_BINDING.json')
    assert code_sha['independent_modular_det'] == binary_binding['binary_sha256']
    immutable(HERE/'FRESH_BINARY_BINDING.json', {'status': 'FRESH_COMPILE_MATCHES_ROOT_AUDITED_BINARY',
              'command': command, 'binary_sha256': code_sha['independent_modular_det'],
              'source_sha256': EXPECTED['independent_modular_det.cpp'],
              'root_binary_binding_sha256': sha(ROOT/'control/FINITE_BINARY_ROOT_BINDING.json')})
    batches = []
    for d in (167, 173):
        for start in range(0, len(candidates[d]), 20):
            index = len(batches)
            vals = candidates[d][start:start+20]
            name = f'batch_{index:03d}_d{d}'
            directory = HERE/name
            directory.mkdir()
            for filename in code_sha:
                if filename == 'run_manifest_group_verification.py':
                    # Its BASE is Path(__file__).resolve().parent, therefore
                    # this Python entry point must be an actual batch file.
                    shutil.copyfile(common/filename, directory/filename)
                else:
                    (directory/filename).symlink_to(Path('..')/'common'/filename)
            raw = ROOT/'new_pointwise'/f'd{d}'/f'd{d}_prefix_rows.jsonl'
            (directory/f'd{d}_prefix_rows.jsonl').symlink_to(os.path.relpath(raw, directory))
            (directory/f'd{d}_candidates.txt').write_text('\n'.join(vals)+'\n')
            manifest = {'status': 'POINTWISE_BATCH_NO_INTERVAL_PROMOTION',
                        'groups': [{'d': d, 'N': 4*d, 'candidate_count': len(vals), 'candidates': vals}],
                        'classification_sha256': sha(DATA),
                        'source_ready_receipt': f'../d{d}_SOURCE_READY.json',
                        'source_identity_rule': 'Run only after the immutable complete source receipt exists and matches the current source bytes.',
                        'purpose': 'Conditional-primality connectivity or exact compositeness only; no interval or uniform167/173 claim.'}
            save(directory/'batch_manifest.json', manifest)
            batches.append({'id': index, 'directory': name, 'worker': index % 2, 'd': d,
                            'candidate_count': len(vals), 'candidates': vals,
                            'manifest_sha256': sha(directory/'batch_manifest.json')})
    immutable(PLAN, {'status': 'SCHEDULE_ONLY_NO_INTERVAL_PROMOTION', 'classification_sha256': sha(DATA),
              'planned_total': 159, 'workers': 2, 'max_batch_size': 20,
              'common_code_sha256': code_sha, 'campaign_code_sha256': sha(__file__),
              'batches': batches,
              'packaging': 'Keep each authoritative raw source once. Batch files are relative symlinks; restore those links, do not duplicate 135MB+ per batch.',
              'source_theory': 'Reviewed general prime-d finite-point rank criterion plus high-half/Frobenius/pure-Vieta bridge remain explicit dependencies.'})
    for d in (167, 173):
        receipt = source_ready(d)
        print(json.dumps({'d': d, 'source_ready': receipt is not None}), flush=True)
    print(json.dumps({'plan': str(PLAN), 'batches': len(batches),
                      'worker_counts': {w: sum(b['candidate_count'] for b in batches if b['worker'] == w) for w in (0,1)}}), flush=True)


def verify_batch(batch, plan):
    directory, d = HERE/batch['directory'], batch['d']
    rp = directory/f'd{d}_complete_verification.json'
    if not rp.exists():
        return None
    ready = load(HERE/f'd{d}_SOURCE_READY.json')
    assert sha(ROOT/ready['source_path']) == ready['source_sha256']
    assert sha(directory/'batch_manifest.json') == batch['manifest_sha256']
    for name, digest in plan['common_code_sha256'].items():
        assert sha(directory/name) == digest
    r = load(rp)
    assert r['status'] == 'COMPLETE' and r['d'] == d and r['N'] == 4*d
    assert r['source_sha256'] == ready['source_sha256']
    assert r['gmp_code_sha256'] == EXPECTED['independent_modular_det.cpp']
    assert r['python_code_sha256'] == EXPECTED['run_manifest_group_verification.py']
    audit_path = directory/f'd{d}_source_audit.json'
    audit = load(audit_path)
    assert audit['source_sha256'] == ready['source_sha256']
    assert audit['generator_sha256'] == EXPECTED['phi_stream_selected.cpp']
    assert audit['candidate_manifest_sha256'] == batch['manifest_sha256']
    assert audit['source_row_count'] == 4*d-2
    assert audit['all_rows_legal_and_complete'] and audit['all_three_moment_relations']
    expected_keys = [(k*d,m) for k in range(1,5) for m in range(k*(d-1)//2 if k<4 else d+1)]
    assert [(x['n'],x['m']) for x in audit['row_order']] == expected_keys
    assert [x['matrix_row'] for x in audit['row_order']] == list(range(4*d-2))
    assert sorted(x['source_line_zero_based'] for x in audit['row_order']) == list(range(4*d-2))
    gp = directory/f'd{d}_gmp_determinants.jsonl'
    xp = directory/f'd{d}_cross_verified_certificates.jsonl'
    gmp, cross = [], []
    for p, output in ((gp,gmp),(xp,cross)):
        with p.open() as f:
            for line in f:
                assert line.endswith('\n')
                output.append(json.loads(line))
    assert cross == r['records']
    assert len(cross) == len(gmp) == r['candidate_count'] == batch['candidate_count']
    assert [x['candidate'] for x in cross] == [x['candidate'] for x in gmp] == batch['candidates']
    counts = collections.Counter()
    for item, original in zip(cross,gmp):
        assert all(item[k] == v for k,v in original.items())
        n = int(item['candidate'])
        assert item['d'] == d and item['N'] == 4*d
        assert item['source_sha256'] == ready['source_sha256']
        assert n > 8*d and n % (2*d) not in (1,2*d-1)
        assert item['characteristic_gate'] and item['congruence_gate']
        status = item['accepted_status']
        if status == 'IF_PRIME_THEN_CONNECTED':
            det = int(item['determinant_mod_n'])
            assert item['status'] == 'nonzero_determinant' and 0 < det < n
            assert math.gcd(det,n) == 1 and det == int(item['flint_determinant_mod_n'])
            assert item['unit_pivots'] == 4*d-2 and item['independent_determinants_agree'] is True
        elif status == 'EXACT_COMPOSITE_FACTOR':
            factor = int(item['factor'])
            assert item['status'] == 'composite_factor' and 1 < factor < n and n % factor == 0
        else:
            assert status == 'UNRESOLVED_ZERO_SELECTED_MINOR'
            assert item['status'] == 'zero_determinant'
            assert int(item['determinant_mod_n']) == int(item['flint_determinant_mod_n']) == 0
        counts[status] += 1
    for key, status in (('conditional_connected','IF_PRIME_THEN_CONNECTED'),('composite','EXACT_COMPOSITE_FACTOR'),('unresolved','UNRESOLVED_ZERO_SELECTED_MINOR')):
        assert r[key] == counts[status]
    return {'batch': batch['directory'], 'count': len(cross), 'unresolved': r['unresolved'],
            'records': cross, 'bindings': {str(p.relative_to(ROOT)): sha(p) for p in (rp,gp,xp,audit_path,HERE/f'd{d}_SOURCE_READY.json')}}


def run(worker, resume=False, stop=False):
    stop_path = HERE/f'STOP_AFTER_BATCH_{worker}'
    if stop:
        stop_path.write_text(now()+'\n')
        print('Stop request saved; the current batch may finish.')
        return
    lock = (HERE/f'pipeline_{worker}.lock').open('a+')
    fcntl.flock(lock.fileno(), fcntl.LOCK_EX|fcntl.LOCK_NB)
    if resume and stop_path.exists():
        stop_path.unlink()
    plan = load(PLAN)
    assert plan['campaign_code_sha256'] == sha(__file__)
    assert plan['classification_sha256'] == sha(DATA)
    batches = [b for b in plan['batches'] if b['worker'] == worker]
    started, child, completed, waiting = time.monotonic(), None, [], []
    def event(stage, **extra):
        record = {'stage': stage, 'worker': worker, 'utc': now(), 'seconds': time.monotonic()-started,
                  'plan_sha256': sha(PLAN), 'completed_batches': completed, **extra}
        save(HERE/f'pipeline_{worker}.state.json',record)
        with (HERE/f'pipeline_{worker}.events.jsonl').open('a') as f:
            f.write(json.dumps(record)+'\n')
        print(json.dumps(record),flush=True)
    def interrupted(sig, frame):
        if child is not None and child.poll() is None:
            os.killpg(child.pid,signal.SIGTERM)
        event('INTERRUPTED',signal=sig)
        raise SystemExit(128+sig)
    signal.signal(signal.SIGTERM,interrupted)
    signal.signal(signal.SIGINT,interrupted)
    environment = os.environ.copy()
    deps = str(ROOT/'agents/attachment_audit/d58/deps')
    environment['PYTHONPATH'] = deps+(os.pathsep+environment['PYTHONPATH'] if environment.get('PYTHONPATH') else '')
    event('STARTED',candidate_count=sum(b['candidate_count'] for b in batches))
    for batch in batches:
        result = verify_batch(batch,plan)
        if result:
            completed.append(batch['directory'])
            event('RESUMED_COMPLETE',batch=batch['directory'])
            if result['unresolved']:
                event('STOPPED_UNRESOLVED',batch=batch['directory'])
                return
            continue
        if stop_path.exists():
            event('STOPPED_AFTER_BATCH',next_batch=batch['directory'])
            return
        d = batch['d']
        if source_ready(d) is None:
            waiting.append(batch['directory'])
            event('WAITING_SOURCE',batch=batch['directory'],d=d)
            continue
        directory = HERE/batch['directory']
        assert sha(directory/'batch_manifest.json') == batch['manifest_sha256']
        for name, digest in plan['common_code_sha256'].items():
            assert sha(directory/name) == digest
        # Preserve any incomplete prior attempt before the unmodified legacy
        # wrapper starts its full batch again. Completed receipts never move.
        outputs = [p for p in directory.iterdir() if p.name in {
            f'd{d}_source_audit.json', f'd{d}_gmp_determinants.jsonl', f'd{d}_gmp_determinants.log',
            f'd{d}_cross_verified_certificates.jsonl', f'd{d}_verification_state.json', 'run.log'}]
        if outputs:
            saved = directory/('incomplete_attempt_'+datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%f'))
            saved.mkdir()
            for p in outputs:
                p.rename(saved/p.name)
        event('BATCH_RUNNING',batch=batch['directory'],candidate_count=batch['candidate_count'])
        with (directory/'run.log').open('w') as log:
            child = subprocess.Popen([sys.executable,str(directory/'run_manifest_group_verification.py'),
                    str(d),'--manifest',str(directory/'batch_manifest.json')],cwd=directory,env=environment,
                    stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            rc = child.wait()
        child = None
        if rc:
            event('BATCH_FAILED',batch=batch['directory'],returncode=rc)
            raise SystemExit(rc)
        result = verify_batch(batch,plan)
        assert result is not None
        completed.append(batch['directory'])
        event('BATCH_COMPLETE',batch=batch['directory'],unresolved=result['unresolved'])
        if result['unresolved']:
            event('STOPPED_UNRESOLVED',batch=batch['directory'])
            return
    event('POINTWISE_PIPELINE_COMPLETE' if not waiting else 'PARTIAL_WAITING_SOURCE',waiting_batches=waiting,
          scope='Pointwise conditional evidence only; no full interval acceptance.')


def collect(require_complete):
    plan = load(PLAN)
    assert plan['campaign_code_sha256'] == sha(__file__)
    assert plan['classification_sha256'] == sha(DATA)
    results, bindings = {}, []
    for batch in plan['batches']:
        result = verify_batch(batch,plan)
        if result is None:
            continue
        bindings.append(result['bindings'])
        for r in result['records']:
            n = int(r['candidate'])
            assert n not in results
            results[n] = r
    expected = {int(n) for n in load(DATA)['pointwise_candidates']}
    assert set(results) <= expected
    missing = sorted(expected-set(results))
    unresolved = sorted(n for n,r in results.items() if r['accepted_status']=='UNRESOLVED_ZERO_SELECTED_MINOR')
    complete = not missing and not unresolved
    out = {'status':'ALL_159_POINTWISE_CERTIFICATES_COMPLETE' if complete else 'PARTIAL_POINTWISE_CHECKPOINT',
           'utc': now(),'interval_accepted':False,'scope':'General criteria and independent arithmetic closure still required.',
           'plan_sha256':sha(PLAN),'classification_sha256':sha(DATA),'completed_count':len(results),
           'expected_count':159,'missing_count':len(missing),'unresolved_count':len(unresolved),
           'missing':list(map(str,missing)),'unresolved':list(map(str,unresolved)),
           'counts':dict(collections.Counter(r['accepted_status'] for r in results.values())),
           'receipt_bindings':bindings,'records':[results[n] for n in sorted(results)]}
    save(HERE/'COLLECTED_POINTWISE_RESULT.json',out)
    print(json.dumps({k:v for k,v in out.items() if k not in {'missing','unresolved','records','receipt_bindings'}},indent=2))
    if require_complete:
        assert complete, 'Missing or unresolved pointwise evidence.'


def main():
    ap=argparse.ArgumentParser()
    sub=ap.add_subparsers(dest='command',required=True)
    sub.add_parser('prepare')
    p=sub.add_parser('run');p.add_argument('--worker',type=int,choices=(0,1),required=True)
    p.add_argument('--resume',action='store_true');p.add_argument('--stop',action='store_true')
    p=sub.add_parser('collect');p.add_argument('--require-complete',action='store_true')
    a=ap.parse_args()
    if a.command=='prepare':prepare()
    elif a.command=='run':run(a.worker,a.resume,a.stop)
    else:collect(a.require_complete)


if __name__=='__main__':main()
