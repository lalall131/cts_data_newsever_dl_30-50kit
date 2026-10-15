"""Regression: real prefetch nested layout and failed-state cache recovery."""
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import threading
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
import prepare_layout as layout
import run_download_merge as d

source=(ROOT/'tests/test_download_execution.py').read_text()
fake=source.split("FAKE=r'''",1)[1].split("'''",1)[0]

with tempfile.TemporaryDirectory(prefix='cts-sra-path-') as temp:
    base=Path(temp);bin_dir=base/'bin';bin_dir.mkdir()
    for name in ['prefetch','vdb-validate','fasterq-dump']:
        p=bin_dir/name;p.write_text(fake);p.chmod(0o755)
    old_env=dict(os.environ)
    os.environ['PATH']=str(bin_dir)+os.pathsep+os.environ['PATH']
    calls=base/'prefetch_calls.txt'
    os.environ['MOCK_PREFETCH_CALLS']=str(calls)
    os.environ['MOCK_PREFETCH_NESTED']='1'
    site=layout.read_json(ROOT/'config/site.json')
    site.update(data_root=str(base/'data'),download_attempts=1,min_free_gb=0)
    registry=layout.read_json(ROOT/'config/kit_registry.json')
    refs=layout.read_json(ROOT/'config/reference_registry.json')
    def plan(run):
        row={'目标Kit':'SPI-ISR','归一技术':'TT-seq','Run ID':run,'Source data ID':'GSM'+run[3:],
             'Species':'Homo sapiens','官方物种':'Homo sapiens','官方平台':'ILLUMINA',
             '官方布局':'PAIRED','样本表内Run数':'1'}
        item=layout.plan_run(row,site,registry,refs)
        layout.materialize([item],layout.build_libraries([item],site,registry),site)
        return item
    slots=threading.Semaphore(1)

    # The exact reported server case: nested SRA exists, failed state, no FASTQ.
    old=plan('SRR900009');cache=Path(old['paths']['sra_accession_dir'])
    payload=cache/'SRR900009'/'SRR900009.sra'
    payload.parent.mkdir(parents=True);payload.write_text(json.dumps({'run':'SRR900009'}))
    before=payload.read_bytes()
    state=d.read(old['paths']['state']);state['download']['status']='FAILED'
    state['error']={'message':'Expected one SRA payload'};d.save(old['paths']['state'],state)
    result=d.process_run(old,site,'all',slots)
    assert result['status']=='RUN_INPUT_READY',result
    assert not calls.exists(),'a downloaded nested SRA must not be downloaded again'
    assert payload.read_bytes()==before
    state=d.read(old['paths']['state'])
    assert state['download']['vdb_validate_pass'] and state['validation']['gzip']
    assert 'error' not in state,'historical path error must not remain as a current failure'
    assert state['download']['path']==str(payload)

    # New prefetch downloads may also append the accession to the requested -O.
    fresh=plan('SRR900008')
    assert d.process_run(fresh,site,'all',slots)['status']=='RUN_INPUT_READY'
    assert calls.read_text().splitlines()==['SRR900008']
    assert d.process_run(fresh,site,'all',slots)['status']=='RUN_INPUT_READY_REUSED'
    assert calls.read_text().splitlines()==['SRR900008']

    # Two possible payloads are ambiguous; do not silently choose or download.
    ambiguous=plan('SRR900002');cache=Path(ambiguous['paths']['sra_accession_dir'])
    (cache/'SRR900002').mkdir(parents=True)
    for path in [cache/'SRR900002.sra',cache/'SRR900002'/'SRR900002.sra']:
        path.write_text(json.dumps({'run':'SRR900002'}))
    result=d.process_run(ambiguous,site,'all',slots)
    assert result['status']=='FAILED_OR_NEEDS_REVIEW' and 'Ambiguous' in result['error'],result
    assert len(list(cache.rglob('*.sra')))==2 and not Path(ambiguous['paths']['ready']).exists()
    assert calls.read_text().splitlines()==['SRR900008']

    # Finding a file alone is insufficient: existing SRA must pass vdb-validate.
    invalid=plan('SRR900010');cache=Path(invalid['paths']['sra_accession_dir'])
    (cache/'SRR900010').mkdir(parents=True)
    (cache/'SRR900010'/'SRR900010.sra').write_bytes(b'bad mock SRA')
    result=d.process_run(invalid,site,'all',slots)
    assert result['status']=='FAILED_OR_NEEDS_REVIEW' and 'vdb-validate' in result['error'],result
    assert list(cache.parent.glob(cache.name+'.validation_failed_*'))
    assert not Path(invalid['paths']['ready']).exists()

    # Preserve support for direct .sralite files; ignore unrelated accessions.
    direct=base/'direct';direct.mkdir();lite=direct/'SRR900020.sralite';lite.write_bytes(b'x')
    assert d.sra_payload(direct,'SRR900020')==lite
    try:d.sra_payload(direct,'SRR900021')
    except FileNotFoundError:pass
    else:raise AssertionError('must not select a different accession')
    escape=base/'escape';escape.mkdir()
    (escape/'SRR900020.sra').symlink_to(lite)
    try:d.sra_payload(escape,'SRR900020')
    except RuntimeError as error:assert 'outside' in str(error)
    else:raise AssertionError('must reject an out-of-cache symlink')
    os.environ.clear();os.environ.update(old_env)

print('PASS: nested SRA recovery without prefetch, vdb validation before conversion, actual payload parent, new nested download, resume, ambiguity and quarantine guards')
