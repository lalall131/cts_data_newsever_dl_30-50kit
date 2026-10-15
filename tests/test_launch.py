"""Exercise the Linux entrypoint and local installer without network or real data."""
import csv
import hashlib
import io
import json
import os
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path

root=Path(__file__).resolve().parents[1]
source=(root/'tests/test_download_execution.py').read_text()
fake=source.split("FAKE=r'''",1)[1].split("'''",1)[0]

with tempfile.TemporaryDirectory(prefix='cts-launch-test-') as temp:
    base=Path(temp);pkg=base/'package';pkg.mkdir()
    for name in ['tools','config']:
        shutil.copytree(root/name,pkg/name,ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copy2(root/'start.sh',pkg/'start.sh')
    (pkg/'manifests').mkdir()
    site=json.loads((pkg/'config/site.json').read_text())
    site.update(data_root=str(base/'data'),min_free_gb=0,download_attempts=1)
    (pkg/'config/site.json').write_text(json.dumps(site))
    rows=[{'目标Kit':'SPI-ISR','归一技术':'TT-seq','Run ID':run,'Source data ID':lib,
           'Species':'Homo sapiens','官方物种':'Homo sapiens','官方平台':'ILLUMINA',
           '官方布局':layout,'样本表内Run数':'1'}
          for run,lib,layout in [('SRR900001','GSM900001','SINGLE'),('SRR900004','GSM900004','PAIRED')]]
    with (pkg/'manifests/selected_runs.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    (pkg/'config/library_decisions.json').write_text(json.dumps({r['Source data ID']:
        {'action':'SINGLE_RUN','approved':True,'included_runs':[r['Run ID']]} for r in rows}))
    bin_dir=base/'bin';bin_dir.mkdir()
    for name in ['prefetch','vdb-validate','fasterq-dump']:
        p=bin_dir/name;p.write_text(fake);p.chmod(0o755)
    env=dict(os.environ,PATH=str(bin_dir)+os.pathsep+os.environ['PATH'])
    def launch(*args):
        r=subprocess.run(['bash','start.sh',*args],cwd=pkg,env=env,capture_output=True,text=True)
        assert r.returncode==0,r.stdout+'\n'+r.stderr
        return r
    launch('--check')
    assert not (base/'data').exists(),'--check must not materialize raw data directories'
    launch('--stage','all','--workers','2')
    result=json.loads((pkg/'generated/execution_summary.json').read_text())
    assert result['run_statuses']=={'RUN_INPUT_READY':2},result
    assert result['library_statuses']=={'LIBRARY_INPUT_READY':2},result
    assert result['calculation_jobs_submitted']==0
    assert (pkg/'generated/run_storage.tsv').is_file()
    launch('--stage','all')
    result=json.loads((pkg/'generated/execution_summary.json').read_text())
    assert result['run_statuses']=={'RUN_INPUT_READY_REUSED':2},result

    # Install mock official-shaped binaries through the actual bootstrap path.
    import sys
    sys.path.insert(0,str(root/'tools'))
    import bootstrap_sra as boot
    payload=io.BytesIO()
    with tarfile.open(fileobj=payload,mode='w:gz') as archive:
        for name in ['prefetch','vdb-validate','fasterq-dump']:
            data=b'#!/bin/sh\nprintf "mock-tool 3.0.0\\n"\n'
            item=tarfile.TarInfo('sratoolkit.mock-ubuntu64/bin/'+name)
            item.size=len(data);item.mode=0o755;archive.addfile(item,io.BytesIO(data))
    raw=payload.getvalue();original=boot.urllib.request.urlopen
    boot.ROOT=pkg
    boot.urllib.request.urlopen=lambda *a,**k:io.BytesIO(raw)
    try:boot.main()
    finally:boot.urllib.request.urlopen=original
    receipt=json.loads((pkg/'runtime/runtime_bin.json').read_text())
    assert receipt['archive_sha256']==hashlib.sha256(raw).hexdigest()
    assert all('3.0.0' in v for v in receipt['versions'].values())
    assert Path(receipt['bin_dir']).is_dir()

print('PASS: Linux start.sh preflight, two-Run execution, library readiness, inventory, resumable rerun, local toolkit extraction/receipt; no real network data downloaded')
