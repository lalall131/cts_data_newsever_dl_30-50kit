import copy
import importlib.util
import json
import os
import sys
import tempfile
import threading
from pathlib import Path

root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'tools'))
import prepare_layout as layout
import run_download_merge as d

FAKE=r'''#!/usr/bin/env python3
import json,sys
from pathlib import Path
name=Path(sys.argv[0]).name;args=sys.argv[1:]
if '--version' in args:print(name+' 3.0.0 mock-test-only');raise SystemExit(0)
if name=='prefetch':
    run=args[0];cache=Path(args[args.index('-O')+1]);cache.mkdir(parents=True,exist_ok=True)
    assert cache.name==run
    (cache/(run+'.sra')).write_text(json.dumps({'run':run}))
elif name=='vdb-validate':
    p=Path(args[0]);json.loads(p.read_text());assert p.suffix=='.sra'
elif name=='fasterq-dump':
    cache=Path(args[0]);assert cache.is_dir(),'fasterq must receive accession directory'
    run=cache.name;out=Path(args[args.index('--outdir')+1]);out.mkdir(parents=True,exist_ok=True)
    assert args[args.index('--seq-defline')+1]=='@$ac.$si/$ri'
    paired=run not in ['SRR900001','SRR900007']
    suffixes=['_1','_2'] if paired else ['']
    if run=='SRR900003':suffixes+=['']
    for suffix in suffixes:
        role=2 if suffix=='_2' else 1
        seq='ACGT' if int(run[-1])%2 else 'TGCA'
        if role==2:seq=seq[::-1]
        quality='III' if run=='SRR900007' else 'IIII'
        n=1 if run=='SRR900003' and suffix=='' else 2
        with (out/(run+suffix+'.fastq')).open('w') as f:
            for i in range(1,n+1):f.write('@'+run+'.'+str(i)+'/'+str(role)+'\n'+seq+'\n+\n'+quality+'\n')
else:raise SystemExit(1)
'''

with tempfile.TemporaryDirectory(prefix='cts-exec-test-') as temp:
    base=Path(temp);bin_dir=base/'bin';bin_dir.mkdir()
    for name in ['prefetch','vdb-validate','fasterq-dump']:
        p=bin_dir/name;p.write_text(FAKE);p.chmod(0o755)
    old_path=os.environ.get('PATH','');os.environ['PATH']=str(bin_dir)+os.pathsep+old_path
    site=layout.read_json(root/'config/site.json');site.update(data_root=str(base/'data'),download_attempts=1,min_free_gb=0)
    registry=layout.read_json(root/'config/kit_registry.json');refs=layout.read_json(root/'config/reference_registry.json')
    slots=threading.Semaphore(1)
    def make(run,library,layout_name='PAIRED',count=1):
        row={'目标Kit':'SPI-ISR','归一技术':'TT-seq','Run ID':run,'Source data ID':library,
             'Species':'Homo sapiens','官方物种':'Homo sapiens','官方平台':'ILLUMINA','官方布局':layout_name,
             '样本表内Run数':str(count)}
        p=layout.plan_run(row,site,registry,refs);ll=layout.build_libraries([p],site,registry)
        layout.materialize([p],ll,site)
        return p
    one=make('SRR900001','GSM900001','SINGLE')
    result=d.process_run(one,site,'all',slots);assert result['status']=='RUN_INPUT_READY',result
    assert (Path(one['paths']['merged'])/'SRR900001.fastq.gz').is_symlink()
    assert not list(Path(one['paths']['tmp']).glob('*/raw/*.fastq'))
    result=d.process_run(one,site,'all',slots);assert result['status']=='RUN_INPUT_READY_REUSED',result
    raw_file=Path(d.read(one['paths']['state'])['validation']['input_files'][0]['path'])
    raw_file.write_bytes(b'damaged')
    result=d.process_run(one,site,'all',slots);assert result['status']=='RUN_INPUT_READY',result
    assert list(raw_file.parent.glob('*.saved_*'))
    extra=make('SRR900003','GSM900003')
    result=d.process_run(extra,site,'all',slots)
    assert result['status']=='FAILED_OR_NEEDS_REVIEW' and 'READ_ROLE_REVIEW' in result['error']
    assert len(list(Path(extra['paths']['fastq']).glob('*.fastq.gz')))==3
    assert not Path(extra['paths']['ready']).exists()
    bad=make('SRR900007','GSM900007','SINGLE')
    result=d.process_run(bad,site,'all',slots)
    assert result['status']=='FAILED_OR_NEEDS_REVIEW' and 'length mismatch' in result['error']
    assert not Path(bad['paths']['ready']).exists()
    a=make('SRR900004','GSM900100',count=2);b=make('SRR900005','GSM900100',count=2)
    for p in [a,b]:assert d.process_run(p,site,'all',slots)['status']=='RUN_INPUT_READY'
    lib=layout.build_libraries([a,b],site,registry)[0]
    layout.materialize([a,b],[lib],site)
    mapping={p['run_id']:p for p in [a,b]}
    assert d.merge_library(lib,mapping,{'action':'HOLD'},site)['status']=='MERGE_REVIEW_REQUIRED'
    decision={'action':'CONCAT','approved':True,'included_runs':lib['run_ids'],'reason':'test-only distinct Runs'}
    result=d.merge_library(lib,mapping,decision,site);assert result['status']=='LIBRARY_INPUT_READY',result
    ready=d.read(Path(lib['library_root'])/'library_ready.json')
    assert all(r['records']==4 for r in ready['input_files'])
    assert not ready['calculation_submission_allowed']
    result=d.merge_library(lib,mapping,decision,site);assert result['status']=='LIBRARY_INPUT_READY_REUSED',result
    subset={'action':'USE_RUNS','approved':True,'included_runs':[a['run_id']],'reason':'test-only duplicate exclusion decision'}
    result=d.merge_library(lib,mapping,subset,site);assert result['status']=='LIBRARY_INPUT_READY',result
    ready=d.read(Path(lib['library_root'])/'library_ready.json');assert ready['excluded_runs']==[b['run_id']]
    dup=make('SRR900006','GSM900100',count=2)
    assert d.process_run(dup,site,'all',slots)['status']=='RUN_INPUT_READY'
    dl=layout.build_libraries([a,dup],site,registry)[0];layout.materialize([a,dup],[dl],site)
    result=d.merge_library(dl,{a['run_id']:a,dup['run_id']:dup},{'action':'CONCAT','approved':True,'included_runs':dl['run_ids'],'reason':'test duplicate guard'},site)
    assert result['status']=='FAILED_OR_NEEDS_REVIEW' and 'DUPLICATE_RUN_SEQUENCE_CONTENT' in result['error']
    gsa=make('CRR900001','CRX900001')
    assert d.process_run(gsa,site,'all',slots)['status']=='MANUAL_SOURCE_REQUIRED'
    os.environ['PATH']=old_path
print('PASS: mock CLI download/validation/conversion, resume, corrupt FASTQ replacement, orphan preservation, invalid format hold, approved merge, explicit exclusion, duplicate detection, GSA hold')
