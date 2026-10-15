import copy,csv,gzip,hashlib,json,os,subprocess,sys,tempfile,threading
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
import prepare_layout as layout
import run_download_merge as d

FAKE='''#!/usr/bin/env python3
import os,sys
from pathlib import Path
if '--version' in sys.argv:print('curl 8.0.0 mock');raise SystemExit(0)
args=sys.argv[1:];dest=Path(args[args.index('--output')+1]);name=args[-1].rsplit('/',1)[-1]
source=Path(os.environ['MOCKENA_DIR'])/name
begin=dest.stat().st_size if dest.exists() else 0
with dest.open('ab') as f:f.write(source.read_bytes()[begin:])
'''
with tempfile.TemporaryDirectory(prefix='cts-archive-test-') as temp:
    base=Path(temp);bin_dir=base/'bin';bin_dir.mkdir();(bin_dir/'curl').write_text(FAKE);(bin_dir/'curl').chmod(0o755)
    sources=base/'sources';sources.mkdir()
    oldpath=os.environ.get('PATH','');os.environ['PATH']=str(bin_dir)+os.pathsep+oldpath
    os.environ['MOCKENA_DIR']=str(sources)
    site=layout.read_json(ROOT/'config/site.json');site.update(data_root=str(base/'data'),min_free_gb=0,download_attempts=1)
    registry=layout.read_json(ROOT/'config/kit_registry.json');refs=layout.read_json(ROOT/'config/reference_registry.json')
    row={'目标Kit':'scFLUENT-seq SPI','归一技术':'scFLUENT-seq','Run ID':'SRR900001','Source data ID':'GSM900001',
         'Species':'Mus musculus','官方物种':'Mus musculus','官方平台':'ILLUMINA','官方布局':'PAIRED','样本表内Run数':'1',
         'download_execution_mode':'ENA_ARCHIVE_FASTQ','selection_batch':'added_three_kits_20261014',
         'special_input_profile':'test_barcode','special_calculation_hold':'test dedicated input adapter required'}
    original={}
    for role in [1,2]:
        name=f'SRR900001_{role}.fastq.gz'
        content=f'@instrument:1:2:3:4/ {role}:N:0:TAG\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n'.encode()
        # Original complete headers must survive; only first tokens are used for pairing checks.
        original[name]=gzip.compress(content,mtime=0);(sources/name).write_bytes(original[name])
    row['ENA_FASTQ地址']=';'.join('ftp.sra.ebi.ac.uk/test/'+name for name in original)
    row['ENA_FASTQ字节数']=';'.join(str(len(b)) for b in original.values())
    row['ENA_FASTQ_MD5']=';'.join(hashlib.md5(b).hexdigest() for b in original.values())
    plan=layout.plan_run(row,site,registry,refs);libs=layout.build_libraries([plan],site,registry)
    layout.materialize([plan],libs,site)
    # Start from a partial prefix; ensure the synchronous transfer resumes correctly.
    Path(plan['download']['ena_files'][0]['destination']+'.part').write_bytes(list(original.values())[0][:12])
    result=d.process_run(plan,site,'all',threading.Semaphore(1));assert result['status']=='ARCHIVE_FASTQ_VALIDATED_KIT_PREPARATION_REQUIRED',result
    for item in plan['download']['ena_files']:
        assert Path(item['destination']).read_bytes()==original[Path(item['destination']).name]
    root=Path(plan['paths']['run_root'])
    assert (root/'raw_fastq_manifest.json').is_file() and not (root/'ready.json').exists()
    assert not (root/'library_ready.json').exists()
    result=d.process_run(plan,site,'all',threading.Semaphore(1));assert 'REUSED' in result['status']
    result=d.merge_library(libs[0],{plan['run_id']:plan},None,site)
    assert result['status']=='KIT_SPECIFIC_PREPARATION_REQUIRED'
    # Existing corrupt final file must be preserved and replaced with official bytes.
    dest=Path(plan['download']['ena_files'][0]['destination']);dest.write_bytes(b'corrupted')
    result=d.process_run(plan,site,'all',threading.Semaphore(1));assert result['status']=='ARCHIVE_FASTQ_VALIDATED_KIT_PREPARATION_REQUIRED',result
    assert list(dest.parent.glob('*.unverified_saved_*'))
    # Wrong official digest must not publish a ready manifest or final file.
    bad=copy.deepcopy(row);bad['Run ID']='SRR900002';bad['Source data ID']='GSM900002';bad['ENA_FASTQ_MD5']='0'*32+';'+'0'*32
    badplan=layout.plan_run(bad,site,registry,refs);badlibs=layout.build_libraries([badplan],site,registry)
    layout.materialize([badplan],badlibs,site)
    result=d.process_run(badplan,site,'all',threading.Semaphore(1))
    assert result['status']=='FAILED_OR_NEEDS_REVIEW' and 'MD5 mismatch' in result['error'],result
    assert list(Path(badplan['paths']['fastq']).glob('*.md5_mismatch_saved_*'))
    os.environ['PATH']=oldpath
print('PASS: exact ENA download, partial-prefix resume, original headers/bytes, official MD5/full gzip/FASTQ, replacement quarantine, bad-MD5 hold, no false barcode/FLNC/library readiness')
