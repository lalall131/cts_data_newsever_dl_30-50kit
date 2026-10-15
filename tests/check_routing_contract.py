import copy
import csv
import importlib.util
import json
import tempfile
from pathlib import Path

root=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('layout',root/'tools/prepare_layout.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
site=m.read_json(root/'config/site.json')
registry=m.read_json(root/'config/kit_registry.json')
references=m.read_json(root/'config/reference_registry.json')
with (root/'manifests/selected_runs.tsv').open(encoding='utf-8-sig',newline='') as f:
    rows=list(csv.DictReader(f,delimiter='\t'))
plans=[m.plan_run(r,site,registry,references) for r in rows]
libs=m.build_libraries(plans,site,registry)
assert len(plans)==339 and len(libs)==278
assert len({p['paths']['run_root'] for p in plans})==339
added=[p for p in plans if p.get('selection_batch')=='added_three_kits_20261014']
assert len(added)==76
assert all(p['download']['mode']=='ENA_ARCHIVE_FASTQ' and p['species']['tax_id']=='10090' for p in added)
assert all('SPECIAL_INPUT_ADAPTER_UNCONFIRMED' in p['routing']['hold_reasons'] for p in added)
assert len([p for p in added if 'TARGETED_HBB_CUSTOM_REFERENCE_REQUIRED' in p['routing']['hold_reasons']])==6
for p in plans:
    parts=Path(p['paths']['run_root']).parts
    assert parts[-1]==p['run_id']
    assert p['library_id'] not in parts[-3:]
    assert p['routing']['route_status']=='HELD'
    assert not p['routing']['calculation_submission_allowed']
assert all(len(l['run_ids'])==1 or 'MULTI_RUN_MERGE_NEEDS_REVIEW' in l['hold_reasons'] for l in libs)
assert sum(p['download']['mode']=='GSA_MANUAL_CHANNEL' for p in plans)==1

with tempfile.TemporaryDirectory(prefix='cts-routing-') as temp:
    base=Path(temp)
    kit=copy.deepcopy(registry)
    active=kit['kits']['SPI-ISR']
    active.update(enabled=True,version='test-only',entrypoint=str(base/'entry.sh'),
                  adapter_config=str(base/'adapter.json'),supported_presets=['ISR_default'])
    (base/'entry.sh').write_text('#!/bin/sh\nexit 0\n')
    (base/'adapter.json').write_text('{}')
    assets={}
    for name,kind in active['required_assets'].items():
        p=base/name
        assets[name]=str(p)
        if kind=='directory':p.mkdir()
        elif kind=='file':p.write_text('fixture\n')
        elif kind=='intron_prefix':
            for ending in ['.bed','_part1.bed','_part2.bed']:Path(str(p)+ending).write_text('fixture\n')
        elif kind=='bowtie2_prefix':
            for ending in ['1','2','3','4','rev.1','rev.2']:Path(str(p)+'.'+ending+'.bt2').write_text('fixture\n')
    ref={'defaults':{'9606':'human_fixture','10090':'mouse_fixture'},'bundles':{}}
    for tax,refid in [('9606','human_fixture'),('10090','mouse_fixture')]:
        ref['bundles'][refid]={'reference_id':refid,'tax_id':tax,'validated':True,
                             'assembly_accession':'fixture-only','annotation_release':'fixture-only',
                             'compatible_kits':['SPI-ISR'],'assets':assets}
    row=copy.deepcopy(rows[0]);row.update(biological_review_status='CONFIRMED',biological_review_evidence_url='fixture://review')
    valid=m.plan_run(row,site,kit,ref,check_assets=True)
    assert valid['routing']['route_status']=='CONFIGURATION_READY'
    assert valid['routing']['reference_id']=='human_fixture'
    assert not valid['routing']['calculation_submission_allowed']
    mouse=copy.deepcopy(row);mouse.update(Species='Mus musculus',官方物种='Mus musculus')
    assert m.plan_run(mouse,site,kit,ref,True)['routing']['reference_id']=='mouse_fixture'
    inferred=copy.deepcopy(row);inferred['目标Kit']=''
    assert m.plan_run(inferred,site,kit,ref,True)['routing']['kit_id']=='SPI-ISR'
    for platform,expected in [('OXFORD_NANOPORE','SPI-NPLR'),('PACBIO_SMRT','SPI-PBLR')]:
        flep=copy.deepcopy(row);flep.update(目标Kit='',归一技术='FLEP-seq',官方平台=platform)
        assert m.plan_run(flep,site,kit,ref,True)['routing']['kit_id']==expected
    wrong=copy.deepcopy(row);wrong['官方平台']='OXFORD_NANOPORE'
    assert 'PLATFORM_KIT_CONFLICT' in m.plan_run(wrong,site,kit,ref,True)['routing']['hold_reasons']
    refwrong=copy.deepcopy(ref);refwrong['bundles']['human_fixture']['tax_id']='10090'
    assert 'REFERENCE_TAXON_MISMATCH' in m.plan_run(row,site,kit,refwrong,True)['routing']['hold_reasons']
    missing=copy.deepcopy(ref);missing['bundles']['human_fixture']['assets']['genes_bed']=str(base/'missing.bed')
    assert 'REFERENCE_FILE_MISSING:genes_bed' in m.plan_run(row,site,kit,missing,True)['routing']['hold_reasons']
    unsafe=copy.deepcopy(row);unsafe['Run ID']='../escape'
    try:m.plan_run(unsafe,site,kit,ref,True)
    except ValueError:pass
    else:raise AssertionError('Path traversal was accepted')

    # Materialization preserves active state and writes refreshed routing separately.
    testsite=copy.deepcopy(site);testsite['data_root']=(base/'data').as_posix()
    p=m.plan_run(row,testsite,registry,references)
    ll=m.build_libraries([p],testsite,registry)
    m.materialize([p],ll,testsite)
    state=Path(p['paths']['state'])
    assert json.loads(state.read_text())['download']['status']=='PLANNED'
    state.write_text(json.dumps({'download':'RUNNING'}))
    m.materialize([p],ll,testsite)
    assert json.loads(state.read_text())['download']=='RUNNING'
    assert Path(p['paths']['routing']).is_file()
    assert Path(ll[0]['library_root'])==Path(p['paths']['run_root'])
print('PASS: 339 Run paths / 278 sample groups; 76 exact mouse FASTQ sources, special adapter/targeted reference holds, reference routing, platform conflicts, missing assets, traversal and state preservation')
