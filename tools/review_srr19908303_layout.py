#!/usr/bin/env python3
"""Review all SRA spots for SRR19908303 and reuse its validated singleton FASTQ.

This tool never downloads or converts reads. It takes this Run's lock, verifies
the complete READ_LEN/READ_TYPE stream, checks the stored FASTQ SHA256 again,
and publishes SINGLE input manifests while retaining declared PAIRED metadata.
"""
import argparse
import copy
import fcntl
import subprocess
from collections import Counter
from pathlib import Path

import run_download_merge as d

RUN='SRR19908303'
LIBRARY='GSM6283344'


def scan_rows(lines,expected_records):
    """Require exactly one nonempty biological and one empty technical read."""
    count=0
    patterns=Counter()
    biological_positions=Counter()
    for raw in lines:
        line=raw.rstrip('\r\n')
        if not line:continue
        count+=1
        columns=line.split('\t')
        if len(columns)!=2:
            raise RuntimeError('Unexpected vdb-dump row format at spot '+str(count))
        try:lengths=[int(x.strip()) for x in columns[0].split(',')]
        except ValueError:raise RuntimeError('Invalid READ_LEN at spot '+str(count))
        types=[x.strip() for x in columns[1].split(',')]
        biological=[i for i,t in enumerate(types) if t=='SRA_READ_TYPE_BIOLOGICAL']
        technical=[i for i,t in enumerate(types) if t=='SRA_READ_TYPE_TECHNICAL']
        if (len(lengths)!=2 or len(types)!=2 or len(biological)!=1 or len(technical)!=1
            or lengths[biological[0]]<=0 or lengths[technical[0]]!=0):
            raise RuntimeError('Different read structure at spot '+str(count)+': '+line[:200])
        patterns[(lengths[0],lengths[1])]+=1
        biological_positions[biological[0]+1]+=1
        if count%1000000==0:d.log('LAYOUT_REVIEW_PROGRESS',RUN,'spots='+str(count))
    if count!=expected_records or count<=0:
        raise RuntimeError('SRA spot count does not equal validated FASTQ records: '+str(count)+' / '+str(expected_records))
    return {'all_spots_checked':count,'single_biological_spots':count,
            'zero_length_technical_spots':count,'invalid_spots':0,
            'biological_read_positions':{str(k):v for k,v in sorted(biological_positions.items())},
            'read_length_patterns':[{'lengths':list(k),'spots':v} for k,v in sorted(patterns.items())]}


def inspect_sra(executable,payload,records,log_path):
    argv=[str(executable),str(payload),'-T','SEQUENCE','-C','READ_LEN,READ_TYPE','-f','tab']
    d.log('LAYOUT_REVIEW_ALL_SPOTS_START',RUN,'expected='+str(records))
    with Path(log_path).open('w') as err:
        proc=subprocess.Popen(argv,stdout=subprocess.PIPE,stderr=err,text=True,
                              encoding='utf-8',errors='strict')
        try:
            result=scan_rows(proc.stdout,records)
            rc=proc.wait()
            if rc:raise RuntimeError('vdb-dump failed rc='+str(rc)+'; see '+str(log_path))
        finally:
            if proc.poll() is None:
                proc.terminate()
                try:proc.wait(timeout=10)
                except subprocess.TimeoutExpired:proc.kill();proc.wait()
            proc.stdout.close()
    result['command']=argv
    return result


def repair(repo):
    repo=Path(repo).resolve()
    d.ROOT=repo
    site=d.read(repo/'config/site.json')
    plans=d.read(repo/'generated/run_plan.json')
    matches=[p for p in plans if p['run_id']==RUN]
    if len(matches)!=1:raise RuntimeError('Expected exactly one target Run in current plan')
    plan=matches[0]
    if (plan['library_id']!=LIBRARY or plan['kit_id']!='SPI-ISR'
        or plan['layout']!='PAIRED' or plan['download']['mode']!='SRA_PREFETCH'):
        raise RuntimeError('Target plan identity differs from reviewed case')
    root=d.guarded(plan['paths']['run_root'],site['data_root'])
    report=root/'stats/read_layout_review.json'
    with (root/'.download.lock').open('a') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise RuntimeError('Target Run is busy; retry after it finishes')
        state=d.read(plan['paths']['state'])
        if state.get('calculation',{}).get('status') in ['RUNNING','SUBMITTED','QUEUED']:
            raise RuntimeError('Target has an active calculation; review held')
        records=d.reusable_fastq(state)
        source=state.get('download',{})
        if (not records or len(records)!=1 or records[0].get('role')!='single'
            or source.get('status')!='VALIDATED_SRA' or not source.get('vdb_validate_pass')):
            raise RuntimeError('Need the already validated SRA and exactly one stable single FASTQ')
        record=records[0]
        payload=d.sra_payload(plan['paths']['sra_accession_dir'],RUN)
        d.guarded(payload,root/'sra')
        fastq=d.guarded(record['path'],root/'fastq')
        if (str(payload)!=source['path'] or d.fp(payload)!=source['fingerprint']
            or not all(record.get(k) for k in ['gzip_pass','format_pass','stable','gzip_sha256'])
            or record['records']<=0):raise RuntimeError('Validation evidence is incomplete or stale')
        before_source=d.fp(payload);before_fastq=d.fp(fastq)
        # Bound the all-spots review to the source and FASTQ it was made for.
        if report.exists():
            effective=d.reviewed_input_plan(plan,state,site)
            review=d.read(report)
            d.log('LAYOUT_REVIEW_REUSED',RUN)
        else:
            runtime=d.read(repo/'runtime/runtime_bin.json')
            executable=Path(runtime['bin_dir'])/'vdb-dump'
            if not executable.is_file():raise RuntimeError('vdb-dump missing in installed toolkit')
            root.joinpath('stats').mkdir(exist_ok=True)
            details=inspect_sra(executable,payload,record['records'],root/'logs/read_layout_review_vdb.err')
            review={'schema_version':'1.0','status':'VERIFIED_SINGLE_BIOLOGICAL_READ',
                    'run_id':RUN,'library_id':LIBRARY,'declared_layout':'PAIRED','effective_layout':'SINGLE',
                    'reviewed_at':d.now(),'source':{'path':str(payload),'fingerprint':before_source},
                    'fastq':{'path':record['path'],'fingerprint':before_fastq,
                             'gzip_sha256':record['gzip_sha256'],'records':record['records']},
                    'reason':'Every SRA spot has a zero-length technical read and one nonempty biological read.',
                    'scope':'Input read layout only; biological eligibility and calculation routing unchanged.',
                    'official_sample_url':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6283344',
                    **details}
            effective=dict(plan,layout='SINGLE',declared_layout='PAIRED',read_layout_review=str(report))
        d.log('FASTQ_SHA256_RECHECK_START',RUN)
        if d.sha(fastq)!=record['gzip_sha256']:raise RuntimeError('FASTQ SHA256 differs from passed validation')
        if d.fp(payload)!=before_source or d.fp(fastq)!=before_fastq:
            raise RuntimeError('Input changed during layout review')
        if not report.exists():d.save(report,review)
        # Validate persisted evidence with the same gate used on future resumes.
        effective=d.reviewed_input_plan(plan,state,site)
        backup=root/'stats'/('state.before_layout_review_'+d.stamp()+'.json')
        d.save(backup,copy.deepcopy(state))
        prior=state.pop('error',None)
        if prior:state.setdefault('resolved_errors',[]).append({'resolved_at':d.now(),'error':prior,'review':str(report)})
        d.prepare_run(effective,d.input_roles(records,'SINGLE'),state,site)
        libraries=d.read(repo/'generated/library_plan.json')
        libraries=[lib for lib in libraries if lib['library_id']==LIBRARY and lib['run_ids']==[RUN]]
        if len(libraries)!=1:raise RuntimeError('Expected target single-Run library; Run ready retained')
        result=d.merge_library(libraries[0],{RUN:plan},None,site)
        if result['status']!='LIBRARY_INPUT_READY':raise RuntimeError('Library preparation incomplete: '+str(result))
        # The running dispatcher owns generated/run_results.json. Do not rewrite it.
        receipt={'finished_at':d.now(),'run_id':RUN,'library_id':LIBRARY,
                 'status':'RUN_INPUT_READY','effective_layout':'SINGLE',
                 'declared_layout':'PAIRED','records':record['records'],
                 'review':str(report),'library_result':result,
                 'run_results_note':'Current dispatcher result remains historical; state.json and ready manifests are repaired.',
                 'calculation_submission_allowed':False}
        d.save(root/'stats/read_layout_repair_result.json',receipt)
        d.log('LAYOUT_REPAIR_FINISHED',RUN,'records='+str(record['records']),result['status'])
        d.log('REVIEW_REPORT',report)
        return receipt


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo-root',type=Path,default=Path(__file__).resolve().parents[1])
    args=ap.parse_args()
    repair(args.repo_root)


if __name__=='__main__':main()
