#!/usr/bin/env python3
"""Execute prefetch -> validation -> FASTQ -> approved library input preparation."""
import argparse
import concurrent.futures
import fcntl
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
import urllib.parse
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import prepare_layout as layout

ROOT=Path(__file__).resolve().parents[1]
PRINT_LOCK=threading.Lock()
COMMAND_CONTEXT=threading.local()
GLOBAL_LOCK_FD=None

def now():return datetime.now(timezone.utc).isoformat()
def stamp():return str(time.time_ns())
def log(*args):
    with PRINT_LOCK:print(now(),*args,flush=True)
def read(path,default=None):
    return json.loads(Path(path).read_text(encoding='utf-8')) if Path(path).exists() else default
def save(path,value):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    q=p.with_name(p.name+'.writing_'+stamp())
    q.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    os.replace(q,p)
def fp(path):
    p=Path(path);s=p.stat()
    return {'bytes':s.st_size,'mtime_ns':s.st_mtime_ns,'inode':s.st_ino}
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()
def guarded(path,root):
    p=Path(path).resolve();r=Path(root).resolve()
    if p!=r and r not in p.parents:raise RuntimeError('Path outside data root: '+str(p))
    return p
def tool(name):
    p=shutil.which(name)
    if not p:
        raise RuntimeError('Missing curl; install curl on Linux before downloading ENA FASTQ' if name=='curl' else 'Missing tool '+name+'; run: bash start.sh --install-tools --check')
    return p
def command(argv,path,run_id):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    log('COMMAND',run_id,Path(argv[0]).name)
    with p.open('a',encoding='utf-8') as out:
        out.write('\n'+now()+' '+json.dumps(argv,ensure_ascii=False)+'\n');out.flush()
        inherited=tuple(fd for fd in [getattr(COMMAND_CONTEXT,'run_lock_fd',None),GLOBAL_LOCK_FD] if fd is not None)
        rc=subprocess.run(argv,stdout=out,stderr=subprocess.STDOUT,pass_fds=inherited).returncode
    if rc:raise RuntimeError(Path(argv[0]).name+' exit='+str(rc)+' log='+str(p))
def available(path,site):
    free=shutil.disk_usage(path).free
    if free<float(site.get('min_free_gb',10))*1e9:raise RuntimeError('INSUFFICIENT_FREE_SPACE bytes='+str(free))
    return free
def withdraw_ready(root,reason):
    for name in ['ready.json','library_ready.json']:
        p=Path(root)/name
        if p.exists():
            q=p.with_name(p.stem+'.invalidated_'+stamp()+'.json');os.replace(p,q)
            log('READY_WITHDRAWN',root,reason)
def sra_payload(cache,run):
    # prefetch may append the accession to -O, even when -O already ends in it.
    # Accept only the two known layouts and this exact accession, not an
    # unrestricted recursive search which could select quarantined data.
    cache=Path(cache)
    files=[]
    for directory in [cache,cache/run]:
        for suffix in ['.sra','.sralite']:
            p=directory/(run+suffix)
            if p.is_file() and p.stat().st_size:
                guarded(p,cache)
                files.append(p)
    if not files:raise FileNotFoundError('No SRA payload in '+str(cache))
    if len(files)!=1:raise RuntimeError('Ambiguous SRA payloads: '+str([str(p) for p in files]))
    return files[0]

def scan_fastq(path,run_id):
    p=Path(path);before=fp(p);count=0;ids=hashlib.sha256();content=hashlib.sha256()
    with gzip.open(p,'rb') as f:
        while True:
            header=f.readline()
            if not header:break
            seq=f.readline();plus=f.readline();qual=f.readline()
            if not seq or not plus or not qual or not header.startswith(b'@') or not plus.startswith(b'+'):
                raise RuntimeError('Invalid/truncated FASTQ '+str(p)+' record='+str(count+1))
            seq=seq.rstrip(b'\r\n');qual=qual.rstrip(b'\r\n')
            if not seq or len(seq)!=len(qual):raise RuntimeError('Sequence/quality length mismatch '+str(p))
            name=header.split()[0][1:]
            if name.endswith((b'/1',b'/2')):name=name[:-2]
            ids.update(name+b'\n');content.update(seq+b'\t'+qual+b'\n');count+=1
            if count%10000000==0:log('FASTQ_PROGRESS',run_id,p.name,count)
    if not count:raise RuntimeError('FASTQ has zero records: '+str(p))
    if before!=fp(p):raise RuntimeError('FASTQ changed during validation: '+str(p))
    digest=sha(p)
    if before!=fp(p):raise RuntimeError('FASTQ changed during fingerprint: '+str(p))
    return {'path':str(p),'records':count,'read_id_sha256':ids.hexdigest(),
            'sequence_quality_sha256':content.hexdigest(),'gzip_sha256':digest,
            'fingerprint':before,'gzip_pass':True,'format_pass':True,'stable':True}

def reusable_fastq(state):
    validation=state.get('validation',{})
    files=validation.get('input_files',[]) if isinstance(validation,dict) else []
    if not files or validation.get('status')!='FASTQ_VALIDATED':return None
    try:
        if all(Path(r['path']).is_file() and fp(r['path'])==r['fingerprint'] for r in files):return files
    except (OSError,KeyError):pass
    return None

def compress(raw,dest,site):
    pigz=shutil.which('pigz')
    if pigz:
        with Path(dest).open('wb') as f:
            inherited=tuple(fd for fd in [getattr(COMMAND_CONTEXT,'run_lock_fd',None),GLOBAL_LOCK_FD] if fd is not None)
            rc=subprocess.run([pigz,'-p',str(site.get('conversion_threads',2)),'-1','-c',str(raw)],stdout=f,pass_fds=inherited).returncode
        if rc:raise RuntimeError('pigz failed: '+str(raw))
    else:
        with Path(raw).open('rb') as src,Path(dest).open('wb') as out:
            with gzip.GzipFile(filename='',mode='wb',fileobj=out,compresslevel=int(site.get('compression_level',1)),mtime=0) as gz:
                shutil.copyfileobj(src,gz,8*1024*1024)

def input_roles(records,expected_layout):
    by_role={r['role']:r for r in records}
    if len(by_role)!=len(records):raise RuntimeError('NEEDS_READ_ROLE_REVIEW duplicate read roles')
    if expected_layout=='PAIRED' and set(by_role)=={'r1','r2'}:
        a,b=by_role['r1'],by_role['r2']
        if a['records']!=b['records'] or a['read_id_sha256']!=b['read_id_sha256']:
            raise RuntimeError('Paired read counts/IDs differ')
        return [a,b]
    if expected_layout=='SINGLE' and len(records)==1:
        records[0]['role']='single';return records
    raise RuntimeError('NEEDS_READ_ROLE_REVIEW layout='+str(expected_layout)+' files='+str([r['role'] for r in records]))

def reviewed_input_plan(plan,state,site):
    """Apply only a complete, file-bound SRA single-biological-read review."""
    path=Path(plan['paths']['stats'])/'read_layout_review.json'
    if not path.exists():return plan
    review=read(path)
    records=reusable_fastq(state)
    source=state.get('download',{})
    if (review.get('status')!='VERIFIED_SINGLE_BIOLOGICAL_READ'
        or review.get('run_id')!=plan['run_id'] or plan['layout']!='PAIRED'
        or review.get('declared_layout')!='PAIRED' or review.get('effective_layout')!='SINGLE'
        or not records or len(records)!=1 or records[0].get('role')!='single'
        or source.get('status')!='VALIDATED_SRA' or not source.get('vdb_validate_pass')):
        raise RuntimeError('READ_LAYOUT_REVIEW_INVALID; retain layout hold')
    record=records[0]
    expected={'path':record['path'],'fingerprint':record['fingerprint'],
              'gzip_sha256':record.get('gzip_sha256'),'records':record['records']}
    source_expected={'path':source['path'],'fingerprint':source['fingerprint']}
    guarded(source['path'],Path(plan['paths']['run_root'])/'sra')
    guarded(record['path'],Path(plan['paths']['run_root'])/'fastq')
    if (review.get('fastq')!=expected or review.get('source')!=source_expected
        or source['fingerprint']!=fp(source['path'])
        or review.get('all_spots_checked')!=record['records']
        or review.get('single_biological_spots')!=record['records']
        or review.get('zero_length_technical_spots')!=record['records']
        or review.get('invalid_spots')!=0 or not record.get('gzip_pass')
        or not record.get('format_pass') or not record.get('stable')
        or not record.get('gzip_sha256') or record['records']<=0):
        raise RuntimeError('READ_LAYOUT_REVIEW_STALE_OR_INCOMPLETE; retain layout hold')
    effective=dict(plan)
    effective.update(layout='SINGLE',declared_layout='PAIRED',read_layout_review=str(path))
    return effective

def prepare_run(plan,records,state,site):
    root=Path(plan['paths']['run_root']);merged=Path(plan['paths']['merged'])
    files=[]
    for r in records:
        src=guarded(r['path'],site['data_root'])
        if fp(src)!=r['fingerprint']:raise RuntimeError('Validated FASTQ changed before input publication')
        name=plan['run_id']+('_1' if r['role']=='r1' else '_2' if r['role']=='r2' else '')+'.fastq.gz'
        dest=merged/name
        if dest.exists() or dest.is_symlink():
            if dest.is_symlink() and dest.resolve()==src:pass
            else:os.replace(dest,dest.with_name(dest.name+'.saved_'+stamp()))
        if not dest.exists():
            tmp=dest.with_name(dest.name+'.link_'+stamp())
            os.symlink(os.path.relpath(src,dest.parent),tmp);os.replace(tmp,dest)
        rr=dict(r);rr.update(path=str(dest),source_path=str(src),fingerprint=fp(dest));files.append(rr)
    ready={'schema_version':'1.0','status':'INPUT_READY','ready_scope':'run',
           'run_id':plan['run_id'],'library_id':plan['library_id'],'kit_id':plan['kit_id'],
           'technology':plan['technology'],'tax_id':plan['species'].get('tax_id'),
           'layout':plan['layout'],'source_mode':'SRA_CONVERTED_FASTQ',
           'prepared_at':now(),'input_files':files,'run_ids':[plan['run_id']],
           'comparison_to_ENA_gzip_MD5':'NOT_APPLICABLE',
           'calculation_submission_allowed':False}
    if plan.get('read_layout_review'):
        ready.update(declared_layout=plan['declared_layout'],read_layout_review=plan['read_layout_review'])
    save(root/'ready.json',ready)
    state['input_preparation']={'status':'RUN_INPUT_READY','merge_manifest':str(root/'ready.json'),'tmp_bytes':None}
    save(plan['paths']['state'],state)
    log('RUN_INPUT_READY',plan['run_id'])

def md5_file(path):
    h=hashlib.md5()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

def process_ena_archive(plan,site,state):
    """Preserve exact public FASTQ bytes/headers; no barcode trim or FLNC claim."""
    root=Path(plan['paths']['run_root']);run=plan['run_id'];files=plan['download']['ena_files']
    old=state.get('validation',{});stored=old.get('input_files',[]) if isinstance(old,dict) else []
    reusable=False
    try:
        reusable=old.get('status')=='ARCHIVE_FASTQ_VALIDATED' and len(stored)==len(files) and all(
            record['source_expected_md5']==item['expected_md5'] and record['path']==item['destination'] and
            fp(item['destination'])==record['fingerprint'] and record['fingerprint']['bytes']==item['expected_bytes']
            for record,item in zip(stored,files))
    except (OSError,KeyError):pass
    if reusable:
        log('ARCHIVE_FASTQ_REUSED',run)
        return {'run_id':run,'status':'ARCHIVE_FASTQ_VALIDATED_REUSED_KIT_PREPARATION_REQUIRED'}
    withdraw_ready(root,'special Kit needs raw archive validation and dedicated input adapter')
    records=[]
    for item in files:
        dest=guarded(item['destination'],site['data_root']);part=Path(str(dest)+'.part')
        parsed=urllib.parse.urlsplit(item['url'])
        if parsed.scheme!='https' or parsed.hostname not in ['ftp.sra.ebi.ac.uk','ftp.ebi.ac.uk']:
            raise RuntimeError('Unapproved ENA source URL '+item['url'])
        expected=int(item['expected_bytes']);digest=item['expected_md5'].lower()
        if not expected or not re.fullmatch('[a-f0-9]{32}',digest):raise RuntimeError('Invalid ENA source fingerprint')
        last=None;record=None
        for attempt in range(1,int(site.get('download_attempts',3))+1):
            available(root,site)
            state['download']={'status':'DOWNLOADING_ARCHIVE_FASTQ','source_mode':'ENA_ARCHIVE_FASTQ',
                               'file':dest.name,'attempt':attempt,'started_at':now()};save(plan['paths']['state'],state)
            try:
                if dest.exists():
                    before=fp(dest);good=before['bytes']==expected and md5_file(dest)==digest and fp(dest)==before
                    if not good:os.replace(dest,dest.with_name(dest.name+'.unverified_saved_'+stamp()))
                if not dest.exists():
                    if part.exists() and part.stat().st_size>expected:
                        os.replace(part,part.with_name(part.name+'.oversize_saved_'+stamp()))
                    if not part.exists() or part.stat().st_size!=expected:
                        command([tool('curl'),'-q','-L','--fail','--retry','0','--connect-timeout','30',
                                 '--speed-time','180','--speed-limit','1024','--continue-at','-',
                                 '--output',str(part),item['url']],root/'logs/ena_download.log',run)
                    before=fp(part)
                    if before['bytes']!=expected:raise RuntimeError('ENA size mismatch '+dest.name)
                    actual=md5_file(part)
                    if fp(part)!=before:raise RuntimeError('ENA partial changed during MD5 check')
                    if actual!=digest:
                        os.replace(part,part.with_name(part.name+'.md5_mismatch_saved_'+stamp()))
                        raise RuntimeError('ENA MD5 mismatch '+dest.name+' expected='+digest+' actual='+actual)
                    record=scan_fastq(part,run)
                    os.replace(part,dest);record.update(path=str(dest),fingerprint=fp(dest))
                else:record=scan_fastq(dest,run)
                record.update(source_expected_md5=digest,source_MD5_verified=True,source_url=item['url'],
                    role='r1' if dest.name==run+'_1.fastq.gz' else 'r2' if dest.name==run+'_2.fastq.gz' else 'archive_single')
                records.append(record);log('ENA_ARCHIVE_VERIFIED',run,dest.name,'bytes='+str(expected));break
            except Exception as error:
                last=error;log('ENA_ARCHIVE_RETRY',run,dest.name,attempt,str(error))
                if attempt==int(site.get('download_attempts',3)):raise
                time.sleep(min(attempt*5,15))
        if record is None:raise RuntimeError(str(last))
    # An archive layout is recorded, not substituted for cell read-role/FLNC interpretation.
    if plan['layout']=='PAIRED' and len(records)==2 and {r['role'] for r in records}=={'r1','r2'}:
        input_roles(records,'PAIRED')
    state['download']={'status':'ENA_ARCHIVE_FASTQ_DOWNLOADED','source_mode':'ENA_ARCHIVE_FASTQ',
                       'finished_at':now(),'bytes':sum(r['fingerprint']['bytes'] for r in records)}
    state['validation']={'status':'ARCHIVE_FASTQ_VALIDATED','official_ENA_MD5':True,
                         'gzip':True,'FASTQ_format':True,'input_files':records,'checked_at':now()}
    state['input_preparation']={'status':'KIT_SPECIFIC_PREPARATION_REQUIRED',
                               'special_input_profile':plan.get('special_input_profile'),
                               'reason':plan.get('special_calculation_hold'),'barcode_or_FLNC_certified':False}
    save(plan['paths']['state'],state)
    save(root/'raw_fastq_manifest.json',{'schema_version':'1.0','run_id':run,'library_id':plan['library_id'],
        'status':'RAW_ARCHIVE_FASTQ_VALIDATED','input_files':records,'kit_id':plan['kit_id'],
        'species':plan['species'],'special_input_profile':plan.get('special_input_profile'),
        'dedicated_preparation_required':True,'calculation_submission_allowed':False})
    return {'run_id':run,'status':'ARCHIVE_FASTQ_VALIDATED_KIT_PREPARATION_REQUIRED'}

def process_run(plan,site,stage,conversion_slots):
    run=plan['run_id'];root=Path(plan['paths']['run_root']);state_path=Path(plan['paths']['state'])
    guarded(root,site['data_root'])
    with (root/'.download.lock').open('a') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:return {'run_id':run,'status':'BUSY'}
        COMMAND_CONTEXT.run_lock_fd=lock.fileno()
        state=read(state_path,layout.initial_run_state(run))
        for section in ['download','validation','input_preparation','calculation']:
            if not isinstance(state.get(section),dict):state[section]={'status':state.get(section) or 'UNKNOWN'}
        try:
            current=state.get('calculation',{})
            if isinstance(current,dict) and current.get('status') in ['RUNNING','SUBMITTED','QUEUED']:
                return {'run_id':run,'status':'ACTIVE_CALCULATION_HELD'}
            if plan['download']['mode']=='GSA_MANUAL_CHANNEL':
                state['download']={'status':'MANUAL_SOURCE_REQUIRED','error':'CRR requires a GSA channel'}
                save(state_path,state);return {'run_id':run,'status':'MANUAL_SOURCE_REQUIRED'}
            available(root,site)
            if plan['download']['mode']=='ENA_ARCHIVE_FASTQ':return process_ena_archive(plan,site,state)
            cache=Path(plan['paths']['sra_accession_dir']);cache.parent.mkdir(parents=True,exist_ok=True)
            old=state.get('download',{})
            payload=None
            if isinstance(old,dict) and old.get('status')=='VALIDATED_SRA':
                try:
                    candidate=sra_payload(cache,run)
                    if old.get('fingerprint')==fp(candidate):payload=candidate
                except OSError:pass
                except RuntimeError:pass
            if payload is None:
                last=None
                for attempt in range(1,int(site.get('download_attempts',3))+1):
                    state['download']={'status':'DOWNLOADING','attempt':attempt,'source_mode':'prefetch','started_at':now()}
                    save(state_path,state)
                    try:
                        try:
                            payload=sra_payload(cache,run)
                            log('SRA_CACHE_FOUND',run,payload)
                        except FileNotFoundError:
                            command([tool('prefetch'),run,'--max-size',str(site.get('prefetch_max_size','100G')),'-O',str(cache)],root/'logs/prefetch.log',run)
                        payload=sra_payload(cache,run);before=fp(payload)
                        try:command([tool('vdb-validate'),str(payload)],root/'logs/vdb_validate.log',run)
                        except RuntimeError:
                            guarded(cache,root/'sra');os.rename(cache,cache.with_name(cache.name+'.validation_failed_'+stamp()));payload=None;raise
                        if fp(payload)!=before:raise RuntimeError('SRA changed during validation')
                        state['download']={'status':'VALIDATED_SRA','source_mode':'prefetch','finished_at':now(),
                                           'path':str(payload),'fingerprint':fp(payload),'bytes':payload.stat().st_size,
                                           'vdb_validate_pass':True,'attempt':attempt}
                        save(state_path,state);break
                    except Exception as error:
                        last=error;log('DOWNLOAD_RETRY',run,attempt,str(error))
                        if attempt==int(site.get('download_attempts',3)):raise
                        time.sleep(min(5*attempt,15))
                if payload is None:raise RuntimeError(str(last))
            log('SRA_VALIDATED',run)
            if state.pop('error',None) is not None:save(state_path,state)
            if stage=='download':return {'run_id':run,'status':'VALIDATED_SRA'}
            if plan['kit_id'] in ['SPI-PBLR','scFLUENT-seq SPI','scGRO-seq SR']:
                state['input_preparation']={'status':'KIT_SPECIFIC_PREPARATION_REQUIRED','error':'FLNC/barcode/UMI preparation needs the dedicated adapter'}
                save(state_path,state);return {'run_id':run,'status':'KIT_SPECIFIC_PREPARATION_REQUIRED'}
            plan=reviewed_input_plan(plan,state,site)
            cached=reusable_fastq(state)
            if cached:
                records=input_roles(cached,plan['layout']);prepare_run(plan,records,state,site)
                return {'run_id':run,'status':'RUN_INPUT_READY_REUSED'}
            with conversion_slots:
                available(root,site);withdraw_ready(root,'preparing new validated FASTQ')
                work=root/'tmp'/('conversion_'+stamp());raw=work/'raw';scratch=work/'scratch';candidate=work/'gzip'
                for p in [raw,scratch,candidate]:p.mkdir(parents=True,exist_ok=True)
                state['validation']={'status':'CONVERTING','vdb_validate':True,'input_files':[]};save(state_path,state)
                before=fp(payload)
                command([tool('fasterq-dump'),str(payload.parent),'--threads',str(site.get('conversion_threads',2)),
                         '--seq-defline','@$ac.$si/$ri','--qual-defline','+',
                         '--outdir',str(raw),'--temp',str(scratch)],root/'logs/fasterq_dump.log',run)
                if fp(payload)!=before:raise RuntimeError('SRA changed during conversion')
                raws=sorted(p for p in raw.glob(run+'*.fastq') if p.stat().st_size)
                if not raws:raise RuntimeError('fasterq-dump produced no nonempty FASTQ')
                records=[]
                for p in raws:
                    dest=candidate/(p.name+'.gz');compress(p,dest,site)
                    r=scan_fastq(dest,run)
                    r['role']='r1' if p.name==run+'_1.fastq' else 'r2' if p.name==run+'_2.fastq' else 'single'
                    records.append(r)
                # Preserve every biological output; ambiguous roles remain a review case.
                for r in records:
                    dest=root/'fastq'/Path(r['path']).name
                    if dest.exists():os.replace(dest,dest.with_name(dest.name+'.saved_'+stamp()))
                    os.replace(r['path'],dest);r.update(path=str(dest),fingerprint=fp(dest))
                state['validation']={'status':'FASTQ_VALIDATED','vdb_validate':True,
                                     'gzip':True,'FASTQ_format':True,'input_files':records,
                                     'official_MD5_for_ENA_file':'NOT_APPLICABLE_TO_SRA_CONVERSION','checked_at':now()}
                save(state_path,state)
                # Only generated uncompressed temporary FASTQs are removed after their gzip copies passed full validation.
                for p in raws:guarded(p,work);p.unlink()
                records=input_roles(records,plan['layout']);prepare_run(plan,records,state,site)
            return {'run_id':run,'status':'RUN_INPUT_READY'}
        except Exception as error:
            withdraw_ready(root,str(error));state['error']={'at':now(),'message':str(error)}
            state['input_preparation']={'status':'FAILED_OR_NEEDS_REVIEW','error':str(error)}
            if state.get('download',{}).get('status')!='VALIDATED_SRA':state['download']['status']='FAILED'
            save(state_path,state);log('RUN_ERROR',run,str(error))
            return {'run_id':run,'status':'FAILED_OR_NEEDS_REVIEW','error':str(error)}

def checked_ready(plan,site):
    value=read(plan['paths']['ready'])
    if not value or value.get('status')!='INPUT_READY':raise RuntimeError('Run input not ready: '+plan['run_id'])
    if value.get('run_id')!=plan['run_id'] or value.get('library_id')!=plan['library_id']:raise RuntimeError('Run ready identity mismatch')
    if value.get('schema_version')!='1.0' or not value.get('input_files'):raise RuntimeError('Incomplete ready manifest')
    input_roles(value['input_files'],value['layout'])
    for record in value['input_files']:
        if not record.get('gzip_pass') or not record.get('format_pass') or not record.get('stable') or record.get('records',0)<=0:
            raise RuntimeError('Ready manifest lacks passed validation')
        p=guarded(record['path'],site['data_root'])
        if fp(p)!=record['fingerprint']:raise RuntimeError('Run input changed: '+str(p))
    return value

def merge_library(library,plans,decision,site):
    ident=library['library_id'];root=Path(library['library_root']);root.mkdir(parents=True,exist_ok=True)
    if any(plans[run]['download']['mode']=='ENA_ARCHIVE_FASTQ' and plans[run].get('special_calculation_hold') for run in library['run_ids']):
        withdraw_ready(root,'dedicated Kit input adapter is unconfirmed')
        return {'library_id':ident,'status':'KIT_SPECIFIC_PREPARATION_REQUIRED'}
    if len(library['run_ids'])==1:
        try:
            ready=checked_ready(plans[library['run_ids'][0]],site)
            value=dict(ready);value.update(ready_scope='library',library_run_ids=library['run_ids'],merge_action='SINGLE_RUN')
            save(root/'library_ready.json',value)
            return {'library_id':ident,'status':'LIBRARY_INPUT_READY'}
        except Exception as error:
            withdraw_ready(root,str(error))
            return {'library_id':ident,'status':'WAITING_FOR_RUN_INPUT','error':str(error)}
    if not decision or not decision.get('approved') or decision.get('action') not in ['CONCAT','USE_RUNS']:
        return {'library_id':ident,'status':'MERGE_REVIEW_REQUIRED'}
    selected=decision.get('included_runs',[])
    if not selected or len(set(selected))!=len(selected) or not set(selected).issubset(library['run_ids']):
        return {'library_id':ident,'status':'MERGE_DECISION_INVALID'}
    if decision['action']=='CONCAT' and set(selected)!=set(library['run_ids']):
        return {'library_id':ident,'status':'CONCAT_REQUIRES_COMPLETE_RUN_SET'}
    if not decision.get('reason'):
        return {'library_id':ident,'status':'MERGE_REASON_REQUIRED'}
    with (root/'.merge.lock').open('a') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:return {'library_id':ident,'status':'BUSY'}
        try:
            inputs=[checked_ready(plans[run],site) for run in selected]
            layouts={r['layout'] for r in inputs};taxids={r['tax_id'] for r in inputs}
            if len(layouts)!=1 or len(taxids)!=1:raise RuntimeError('Incompatible layout/species across Runs')
            ordered=[{r['role']:r for r in inp['input_files']} for inp in inputs]
            roles=set(ordered[0])
            if any(set(v)!=roles for v in ordered):raise RuntimeError('Incompatible read roles across Runs')
            signatures=[tuple((role,v[role]['records'],v[role]['sequence_quality_sha256']) for role in sorted(roles)) for v in ordered]
            if len(signatures)!=len(set(signatures)):raise RuntimeError('DUPLICATE_RUN_SEQUENCE_CONTENT; review before concatenation')
            prior=read(root/'library_ready.json')
            source_keys=[{'run':run,'files':inputs[i]['input_files']} for i,run in enumerate(selected)]
            if prior and prior.get('source_inputs')==source_keys:
                if all(Path(r['path']).is_file() and fp(r['path'])==r['fingerprint'] for r in prior['input_files']):
                    return {'library_id':ident,'status':'LIBRARY_INPUT_READY_REUSED'}
            if (root/'library_ready.json').exists():os.replace(root/'library_ready.json',root/('library_ready.invalidated_'+stamp()+'.json'))
            available(root,site);records=[]
            for role in sorted(roles):
                name=ident+('_1' if role=='r1' else '_2' if role=='r2' else '')+'.fastq.gz'
                candidate=root/'tmp'/(name+'.assembling_'+stamp())
                with candidate.open('wb') as out:
                    for v in ordered:
                        source=guarded(v[role]['path'],site['data_root'])
                        with source.open('rb') as f:shutil.copyfileobj(f,out,8*1024*1024)
                r=scan_fastq(candidate,ident);r['role']=role
                if r['records']!=sum(v[role]['records'] for v in ordered):raise RuntimeError('Merged record count differs')
                records.append(r)
            input_roles(records,next(iter(layouts)))
            for run in selected:checked_ready(plans[run],site)
            for r in records:
                role=r['role'];name=ident+('_1' if role=='r1' else '_2' if role=='r2' else '')+'.fastq.gz'
                dest=root/'merged'/name
                if dest.exists():os.replace(dest,dest.with_name(dest.name+'.saved_'+stamp()))
                os.replace(r['path'],dest);r.update(path=str(dest),fingerprint=fp(dest))
            value={'schema_version':'1.0','status':'INPUT_READY','ready_scope':'library','library_id':ident,
                   'kit_id':library['kit_id'],'technology':library['technology'],'tax_id':next(iter(taxids)),
                   'layout':next(iter(layouts)),'prepared_at':now(),'input_files':records,
                   'run_ids':selected,'library_run_ids':library['run_ids'],
                   'excluded_runs':sorted(set(library['run_ids'])-set(selected)),'merge_action':decision['action'],
                   'merge_reason':decision.get('reason'),'source_inputs':source_keys,'calculation_submission_allowed':False}
            save(root/'library_ready.json',value);log('LIBRARY_INPUT_READY',ident)
            return {'library_id':ident,'status':'LIBRARY_INPUT_READY'}
        except Exception as error:
            withdraw_ready(root,str(error))
            save(root/'stats/merge_error.json',{'at':now(),'error':str(error)})
            log('MERGE_ERROR',ident,str(error));return {'library_id':ident,'status':'FAILED_OR_NEEDS_REVIEW','error':str(error)}

def activate_runtime():
    receipt=read(ROOT/'runtime/runtime_bin.json',{})
    if receipt.get('bin_dir'):os.environ['PATH']=receipt['bin_dir']+os.pathsep+os.environ.get('PATH','')

def main():
    global GLOBAL_LOCK_FD
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stage',choices=['download','fastq','merge','all'],default='all')
    ap.add_argument('--check',action='store_true')
    ap.add_argument('--install-tools',action='store_true')
    ap.add_argument('--workers',type=int)
    ap.add_argument('--limit',type=int)
    ap.add_argument('--only-kit')
    ap.add_argument('--only-run',action='append')
    ap.add_argument('--only-added',action='store_true',help='Download only the 76 added mouse Runs')
    args=ap.parse_args();os.chdir(ROOT)
    if os.name!='posix':raise SystemExit('Run this package on the Linux server')
    activate_runtime()
    site=read(ROOT/'config/site.json')
    subprocess.run([sys.executable,str(ROOT/'tools/prepare_layout.py'),'--out-dir','generated'],check=True)
    plans=read(ROOT/'generated/run_plan.json');libraries=read(ROOT/'generated/library_plan.json')
    if args.only_added:plans=[p for p in plans if p.get('selection_batch')=='added_three_kits_20261014']
    if args.only_kit:plans=[p for p in plans if p['kit_id']==args.only_kit]
    if args.only_run:plans=[p for p in plans if p['run_id'] in set(args.only_run)]
    if args.limit is not None:
        if args.limit<1:raise SystemExit('--limit must be positive')
        plans=plans[:args.limit]
    if not plans:raise SystemExit('No Runs selected')
    ids={p['library_id'] for p in plans};libraries=[l for l in libraries if l['library_id'] in ids]
    needs_sra=any(p['download']['mode']=='SRA_PREFETCH' for p in plans) and args.stage!='merge'
    required_sra=(['prefetch','vdb-validate']+(['fasterq-dump'] if args.stage in ['fastq','all'] else [])) if needs_sra else []
    required=required_sra+(['curl'] if args.stage!='merge' and any(p['download']['mode']=='ENA_ARCHIVE_FASTQ' for p in plans) else [])
    def outdated(name):
        p=shutil.which(name)
        if not p:return True
        result=subprocess.run([p,'--version'],capture_output=True,text=True)
        match=re.search(r'\b(\d+)\.(\d+)\.\d+\b',result.stdout+' '+result.stderr)
        return result.returncode!=0 or not match or int(match.group(1))<3
    if args.install_tools and any(outdated(n) for n in required_sra):
        subprocess.run([sys.executable,str(ROOT/'tools/bootstrap_sra.py')],check=True);activate_runtime()
    paths={name:tool(name) for name in required}
    for name,p in paths.items():
        result=subprocess.run([p,'--version'],capture_output=True,text=True)
        if result.returncode:raise SystemExit(name+' cannot execute: '+result.stderr)
        match=re.search(r'\b(\d+)\.(\d+)\.\d+\b',result.stdout+' '+result.stderr)
        if name!='curl' and (not match or int(match.group(1))<3):raise SystemExit(name+' needs SRA Toolkit 3.x or newer; run bash start.sh --install-tools --check')
        log('TOOL',name,p,result.stdout.strip())
    workers=args.workers or int(site.get('download_workers',2))
    if workers<1 or workers>8:raise SystemExit('--workers must be between 1 and 8')
    site['conversion_threads']=int(site.get('conversion_threads',2))
    if site['conversion_threads']<1:raise SystemExit('conversion_threads must be positive')
    if not 1<=int(site.get('conversion_workers',1))<=8:raise SystemExit('conversion_workers must be between 1 and 8')
    if not 1<=int(site.get('download_attempts',3))<=9:raise SystemExit('download_attempts must be between 1 and 9')
    if args.check:
        parent=Path(site['data_root'])
        while not parent.exists() and parent!=parent.parent:parent=parent.parent
        print(json.dumps({'check':'PASS','selected_runs':len(plans),'selected_libraries':len(libraries),
                          'workers':workers,'conversion_workers':site.get('conversion_workers',1),
                          'free_bytes':shutil.disk_usage(parent).free,'data_root':site['data_root'],
                          'stage':args.stage,'download_started':False,'tool_paths':paths},ensure_ascii=False,indent=2));return
    lock_path=Path(site['data_root'])/'.cts_download_merge.lock'
    lock_path.parent.mkdir(parents=True,exist_ok=True)
    with lock_path.open('a') as global_lock:
        try:fcntl.flock(global_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise SystemExit('A download/merge process already holds the project lock')
        GLOBAL_LOCK_FD=global_lock.fileno()
        subprocess.run([sys.executable,str(ROOT/'tools/prepare_layout.py'),'--out-dir','generated','--materialize'],check=True)
        run_results=[];merge_results=[]
        if args.stage!='merge':
            slots=threading.Semaphore(int(site.get('conversion_workers',1)))
            with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
                futures=[executor.submit(process_run,p,site,args.stage,slots) for p in plans]
                for future in concurrent.futures.as_completed(futures):
                    result=future.result();run_results.append(result)
                    save(ROOT/'generated/run_results.json',run_results)
        if args.stage in ['all','merge']:
            all_plans={p['run_id']:p for p in read(ROOT/'generated/run_plan.json')}
            decisions=read(ROOT/'config/library_decisions.json',{})
            for library in libraries:
                result=merge_library(library,all_plans,decisions.get(library['library_id']),site)
                merge_results.append(result);save(ROOT/'generated/merge_results.json',merge_results)
        summary={'finished_at':now(),'run_statuses':dict(Counter(r['status'] for r in run_results)),
                 'library_statuses':dict(Counter(r['status'] for r in merge_results)),
                 'selected_runs':len(plans),'selected_libraries':len(libraries),'calculation_jobs_submitted':0}
        save(ROOT/'generated/execution_summary.json',summary)
        subprocess.run([sys.executable,str(ROOT/'tools/inventory_run_storage.py')],check=False)
        print(json.dumps(summary,ensure_ascii=False,indent=2))
        failed=any(r['status']=='FAILED_OR_NEEDS_REVIEW' for r in run_results+merge_results)
        if failed:raise SystemExit(2)

if __name__=='__main__':main()
