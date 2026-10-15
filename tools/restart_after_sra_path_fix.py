#!/usr/bin/env python3
"""Restart this package after the SRA path fix, letting old transfers finish.

Matches this user's exact Python entrypoint via /proc, never a broad pkill.
Existing child downloads/conversions are allowed to finish. The project lock
and recorded child process start times prevent overlapping writers.
"""
import fcntl
import hashlib
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
RUNNER=ROOT/'tools/run_download_merge.py'

def log(*parts):
    print(time.strftime('%Y-%m-%d %H:%M:%S'),*parts,flush=True)

def processes():
    result={}
    for path in Path('/proc').iterdir():
        if not path.name.isdigit():continue
        try:
            if path.stat().st_uid!=os.getuid():continue
            raw=(path/'stat').read_text()
            fields=raw[raw.rfind(')')+2:].split()
            argv=[os.fsdecode(a) for a in (path/'cmdline').read_bytes().split(b'\0') if a]
            result[int(path.name)]={'pid':int(path.name),'ppid':int(fields[1]),
                                   'state':fields[0],'start':fields[19],'argv':argv}
        except (OSError,ValueError,IndexError):continue
    return result

def runner_index(argv):
    # Only a Python command with this exact script as its script argument.
    if not argv or not Path(argv[0]).name.startswith('python'):return None
    index=1
    while index<len(argv) and argv[index] in ['-u','-B','-E','-s','-S']:
        index+=1
    if index<len(argv):
        try:
            if Path(argv[index]).resolve()==RUNNER.resolve():return index
        except OSError:pass
    return None

def active_runners(table):
    found=[]
    for item in table.values():
        index=runner_index(item['argv'])
        if index is not None and '--check' not in item['argv'][index+1:] and item['state']!='Z':
            found.append(item)
    return found

def descendants(table,parent):
    selected=set()
    while True:
        more={pid for pid,item in table.items()
              if item['ppid'] in selected|{parent} and pid!=parent}
        if more<=selected:break
        selected|=more
    return {pid:table[pid] for pid in selected}

def live_record(item,table):
    current=table.get(item['pid'])
    return current is not None and current['start']==item['start'] and current['state']!='Z'

def main():
    if sys.platform!='linux':raise SystemExit('Run this controller on the Linux server')
    os.chdir(ROOT)
    runtime=ROOT/'runtime';runtime.mkdir(exist_ok=True)
    logs=ROOT/'logs';logs.mkdir(exist_ok=True)
    # Prevent two repair controllers from racing each other.
    with (runtime/'.sra_path_fix_restart.lock').open('a') as controller:
        try:fcntl.flock(controller,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise SystemExit('Another path-fix restart controller is active')
        code=RUNNER.read_text(encoding='utf-8')
        if 'cache/run' not in code or "str(payload.parent),'--threads'" not in code:
            raise SystemExit('Updated runner not found; fetch the path-fix update first')
        source_sha=hashlib.sha256(RUNNER.read_bytes()).hexdigest()
        receipt=runtime/'sra_path_fix_restart.json'
        if receipt.exists():
            previous=json.loads(receipt.read_text())
            if previous.get('status')=='STARTED' and previous.get('runner_sha256')==source_sha:
                log('ALREADY_RESTARTED',previous)
                return
        table=processes();roots=active_runners(table)
        if len(roots)>1:raise SystemExit('Multiple package runners found; inspect their PIDs before restarting')
        child_records={}
        old_pid=None
        arguments=['--stage','all']
        if roots:
            old=roots[0];old_pid=old['pid']
            index=runner_index(old['argv']);arguments=old['argv'][index+1:]
            # Freeze the dispatcher briefly to capture all its existing children.
            # Child tools remain running; no downloaded files are moved/deleted.
            paused=False
            try:
                if not live_record(old,processes()):raise ProcessLookupError()
                os.kill(old_pid,signal.SIGSTOP);paused=True
                time.sleep(0.2)
                child_records=descendants(processes(),old_pid)
                log('OLD_DISPATCHER',old_pid,'arguments='+repr(arguments))
                log('FINISH_EXISTING_CHILDREN',sorted(child_records))
                if not live_record(old,processes()):raise ProcessLookupError()
                os.kill(old_pid,signal.SIGTERM)
            except ProcessLookupError:
                log('OLD_DISPATCHER_ALREADY_EXITED',old_pid)
            finally:
                if paused:
                    try:
                        if live_record(old,processes()):os.kill(old_pid,signal.SIGCONT)
                    except ProcessLookupError:pass
            while live_record(old,processes()):time.sleep(0.5)
        data_root=Path(json.loads((ROOT/'config/site.json').read_text())['data_root'])
        data_root.mkdir(parents=True,exist_ok=True)
        started=time.monotonic();last=-100
        # Files inherited by child tools can keep the project lock alive after
        # the Python dispatcher exits. Wait for both children and that lock.
        with (data_root/'.cts_download_merge.lock').open('a') as project:
            while True:
                table=processes()
                # Track any further descendants of the original child tools.
                for item in list(child_records.values()):
                    if live_record(item,table):child_records.update(descendants(table,item['pid']))
                remaining=[pid for pid,item in child_records.items() if live_record(item,table)]
                acquired=False
                if not remaining:
                    try:
                        fcntl.flock(project,fcntl.LOCK_EX|fcntl.LOCK_NB)
                        acquired=True
                    except BlockingIOError:pass
                if acquired:break
                elapsed=int(time.monotonic()-started)
                if elapsed-last>=30:
                    log('WAIT_TRANSFER_FINISH','seconds='+str(elapsed),'child_pids='+repr(remaining))
                    last=elapsed
                time.sleep(3)
            log('OLD_WRITERS_FINISHED','SRA/FASTQ/cache files retained')
            if active_runners(processes()):raise SystemExit('Another runner appeared; restart cancelled')
            subprocess.run([sys.executable,str(RUNNER),'--check'],check=True)
            fcntl.flock(project,fcntl.LOCK_UN)
        suffix=time.strftime('%Y%m%d_%H%M%S')+'_'+str(os.getpid())
        output=logs/('download_pathfixed_'+suffix+'.out')
        with output.open('ab') as out:
            child=subprocess.Popen([sys.executable,'-u',str(RUNNER),*arguments],
                                   stdin=subprocess.DEVNULL,stdout=out,stderr=subprocess.STDOUT,
                                   start_new_session=True,close_fds=True)
        Path(str(output)+'.pid').write_text(str(child.pid)+'\n')
        time.sleep(3)
        if child.poll() is not None:
            log('RESTART_EXITED',child.returncode,'log='+str(output))
            raise SystemExit(1)
        result={'status':'STARTED','old_pid':old_pid,'new_pid':child.pid,
                'log':str(output),'arguments':arguments,'runner_sha256':source_sha}
        temporary=receipt.with_name(receipt.name+'.writing_'+str(os.getpid()))
        temporary.write_text(json.dumps(result,indent=2)+'\n');temporary.replace(receipt)
        log('RESTARTED',json.dumps(result))

if __name__=='__main__':main()
