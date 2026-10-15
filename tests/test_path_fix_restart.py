"""Linux integration: stop only this dispatcher, drain child, resume once."""
import json
import os
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
FAKE='''import fcntl,json,os,subprocess,sys,time
from pathlib import Path
# cache/run
# str(payload.parent),'--threads'
root=Path(__file__).resolve().parents[1]
if '--check' in sys.argv:
    print('MOCK_CHECK_PASS',flush=True);raise SystemExit(0)
lock=(root/'data/.cts_download_merge.lock').open('a')
fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
if not (root/'runtime/fixed').exists():
    program="import time;from pathlib import Path;time.sleep(4);Path("+repr(str(root/'runtime/child_finished'))+").write_text('complete')"
    child=subprocess.Popen([sys.executable,'-c',program],pass_fds=(lock.fileno(),))
    (root/'runtime/child_pid').write_text(str(child.pid))
else:
    assert (root/'runtime/child_finished').exists(),'old writer was interrupted or overlapped'
    (root/'runtime/new_started').write_text(str(os.getpid()))
while True:time.sleep(1)
'''

with tempfile.TemporaryDirectory(prefix='cts-path-restart-') as temp:
    pkg=Path(temp)
    for directory in ['tools','config','data','runtime','logs']:(pkg/directory).mkdir()
    helper=pkg/'tools/restart_after_sra_path_fix.py'
    shutil.copy2(ROOT/'tools/restart_after_sra_path_fix.py',helper)
    runner=pkg/'tools/run_download_merge.py';runner.write_text(FAKE)
    config={'data_root':str(pkg/'data'),'download_workers':8,'conversion_workers':4,'conversion_threads':16}
    (pkg/'config/site.json').write_text(json.dumps(config))
    old=subprocess.Popen([sys.executable,'-u',str(runner),'--stage','all'],stdout=subprocess.DEVNULL)
    new_pid=None
    try:
        deadline=time.monotonic()+10
        while not (pkg/'runtime/child_pid').exists():
            assert old.poll() is None
            if time.monotonic()>deadline:raise AssertionError('mock child startup timed out')
            time.sleep(0.05)
        (pkg/'runtime/fixed').write_text('path fix installed')
        result=subprocess.run([sys.executable,'-u',str(helper)],cwd=pkg,capture_output=True,text=True,timeout=25)
        assert result.returncode==0,result.stdout+'\n'+result.stderr
        receipt=json.loads((pkg/'runtime/sra_path_fix_restart.json').read_text())
        new_pid=receipt['new_pid']
        assert receipt['old_pid']==old.pid and receipt['arguments']==['--stage','all']
        assert old.wait(timeout=2)==-signal.SIGTERM
        assert (pkg/'runtime/child_finished').read_text()=='complete'
        assert int((pkg/'runtime/new_started').read_text())==new_pid
        assert json.loads((pkg/'config/site.json').read_text())==config
        result=subprocess.run([sys.executable,'-u',str(helper)],cwd=pkg,capture_output=True,text=True,timeout=10)
        assert result.returncode==0 and 'ALREADY_RESTARTED' in result.stdout,result.stdout+'\n'+result.stderr
        assert json.loads((pkg/'runtime/sra_path_fix_restart.json').read_text())['new_pid']==new_pid
        print('PASS: exact dispatcher targeting, child transfer completes, project lock respected, parameters retained, detached restart, idempotent repeat')
    finally:
        if old.poll() is None:
            old.kill();old.wait()
        if new_pid:
            try:os.kill(new_pid,signal.SIGTERM)
            except ProcessLookupError:pass
