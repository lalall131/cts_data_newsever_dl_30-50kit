#!/usr/bin/env python3
"""Install official Ubuntu x86_64 SRA Toolkit into this package, without sudo."""
import hashlib
import json
import os
import platform
import subprocess
import tarfile
import time
import urllib.request
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
URL='https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz'

def main():
    if platform.system()!='Linux' or platform.machine().lower() not in ['x86_64','amd64']:
        raise SystemExit('Automatic installer supports Ubuntu Linux x86_64; install the matching SRA Toolkit for this machine')
    runtime=ROOT/'runtime';runtime.mkdir(exist_ok=True)
    archive=runtime/'sratoolkit.official.tar.gz'
    partial=runtime/('sratoolkit.download_'+str(time.time_ns())+'.part')
    print('DOWNLOAD_TOOLKIT',URL,flush=True)
    digest=hashlib.sha256()
    with urllib.request.urlopen(URL,timeout=120) as response,partial.open('wb') as out:
        for block in iter(lambda:response.read(8*1024*1024),b''):
            out.write(block);digest.update(block)
    install=runtime/('install_'+str(time.time_ns()));install.mkdir()
    with tarfile.open(partial,'r:gz') as t:
        members=t.getmembers()
        for member in members:
            p=(install/member.name).resolve()
            if install.resolve() not in p.parents and p!=install.resolve():raise RuntimeError('Unsafe toolkit archive member')
            if member.issym() or member.islnk():
                link=(p.parent/member.linkname).resolve() if member.issym() else (install/member.linkname).resolve()
                if install.resolve() not in link.parents and link!=install.resolve():raise RuntimeError('Unsafe toolkit link')
            if not (member.isfile() or member.isdir() or member.issym() or member.islnk()):raise RuntimeError('Unsupported toolkit archive member')
        if hasattr(tarfile,'data_filter'):t.extractall(install,filter='data')
        else:t.extractall(install,members=members)
    bins=[p.parent for p in install.glob('*/bin/prefetch')]
    if len(bins)!=1:raise RuntimeError('No unique prefetch binary in official archive')
    bin_dir=bins[0];versions={}
    for name in ['prefetch','vdb-validate','fasterq-dump']:
        exe=bin_dir/name
        if not exe.is_file():raise RuntimeError('Missing toolkit binary '+name)
        r=subprocess.run([str(exe),'--version'],capture_output=True,text=True)
        if r.returncode:raise RuntimeError(name+' failed to execute: '+r.stderr)
        versions[name]=r.stdout.strip()
    os.replace(partial,archive)
    receipt={'bin_dir':str(bin_dir.resolve()),'source_url':URL,'archive_sha256':digest.hexdigest(),'versions':versions}
    (runtime/'runtime_bin.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('TOOLKIT_READY',json.dumps(receipt,ensure_ascii=False),flush=True)

if __name__=='__main__':main()
