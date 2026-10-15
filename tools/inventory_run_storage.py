#!/usr/bin/env python3
"""Record Run storage snapshots, including merged/tmp; never infer integrity from size."""
import argparse
import csv
import json
import os
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

def measure(path):
    root=Path(path)
    result={'path':str(root),'exists':None,'files':None,'logical_bytes':None,
            'unique_inode_bytes':None,'symlinks':0,'errors':[],'extension_counts':{}}
    try:result['exists']=root.is_dir()
    except OSError as error:result['errors'].append(str(error));return result
    if not result['exists']:
        return result
    files=0; logical=0; inodes={}; types=Counter()
    def onerror(error):result['errors'].append(str(error))
    for parent,dirs,names in os.walk(root,followlinks=False,onerror=onerror):
        symlink_dirs=[d for d in dirs if (Path(parent)/d).is_symlink()]
        result['symlinks']+=len(symlink_dirs)
        dirs[:]=[d for d in dirs if d not in symlink_dirs]
        for name in names:
            p=Path(parent)/name
            try:
                if p.is_symlink():result['symlinks']+=1;continue
                s=p.stat()
                if not p.is_file():continue
                files+=1;logical+=s.st_size;inodes[(s.st_dev,s.st_ino)]=s.st_size
                ext='.fastq.gz' if name.endswith('.fastq.gz') else p.suffix or '<none>'
                types[ext]+=1
            except OSError as error:result['errors'].append(str(error))
    result.update(files=files,logical_bytes=logical,unique_inode_bytes=sum(inodes.values()),extension_counts=dict(types))
    result['complete']=not result['errors']
    return result

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run-plan',default='generated/run_plan.json')
    ap.add_argument('--out',default='generated/run_storage.tsv')
    args=ap.parse_args()
    plans=json.loads(Path(args.run_plan).read_text(encoding='utf-8'))
    rows=[]
    for plan in plans:
        snapshot={'checked_at_utc':datetime.now(timezone.utc).isoformat(),'run_id':plan['run_id'],
                  'library_id':plan['library_id'],'kit_id':plan['kit_id'],'technology':plan['technology'],
                  'library_context':plan.get('library_context'),
                  'areas':{key:measure(plan['paths'][key]) for key in ['sra','fastq','merged','tmp','output']},
                  'note':'目录大小是当前保留文件量；不等于累计网络下载量，也不证明MD5/gzip/FASTQ校验通过。硬链接仅在同一目录快照内去重，符号链接不跟随。'}
        stats=Path(plan['paths']['stats'])
        if stats.is_dir():
            (stats/'run_storage.json').write_text(json.dumps(snapshot,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        for key,result in snapshot['areas'].items():
            rows.append({'run_id':plan['run_id'],'kit':plan['kit_id'],'technology':plan['technology'],
                         'area':key,'exists':result['exists'],'files':result['files'],
                         'logical_bytes':result['logical_bytes'],'unique_inode_bytes':result['unique_inode_bytes'],
                         'symlinks':result['symlinks'],'complete':result.get('complete',False),
                         'error':';'.join(result['errors'])})
    out=Path(args.out);out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else ['run_id'],delimiter='\t')
        w.writeheader();w.writerows(rows)
    print('REPORT',out,'runs',len(plans))

if __name__=='__main__':main()
