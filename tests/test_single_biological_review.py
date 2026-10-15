"""Regression for zero-length technical + one biological read, not missing mates."""
import copy
import io
import json
import os
import sys
import tempfile
import threading
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
import prepare_layout as layout
import run_download_merge as d
import review_srr19908303_layout as review


class ReviewTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='cts-read-layout-')
        self.base=Path(self.temp.name)
        self.repo=self.base/'repo'
        for name in ['config','generated','runtime','bin']:(self.repo/name).mkdir(parents=True)
        self.site=layout.read_json(ROOT/'config/site.json')
        self.site.update(data_root=str(self.base/'data'),min_free_gb=0)
        registry=layout.read_json(ROOT/'config/kit_registry.json')
        refs=layout.read_json(ROOT/'config/reference_registry.json')
        row={'目标Kit':'SPI-ISR','归一技术':'ChAR-seq','Run ID':review.RUN,
             'Source data ID':review.LIBRARY,'Species':'Homo sapiens','官方物种':'Homo sapiens',
             '官方平台':'ILLUMINA','官方布局':'PAIRED','样本表内Run数':'1'}
        self.plan=layout.plan_run(row,self.site,registry,refs)
        self.libraries=layout.build_libraries([self.plan],self.site,registry)
        layout.materialize([self.plan],self.libraries,self.site)
        d.save(self.repo/'config/site.json',self.site)
        d.save(self.repo/'generated/run_plan.json',[self.plan])
        d.save(self.repo/'generated/library_plan.json',self.libraries)
        self.result=[{'run_id':review.RUN,'status':'FAILED_OR_NEEDS_REVIEW','error':'original role hold'}]
        d.save(self.repo/'generated/run_results.json',self.result)
        self.payload=Path(self.plan['paths']['sra_accession_dir'])/review.RUN/(review.RUN+'.sra')
        self.payload.parent.mkdir(parents=True)
        self.payload.write_bytes(b'mock SRA')
        self.fastq=Path(self.plan['paths']['fastq'])/(review.RUN+'.fastq.gz')
        import gzip
        with gzip.open(self.fastq,'wb') as f:
            for i in range(1,4):f.write(('@'+review.RUN+'.'+str(i)+'/2\nACGT\n+\nIIII\n').encode())
        self.record=d.scan_fastq(self.fastq,review.RUN)
        self.record['role']='single'
        state=d.read(self.plan['paths']['state'])
        state['download']={'status':'VALIDATED_SRA','vdb_validate_pass':True,
                           'path':str(self.payload),'fingerprint':d.fp(self.payload)}
        state['validation']={'status':'FASTQ_VALIDATED','input_files':[self.record]}
        state['input_preparation']={'status':'FAILED_OR_NEEDS_REVIEW','error':'role mismatch'}
        state['error']={'message':'role mismatch'}
        d.save(self.plan['paths']['state'],state)
        self.mock=self.repo/'bin/vdb-dump'
        self.mock.write_text('#!'+sys.executable+'\nimport sys\n'
                             'assert "-R" not in sys.argv, "must check every spot"\n'
                             'for _ in range(2):print("0, 100\\tSRA_READ_TYPE_TECHNICAL, SRA_READ_TYPE_BIOLOGICAL")\n'
                             'print("100, 0\\tSRA_READ_TYPE_BIOLOGICAL, SRA_READ_TYPE_TECHNICAL")\n')
        self.mock.chmod(0o755)
        d.save(self.repo/'runtime/runtime_bin.json',{'bin_dir':str(self.mock.parent)})
        self.old_root=d.ROOT

    def tearDown(self):
        d.ROOT=self.old_root
        self.temp.cleanup()

    def test_repair_and_resume_without_download_or_conversion(self):
        # Unreviewed singleton must retain its existing layout hold.
        with self.assertRaisesRegex(RuntimeError,'READ_ROLE_REVIEW'):
            d.input_roles([copy.deepcopy(self.record)],'PAIRED')
        before=self.fastq.read_bytes(),self.payload.read_bytes()
        receipt=review.repair(self.repo)
        self.assertEqual(receipt['status'],'RUN_INPUT_READY')
        self.assertEqual(d.read(self.plan['paths']['state'])['input_preparation']['status'],'RUN_INPUT_READY')
        self.assertNotIn('error',d.read(self.plan['paths']['state']))
        ready=d.read(self.plan['paths']['ready'])
        self.assertEqual((ready['layout'],ready['declared_layout']),('SINGLE','PAIRED'))
        self.assertFalse(ready['calculation_submission_allowed'])
        audit=d.read(Path(self.plan['paths']['stats'])/'read_layout_review.json')
        self.assertEqual(audit['biological_read_positions'],{'1':1,'2':2})
        self.assertEqual(d.read(Path(self.libraries[0]['library_root'])/'library_ready.json')['layout'],'SINGLE')
        self.assertEqual(d.read(self.repo/'generated/run_results.json'),self.result)
        self.assertEqual(d.read(self.repo/'generated/run_plan.json')[0]['layout'],'PAIRED')
        self.assertEqual((self.fastq.read_bytes(),self.payload.read_bytes()),before)
        # The persisted audit also survives future ordinary dispatches.
        result=d.process_run(self.plan,self.site,'all',threading.Semaphore(1))
        self.assertEqual(result['status'],'RUN_INPUT_READY_REUSED')
        # Idempotent retry reuses proof but checks compressed bytes again.
        self.mock.write_text('#!'+sys.executable+'\nraise SystemExit(99)\n')
        self.assertEqual(review.repair(self.repo)['status'],'RUN_INPUT_READY')

    def test_nonmatching_late_spot_and_tool_failure_block_publication(self):
        self.mock.write_text('#!'+sys.executable+'\n'
                             'print("0, 100\\tSRA_READ_TYPE_TECHNICAL, SRA_READ_TYPE_BIOLOGICAL")\n'
                             'print("100, 100\\tSRA_READ_TYPE_BIOLOGICAL, SRA_READ_TYPE_BIOLOGICAL")\n')
        with self.assertRaisesRegex(RuntimeError,'Different read structure'):
            review.repair(self.repo)
        self.assertFalse(Path(self.plan['paths']['ready']).exists())
        self.mock.write_text('#!'+sys.executable+'\n'
                             'for _ in range(3):print("0, 100\\tSRA_READ_TYPE_TECHNICAL, SRA_READ_TYPE_BIOLOGICAL")\n'
                             'raise SystemExit(7)\n')
        with self.assertRaisesRegex(RuntimeError,'vdb-dump failed'):
            review.repair(self.repo)
        self.assertFalse(Path(self.plan['paths']['ready']).exists())

    def test_changed_source_and_fastq_block_reuse(self):
        review.repair(self.repo)
        state=d.read(self.plan['paths']['state'])
        self.payload.write_bytes(b'changed source')
        with self.assertRaisesRegex(RuntimeError,'STALE_OR_INCOMPLETE'):
            d.reviewed_input_plan(self.plan,state,self.site)
        with self.assertRaisesRegex(RuntimeError,'stale'):
            review.repair(self.repo)

    def test_hash_mismatch_blocks_ready_even_with_forged_validation_hash(self):
        state=d.read(self.plan['paths']['state'])
        state['validation']['input_files'][0]['gzip_sha256']='0'*64
        d.save(self.plan['paths']['state'],state)
        with self.assertRaisesRegex(RuntimeError,'SHA256 differs'):
            review.repair(self.repo)
        self.assertFalse(Path(self.plan['paths']['ready']).exists())

    def test_run_lock_is_respected(self):
        import fcntl
        with (Path(self.plan['paths']['run_root'])/'.download.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            with self.assertRaisesRegex(RuntimeError,'busy'):
                review.repair(self.repo)


class RowTests(unittest.TestCase):
    def test_counts_types_empty_and_unexpected_output(self):
        good='0, 100\tSRA_READ_TYPE_TECHNICAL, SRA_READ_TYPE_BIOLOGICAL\n'
        reverse='100, 0\tSRA_READ_TYPE_BIOLOGICAL, SRA_READ_TYPE_TECHNICAL\n'
        # Mirrors the reported late switch at spot 16,760,307: slot order
        # changes, but every spot still contains exactly one biological read.
        result=review.scan_rows(io.StringIO(good*2+reverse),3)
        self.assertEqual(result['all_spots_checked'],3)
        self.assertEqual(result['biological_read_positions'],{'1':1,'2':2})
        self.assertEqual(review.scan_rows(io.StringIO(reverse*3),3)['biological_read_positions'],{'1':3})
        for bad in ['',good,good+reverse,
                    good+'0, 100\tSRA_READ_TYPE_BIOLOGICAL, SRA_READ_TYPE_BIOLOGICAL\n',
                    good*2+'100, 100\tSRA_READ_TYPE_BIOLOGICAL, SRA_READ_TYPE_BIOLOGICAL\n',
                    good*2+'0, 0\tSRA_READ_TYPE_TECHNICAL, SRA_READ_TYPE_BIOLOGICAL\n',
                    good*2+'100, 100\tSRA_READ_TYPE_TECHNICAL, SRA_READ_TYPE_BIOLOGICAL\n',
                    good*2+'0, 100\tSRA_READ_TYPE_TECHNICAL, SRA_READ_TYPE_TECHNICAL\n',
                    good*2+'100, -1\tSRA_READ_TYPE_BIOLOGICAL, SRA_READ_TYPE_TECHNICAL\n',
                    good*2+'0, 100, 100\tSRA_READ_TYPE_TECHNICAL, SRA_READ_TYPE_BIOLOGICAL, SRA_READ_TYPE_BIOLOGICAL\n',
                    good+'invalid header\n']:
            with self.assertRaises(RuntimeError):review.scan_rows(io.StringIO(bad),3)


if __name__=='__main__':unittest.main()
