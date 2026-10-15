#!/usr/bin/env python3
"""Prepare Run-oriented storage and resolve registered species/Kit routes.

This program writes metadata/plans and optionally creates directories. It does
not download, merge, launch a shell, or submit calculation jobs.
"""
import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


def component(value):
    value = str(value or '').strip()
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', value) or value in {'.', '..'}:
        raise ValueError('Unsafe directory component: ' + repr(value))
    return value


def normalize_technology(value, registry):
    lookup = re.sub(r'[\s_-]+', '', str(value).strip()).casefold()
    matches = []
    for canonical, recipe in registry['technology_routes'].items():
        aliases = [canonical] + recipe.get('aliases', [])
        if any(re.sub(r'[\s_-]+', '', alias).casefold() == lookup for alias in aliases):
            matches.append((canonical, recipe))
    if len(matches) != 1:
        raise ValueError('Technology has no unique route: ' + str(value))
    return matches[0]


def resolve_species(row, registry):
    # Official record takes precedence; original taxonomy is retained as a hint.
    official = str(row.get('官方物种') or '').strip()
    hint = str(row.get('Species') or '').strip()
    found = [(key, item) for key, item in registry['species'].items()
             if official and official in [item['scientific_name']] + item.get('aliases', [])]
    if len(found) == 1:
        taxid, profile = found[0]
        conflict = bool(hint and hint not in [profile['scientific_name']] + profile.get('aliases', []))
        return {'tax_id': taxid, 'scientific_name': profile['scientific_name'],
                'species_hint': hint, 'taxonomy_confirmed': not conflict,
                'taxonomy_source': 'official_metadata', 'conflict': conflict}
    return {'tax_id': None, 'scientific_name': official or hint,
            'species_hint': hint, 'taxonomy_confirmed': False,
            'taxonomy_source': 'unresolved', 'conflict': bool(official and hint and official != hint)}


def required_asset_error(name, value, asset_type, check_assets):
    if not value:
        return 'REFERENCE_ASSET_UNCONFIGURED:' + name
    if not check_assets:
        return None
    p = Path(value)
    if asset_type == 'directory':
        return None if p.is_dir() else 'REFERENCE_DIRECTORY_MISSING:' + name
    if asset_type == 'file':
        return None if p.is_file() and p.stat().st_size else 'REFERENCE_FILE_MISSING:' + name
    if asset_type == 'intron_prefix':
        required = [Path(str(p) + ending) for ending in ['.bed', '_part1.bed', '_part2.bed']]
        return None if all(q.is_file() and q.stat().st_size for q in required) else 'INTRON_BED_SET_INCOMPLETE:' + name
    if asset_type == 'bowtie2_prefix':
        for ext in ['bt2', 'bt2l']:
            required = [Path(str(p) + '.' + ending + '.' + ext)
                        for ending in ['1', '2', '3', '4', 'rev.1', 'rev.2']]
            if all(q.is_file() and q.stat().st_size for q in required):
                return None
        return 'BOWTIE2_INDEX_INCOMPLETE:' + name
    return 'UNKNOWN_REFERENCE_ASSET_TYPE:' + name


def select_reference(species, kit_id, kits, references, check_assets):
    errors = []
    reference = None
    if not species['taxonomy_confirmed']:
        errors.append('TAXONOMY_UNCONFIRMED')
    else:
        ref_id = references.get('defaults', {}).get(species['tax_id'])
        reference = references.get('bundles', {}).get(ref_id) if ref_id else None
        if reference is None:
            errors.append('REFERENCE_NOT_REGISTERED')
        else:
            if str(reference.get('tax_id')) != species['tax_id']:
                errors.append('REFERENCE_TAXON_MISMATCH')
            if not reference.get('validated'):
                errors.append('REFERENCE_NOT_VALIDATED')
            if not reference.get('assembly_accession') or not reference.get('annotation_release'):
                errors.append('REFERENCE_VERSION_INCOMPLETE')
            if kit_id not in reference.get('compatible_kits', []):
                errors.append('REFERENCE_KIT_COMPATIBILITY_UNCONFIRMED')
            for name, kind in kits[kit_id].get('required_assets', {}).items():
                error = required_asset_error(name, reference.get('assets', {}).get(name), kind, check_assets)
                if error:
                    errors.append(error)
    return reference, errors


def plan_run(row, site, registry, references, check_assets=False):
    kits = registry['kits']
    canonical, technology = normalize_technology(row['归一技术'], registry)
    platform = str(row.get('官方平台') or '').strip()
    kit_id = str(row.get('目标Kit') or '').strip()
    selection_source = 'manifest_requested_kit' if kit_id else 'technology_and_platform'
    if not kit_id:
        candidates = [k for k in technology['allowed_kits']
                      if not platform or platform in kits[k]['allowed_platforms']]
        if len(candidates) != 1:
            raise ValueError('No unique Kit for technology/platform: ' + canonical + '/' + platform)
        kit_id = candidates[0]
    if kit_id not in kits:
        raise ValueError('Unknown Kit: ' + kit_id)
    kit = kits[kit_id]
    run = component(row['Run ID'])
    sample = component(row['Source data ID'])
    root = PurePosixPath(site['data_root'])
    if not root.is_absolute():
        raise ValueError('data_root must be an absolute Linux path')
    run_root = root / component(kit['directory_name']) / component(technology['directory_name']) / run
    paths = {key: str(run_root / key) for key in ['sra', 'fastq', 'merged', 'tmp', 'stats', 'logs', 'output']}
    paths['run_root'] = str(run_root)
    paths['sra_accession_dir'] = str(run_root / 'sra' / run)
    paths['metadata'] = str(run_root / 'metadata.json')
    paths['state'] = str(run_root / 'state.json')
    paths['routing'] = str(run_root / 'routing.json')
    paths['ready'] = str(run_root / 'ready.json')
    species = resolve_species(row, registry)
    errors = []
    assay_confirmed = row.get('biological_review_status') == 'CONFIRMED' and bool(row.get('biological_review_evidence_url'))
    if not assay_confirmed:
        errors.append('ASSAY_REVIEW_UNCONFIRMED')
    if row.get('special_calculation_hold'):
        errors.append('SPECIAL_INPUT_ADAPTER_UNCONFIRMED')
    if row.get('special_input_profile')=='MEL_HBB_targeted_custom_reference':
        errors.append('TARGETED_HBB_CUSTOM_REFERENCE_REQUIRED')
    if kit_id not in technology['allowed_kits']:
        errors.append('TECHNOLOGY_KIT_CONFLICT')
    if not platform:
        errors.append('PLATFORM_UNCONFIRMED')
    elif platform not in kit['allowed_platforms']:
        errors.append('PLATFORM_KIT_CONFLICT')
    if row.get('官方布局') not in ['SINGLE', 'PAIRED']:
        errors.append('LAYOUT_UNCONFIRMED')
    if not kit.get('enabled'):
        errors.append('KIT_NOT_ENABLED')
    if not kit.get('version') or not kit.get('entrypoint') or not kit.get('adapter_config'):
        errors.append('KIT_RUNTIME_NOT_REGISTERED')
    if not kit.get('required_assets'):
        errors.append('KIT_REFERENCE_CONTRACT_NOT_REGISTERED')
    preset = technology.get('presets', {}).get(kit_id)
    if not preset or preset not in kit.get('supported_presets', []):
        errors.append('KIT_PRESET_NOT_REGISTERED')
    if check_assets and kit.get('enabled'):
        for field in ['entrypoint', 'adapter_config']:
            if not kit.get(field) or not Path(kit[field]).is_file():
                errors.append('KIT_FILE_MISSING:' + field)
        if kit.get('container_image') and not Path(kit['container_image']).is_file():
            errors.append('CONTAINER_IMAGE_MISSING')
    reference, reference_errors = select_reference(species, kit_id, kits, references, check_assets)
    errors.extend(reference_errors)
    urls = [u for u in str(row.get('ENA_FASTQ地址') or '').split(';') if u]
    hashes = str(row.get('ENA_FASTQ_MD5') or '').split(';')
    sizes = str(row.get('ENA_FASTQ字节数') or '').split(';')
    ena_files = []
    for i, url in enumerate(urls):
        url = 'https://' + url[6:] if url.startswith('ftp://') else url if url.startswith(('http://', 'https://')) else 'https://' + url
        filename = component(url.rsplit('/', 1)[-1])
        ena_files.append({'url': url, 'expected_bytes': int(sizes[i]) if i < len(sizes) and sizes[i].isdigit() else None,
                          'expected_md5': hashes[i] if i < len(hashes) and re.fullmatch('[0-9a-fA-F]{32}', hashes[i]) else None,
                          'destination': str(run_root / 'fastq' / filename),
                          'checksum_scope': 'ENA_archive_generated_gzip_bytes'})
    mode = 'GSA_MANUAL_CHANNEL' if run.startswith('CRR') else 'SRA_PREFETCH'
    if row.get('download_execution_mode')=='ENA_ARCHIVE_FASTQ':
        if run.startswith('CRR') or not ena_files or any(not f['expected_bytes'] or not f['expected_md5'] for f in ena_files):
            raise ValueError('Exact ENA download needs public SRA Run and complete source fingerprints')
        mode='ENA_ARCHIVE_FASTQ'
    return {'schema_version': '1.0', 'run_id': run, 'library_id': sample,
            'kit_id': kit_id, 'technology': canonical,
            'original_technology': row.get('Technology'), 'original_kit': row.get('Kit'),
            'research_id': row.get('Source research ID'), 'species': species,
            'platform': platform or None, 'layout': row.get('官方布局'),
            'expected_library_runs': int(row['样本表内Run数']), 'paths': paths,
            'download': {'mode': mode, 'ena_files': ena_files,
                         'official_metadata_url': row.get('ENA证据URL'),
                         'sra_validation': 'vdb-validate',
                         'converted_fastq_validation': ['full_gzip_read', 'FASTQ_format', 'layout_and_read_role', 'pair_counts_if_applicable'],
                         'compare_converted_fastq_to_ENA_MD5': False},
            'routing': {'kit_id': kit_id, 'kit_version': kit.get('version'),
                        'kit_selection_source': selection_source,
                        'entrypoint': kit.get('entrypoint'), 'adapter_config': kit.get('adapter_config'),
                        'container_image': kit.get('container_image'), 'preset': preset,
                        'reference_id': reference.get('reference_id') if reference else None,
                        'reference_bundle': reference, 'asset_check_performed': check_assets,
                        'route_status': 'CONFIGURATION_READY' if not errors else 'HELD',
                        'hold_reasons': sorted(set(errors)),
                        'calculation_submission_allowed': False},
            'source_status': 'DOWNLOAD_HISTORY_UNKNOWN',
            'input_validation': {'passed': False, 'records': []},
            'review': {'library_merge_approved': False, 'previous_run_duplicates_checked': False,
                       'biological_review_confirmed': assay_confirmed,
                       'biological_review_evidence_url': row.get('biological_review_evidence_url'),
                       'notes': row.get('运行前核对')},
            'selection_batch':row.get('selection_batch'),
            'special_input_profile':row.get('special_input_profile'),
            'special_calculation_hold':row.get('special_calculation_hold'),
            'source_sheet_row': row.get('来源sheet行号')}


def build_libraries(plans, site, registry):
    groups = defaultdict(list)
    for plan in plans:
        groups[plan['library_id']].append(plan)
    libraries = []
    for sample, members in sorted(groups.items()):
        kits = {p['kit_id'] for p in members}
        techs = {p['technology'] for p in members}
        if len(kits) != 1 or len(techs) != 1:
            raise ValueError('Library has conflicting Kit/technology: ' + sample)
        first = members[0]
        kit_dir = registry['kits'][first['kit_id']]['directory_name']
        _, tech = normalize_technology(first['technology'], registry)
        library_root = (PurePosixPath(first['paths']['run_root']) if len(members) == 1 else
                        PurePosixPath(site['data_root']) / kit_dir / tech['directory_name'] / '_libraries' / component(sample))
        reasons = sorted({reason for p in members for reason in p['routing']['hold_reasons']})
        taxids = {p['species']['tax_id'] for p in members}
        reference_ids = {p['routing']['reference_id'] for p in members}
        if len(taxids) != 1 or len(reference_ids) != 1:
            reasons.append('LIBRARY_REFERENCE_CONFLICT')
        if any(p['expected_library_runs'] != len(members) for p in members):
            reasons.append('LIBRARY_RUN_SET_INCOMPLETE')
        if len(members) > 1:
            reasons.append('MULTI_RUN_MERGE_NEEDS_REVIEW')
        reasons.append('SOURCE_INPUTS_NOT_VALIDATED')
        libraries.append({'schema_version': '1.0', 'library_id': sample,
                          'kit_id': first['kit_id'], 'technology': first['technology'],
                          'run_ids': sorted(p['run_id'] for p in members),
                          'run_metadata_paths': [p['paths']['metadata'] for p in members],
                          'tax_id': next(iter(taxids)) if len(taxids) == 1 else None,
                          'reference_id': next(iter(reference_ids)) if len(reference_ids) == 1 else None,
                          'compute_unit': site.get('compute_unit', 'library'),
                          'paths': {k: str(library_root / k) for k in ['merged', 'tmp', 'stats', 'logs', 'output']},
                          'library_root': str(library_root), 'merge_approved': len(members) == 1,
                          'calculation_submission_allowed': False,
                          'hold_reasons': sorted(set(reasons))})
    return libraries


def write_tsv(path, rows, fields):
    with Path(path).open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def initial_run_state(run_id):
    return {'schema_version':'1.0','run_id':run_id,
            'download':{'status':'PLANNED','source_mode':None,'started_at':None,'finished_at':None,'bytes':None,'error':None},
            'validation':{'status':'NOT_RUN','vdb_validate':None,'official_MD5_for_ENA_file':None,
                          'gzip':None,'FASTQ_format':None,'input_files':[],'read_roles':[],
                          'validated_input_fingerprints':[]},
            'input_preparation':{'status':'NOT_RUN','merge_manifest':None,'tmp_bytes':None},
            'calculation':{'status':'NOT_SUBMITTED','kit_version':None,'reference_id':None,'job_id':None,'error':None}}


def materialize(plans, libraries, site):
    data_root = Path(site['data_root']).resolve()
    for plan in plans:
        run_root = Path(plan['paths']['run_root']).resolve()
        if data_root not in run_root.parents:
            raise ValueError('Run path escapes configured data root')
        for key in ['sra', 'fastq', 'merged', 'tmp', 'stats', 'logs', 'output']:
            Path(plan['paths'][key]).mkdir(parents=True, exist_ok=True)
        metadata = Path(plan['paths']['metadata'])
        if metadata.exists():
            existing = read_json(metadata)
            identity = ('run_id', 'library_id', 'kit_id', 'technology')
            if any(existing.get(k) != plan.get(k) for k in identity):
                raise ValueError('Existing Run metadata identity conflicts: ' + str(metadata))
        else:
            write_json(metadata, {k:v for k,v in plan.items() if k not in ['routing', 'input_validation']})
        write_json(plan['paths']['routing'], plan['routing'])
        state = Path(plan['paths']['state'])
        if not state.exists():
            write_json(state, initial_run_state(plan['run_id']))
    for library in libraries:
        library_root = Path(library['library_root']).resolve()
        if data_root not in library_root.parents:
            raise ValueError('Library path escapes configured data root')
        for key in ['merged', 'tmp', 'stats', 'logs', 'output']:
            Path(library['paths'][key]).mkdir(parents=True, exist_ok=True)
        metadata = library_root / 'library_metadata.json'
        if not metadata.exists():
            write_json(metadata, library)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--site', default='config/site.json')
    parser.add_argument('--runs', default='manifests/selected_runs.tsv')
    parser.add_argument('--kit-registry', default='config/kit_registry.json')
    parser.add_argument('--references', default='config/reference_registry.json')
    parser.add_argument('--out-dir', default='generated')
    parser.add_argument('--materialize', action='store_true', help='Create directories and initialize missing metadata/state files only')
    parser.add_argument('--check-assets', action='store_true', help='Check configured Kit/reference paths on this server')
    args = parser.parse_args()
    site, registry, references = map(read_json, [args.site, args.kit_registry, args.references])
    with Path(args.runs).open(encoding='utf-8-sig', newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    ids = [r['Run ID'] for r in rows]
    if len(ids) != len(set(ids)):
        raise ValueError('Duplicate Run IDs in manifest')
    plans, planning_errors = [], []
    for row in rows:
        try:
            plans.append(plan_run(row, site, registry, references, args.check_assets))
        except ValueError as error:
            planning_errors.append({'run_id':row.get('Run ID'), 'library_id':row.get('Source data ID'),
                                    'kit':row.get('目标Kit'), 'technology':row.get('归一技术'), 'error':str(error)})
    libraries = build_libraries(plans, site, registry)
    by_library = {item['library_id']:item for item in libraries}
    for plan in plans:
        lib=by_library[plan['library_id']]
        plan['library_context']={'compute_unit':site.get('compute_unit','library'),
                                 'library_id':lib['library_id'],'library_root':lib['library_root'],
                                 'library_run_ids':lib['run_ids'],'shared_across_runs':len(lib['run_ids'])>1}
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    write_json(out / 'run_plan.json', plans)
    write_json(out / 'library_plan.json', libraries)
    write_tsv(out / 'planning_errors.tsv', planning_errors, ['run_id','library_id','kit','technology','error'])
    write_tsv(out / 'run_plan.tsv', [
        {'run_id': p['run_id'], 'library_id': p['library_id'], 'kit': p['kit_id'],
         'technology': p['technology'], 'tax_id': p['species']['tax_id'],
         'platform': p['platform'], 'layout': p['layout'], 'run_root': p['paths']['run_root'],
         'download_mode': p['download']['mode'], 'reference_id': p['routing']['reference_id'],
         'route_status': p['routing']['route_status'], 'hold_reasons': ';'.join(p['routing']['hold_reasons'])}
        for p in plans], ['run_id', 'library_id', 'kit', 'technology', 'tax_id', 'platform', 'layout', 'run_root', 'download_mode', 'reference_id', 'route_status', 'hold_reasons'])
    summary = {'prepared_at_utc': datetime.now(timezone.utc).isoformat(), 'runs': len(plans),
               'libraries': len(libraries), 'data_root': site['data_root'],
               'run_route_status': dict(Counter(p['routing']['route_status'] for p in plans)),
               'download_channels': dict(Counter(p['download']['mode'] for p in plans)),
               'compute_unit': site.get('compute_unit', 'library'),
               'unresolved_manifest_records': len(planning_errors),
               'run_manifest_sha256': hashlib.sha256(Path(args.runs).read_bytes()).hexdigest(),
               'directories_created': args.materialize, 'download_or_compute_started': False}
    write_json(out / 'summary.json', summary)
    if args.materialize:
        materialize(plans, libraries, site)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
