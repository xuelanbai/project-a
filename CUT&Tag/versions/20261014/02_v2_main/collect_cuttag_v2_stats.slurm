#!/bin/bash -l
#SBATCH --job-name=cuttag_v2_stats
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/mnt/nas/swine/home/jxwu/Projects/A/08.CUT&Tag/00.26.09/analysis_nfcore_cutandrun_3.2.2/logs/cuttag-v2-stats-%j.out
set -euo pipefail
PROJECT="$HOME/Projects/A/08.CUT&Tag/00.26.09/analysis_nfcore_cutandrun_3.2.2"
PYTHON="$HOME/apps/cutandtag/conda/envs/env-0d01150023c68e439ebff3eb959d98bf/bin/python"
[[ -x "$PYTHON" ]] || { echo "Missing Python: $PYTHON" >&2; exit 1; }
[[ -n "${SLURM_JOB_ID:-}" && "$(hostname)" != mn-* ]] || { echo "Submit this script with sbatch." >&2; exit 1; }
export PYTHONIOENCODING=utf-8
"$PYTHON" - "$PROJECT" <<'PY'
import bisect
import csv
import datetime
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import sys

def bed_records(path):
    records = []
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'rt', encoding='utf-8') as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            fields = line.rstrip('\n').split('\t')
            try:
                start, end = int(fields[1]), int(fields[2])
                if start < 0 or end <= start:
                    raise ValueError('invalid interval')
            except (ValueError, IndexError) as exc:
                raise ValueError('{}:{} invalid BED: {}'.format(path, line_no, exc))
            records.append((fields[0], start, end, fields))
    return records

def merged(records):
    result = {}
    for chrom, start, end, _ in sorted(records, key=lambda r: (r[0], r[1], r[2])):
        intervals = result.setdefault(chrom, [])
        if intervals and start <= intervals[-1][1]:
            intervals[-1][1] = max(end, intervals[-1][1])
        else:
            intervals.append([start, end])
    return result

def shared_intervals(a, b):
    am, bm = merged(a), merged(b)
    shared = []
    for chrom in sorted(am.keys() & bm.keys()):
        left, right = am[chrom], bm[chrom]
        i = j = 0
        while i < len(left) and j < len(right):
            start, end = max(left[i][0], right[j][0]), min(left[i][1], right[j][1])
            if start < end:
                shared.append((chrom, start, end, []))
            if left[i][1] < right[j][1]:
                i += 1
            else:
                j += 1
    # Match BEDTools merge's joining of adjacent covered intervals.
    return [(c, s, e) for c, intervals in merged(shared).items() for s, e in intervals]

def supported_count(a, b):
    bm = merged(b)
    starts = {c: [s for s, e in intervals] for c, intervals in bm.items()}
    count = 0
    for chrom, start, end, _ in a:
        if chrom in bm:
            idx = bisect.bisect_left(starts[chrom], end) - 1
            if idx >= 0 and bm[chrom][idx][1] > start:
                count += 1
    return count

def pq_pass(record):
    fields = record[3]
    if len(fields) < 9:
        raise ValueError('Expected MACS2 narrowPeak p/q columns 8 and 9')
    return float(fields[7]) >= 5 and float(fields[8]) >= 2

def main():
    project = Path(sys.argv[1]).expanduser().resolve()
    results = project / 'results_pig_v2'
    if not results.is_dir():
        raise SystemExit('Missing new results: {}'.format(results))
    dest = project / 'logs' / 'cuttag_v2_table_inputs.json.gz'
    document = {
        'schema': 'cuttag_v2_table_inputs_v1',
        'collected_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'project': str(project), 'results': str(results),
        'definitions': {
            'raw_peaks': 'Number of called BED/narrowPeak intervals; empty existing file is zero; missing file is null.',
            'filtered_macs2_peaks': 'narrowPeak column8 >= 5 AND column9 >= 2 (-log10 P and Q).',
            'reproducibility': 'At least one overlapping bp between independent replicate calls; shared genomic bases unioned into non-overlapping intervals.',
            'pooled_peaks': 'MACS2 calls on jointly supplied replicate BAMs; separate from replicate-supported peaks.',
            'TSS': 'Profiles collected if present; no arbitrary enrichment cutoff or center/flank ratio is assigned.',
            'fragment_length': 'Existing QC sources only; not estimated from read length or pooled samples.'
        },
        'samples': [], 'peak_files': [], 'reproducibility': [], 'merged_tracks': [],
        'inventory': [], 'text_sources': {}, 'missing': [], 'skipped': [], 'errors': []
    }
    total_text = 0
    def capture(path, limit=4*1024*1024, tail=False):
        nonlocal total_text
        path = Path(path)
        if not path.is_file():
            return
        key = str(path.relative_to(project)) if path.is_relative_to(project) else str(path)
        if key in document['text_sources']:
            return
        size = path.stat().st_size
        if (size > limit and not tail) or total_text + min(size, limit) > 150*1024*1024:
            document['skipped'].append({'path': key, 'reason': 'text_size_limit', 'bytes': size})
            return
        with path.open('rb') as handle:
            if tail and size > limit:
                handle.seek(size-limit)
            data = handle.read(limit)
        document['text_sources'][key] = {
            'bytes': size, 'tail_only': bool(tail and size > limit),
            'sha256_captured_bytes': hashlib.sha256(data).hexdigest(),
            'text': data.decode('utf-8', errors='replace')
        }
        total_text += len(data)

    samplesheet = project / 'config/samplesheet.csv'
    capture(samplesheet)
    with samplesheet.open(encoding='utf-8-sig', newline='') as handle:
        document['samples'] = list(csv.DictReader(handle))
    for sample in document['samples']:
        for read in ['fastq_1', 'fastq_2']:
            p = Path(sample[read])
            sample[read+'_bytes'] = p.stat().st_size if p.is_file() else None

    for relative in ['config/params.yaml', 'config/run_cutandrun_pig_v2.slurm',
                     'config/postprocess_cutandtag_pig_v2.slurm', 'config/run_cuttag_tss_v2.slurm']:
        capture(project / relative)
    for p in (project / 'logs').glob('pig-v2-gsize-*.txt'):
        capture(p)
    for pattern in ['cuttag-pig-v2-*.out', 'cuttag-tracks-v2-*.out', 'cuttag-tss-v2-*.out']:
        for p in (project / 'logs').glob(pattern):
            capture(p, limit=16000, tail=True)
    for p in (project / 'results/pipeline_info').glob('*params*'):
        if p.suffix.lower() in ['.json', '.yaml', '.yml', '.txt']:
            capture(p)

    # Scan only published new results. Never read FASTQ, BAM, bigWig or genome sequences.
    for directory, subdirs, filenames in os.walk(results, followlinks=False):
        subdirs[:] = [d for d in subdirs if not d.startswith('.') and d != '00_genome']
        for name in sorted(filenames):
            p = Path(directory) / name
            if not p.is_file():
                document['errors'].append({'path': str(p), 'reason': 'broken_symlink_or_missing_file'})
                continue
            rel = p.relative_to(results).as_posix()
            low = name.lower()
            size = p.stat().st_size
            document['inventory'].append({'path': rel, 'bytes': size})
            in_peak_dir = rel.startswith(('03_peak_calling/04_called_peaks/',
                '03_peak_calling/05_consensus_peaks/', 'custom_tracks/reproducible_peaks_anyOverlap/'))
            is_peak = low.endswith(('.narrowpeak', '.narrowpeak.gz', '.broadpeak', '.bed', '.bed.gz'))
            if (in_peak_dir and is_peak) or (rel.startswith('custom_tracks/macs2_ppois/') and low.endswith('.narrowpeak')):
                try:
                    records = bed_records(p)
                    is_macs = '.narrowpeak' in low
                    document['peak_files'].append({'path': rel, 'raw_peaks': len(records),
                        'filtered_p5_q2': sum(pq_pass(r) for r in records) if is_macs else None})
                except Exception as exc:
                    document['errors'].append({'path': rel, 'reason': str(exc)})
            suffix = p.suffix.lower()
            selected_text = suffix in ['.txt', '.tsv', '.csv', '.json', '.yaml', '.yml', '.log', '.out', '.stats', '.flagstat']
            if selected_text and (rel.startswith(('01_', '02_', '03_', '04_', 'pipeline_info/', 'custom_metadata/', 'custom_tracks/')) or 'multiqc' in rel.lower()):
                capture(p, limit=20*1024*1024 if low == 'multiqc_data.json' else 4*1024*1024)
            if low.endswith('_peaks.xls'):
                # MACS2 XLS is plain text. Include arguments only, not duplicate all peaks.
                with p.open(encoding='utf-8') as handle:
                    header = []
                    for line in handle:
                        if line.startswith('#') or not line.strip():
                            header.append(line)
                        else:
                            break
                document['text_sources'][str(p.relative_to(project))] = {'text': ''.join(header), 'header_only': True}

    for state in ['119', '137']:
        for mark in ['Pol2', 'Tag']:
            group = '3D4_{}_{}'.format(state, mark)
            for caller, directory, ending in [('MACS2','macs2','.macs2_peaks.narrowPeak'), ('SEACR','seacr','.seacr.peaks.stringent.bed')]:
                pair = [results/'03_peak_calling/04_called_peaks'/directory/(group+'_R'+str(rep)+ending) for rep in [1,2]]
                if not all(p.is_file() for p in pair):
                    document['missing'].extend(str(p) for p in pair if not p.is_file())
                    continue
                a, b = map(bed_records, pair)
                for label, aa, bb in [(caller,a,b)] + ([(caller+'_p5q2',[r for r in a if pq_pass(r)],[r for r in b if pq_pass(r)])] if caller=='MACS2' else []):
                    shared = shared_intervals(aa, bb)
                    document['reproducibility'].append({'group': group, 'caller': label,
                        'control': '3D4_{}_IgG_R1'.format(state),
                        'rep1_called': len(aa), 'rep2_called': len(bb),
                        'rep1_supported': supported_count(aa,bb), 'rep2_supported': supported_count(bb,aa),
                        'shared_intervals': len(shared), 'shared_bp': sum(e-s for c,s,e in shared),
                        'sources': [str(p) for p in pair]})
            bw = results/'custom_tracks/macs2_ppois'/('{}{}_vs_3D4_{}_IgG_rep1.MACS2_ppois_Sscrofa11.1.bigWig'.format(group,'_rep1-rep2_pooled',state))
            document['merged_tracks'].append({'group':group,'path':str(bw), 'exists':bw.is_file(),
                'bytes':bw.stat().st_size if bw.is_file() else None})

    # The old fragment summary remains separate: it is not asserted to be a new-run measurement.
    old_fragment = project/'logs/fragment_length.tsv'
    if old_fragment.is_file():
        capture(old_fragment)
        document['definitions']['old_fragment_summary'] = 'logs/fragment_length.tsv is a previous-run source, to be reused only after confirming BAM provenance.'
    dest.parent.mkdir(parents=True, exist_ok=True)
    partial = dest.with_name(dest.name+'.partial')
    with gzip.open(partial, 'wt', encoding='utf-8') as handle:
        json.dump(document, handle, ensure_ascii=False, allow_nan=False)
    os.replace(partial, dest)
    print('Saved:', dest)
    print('Peak files:',len(document['peak_files']),'Reproducibility summaries:',len(document['reproducibility']))
    print('Merged bigWigs present:',sum(t['exists'] and bool(t['bytes']) for t in document['merged_tracks']),'/ 4')
    print('Missing:',len(document['missing']),'Errors:',len(document['errors']),'Skipped oversized reports:',len(document['skipped']))
    print('Upload this JSON.GZ file for workbook completion. Missing outputs are not reported as zero.')

if __name__ == '__main__':
    main()
PY
